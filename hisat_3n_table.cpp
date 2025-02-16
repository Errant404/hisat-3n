/*
 * hisat_3n_table.cpp
 *
 * Supported parameters:
 *   --alignments <alignmentFile>   (required) Sorted SAM file (use "-" to read from standard input)
 *   --ref <refFile>                (required) Reference genome FASTA file
 *   --output-name <outputFile>     (optional) Output file name (default output to stdout)
 *   --base-change <char1,char2>    (required) Base-change rule, e.g., "C,T"
 *   -u, --unique-only              Count only uniquely aligned reads
 *   -m, --multiple-only            Count only multiply aligned reads
 *   -c, --CG-only                  Count only CpG sites in the reference (e.g., used for bisulfite sequencing data)
 *       --added-chrname            Chromosome names in SAM have a "chr" prefix (while the reference does not)
 *       --removed-chrname          Chromosome names in SAM do not have the "chr" prefix (while the reference does)
 *   -p, --threads <int>            Number of threads (default: 1)
 *   -h, --help                     Print this help message and exit
 *
 * Output format is TSV with 7 columns (with header):
 *   ref, pos, strand, convertedBaseQualities, convertedBaseCount, unconvertedBaseQualities, unconvertedBaseCount
 *
 * Notes:
 *   1. The program uses a producer–consumer model to stream alignment records from the SAM file.
 *      Worker threads parse each record and directly update the global conversion table (globalConvTable),
 *      while also updating the current maximum reference position (globalMaxPos) for each chromosome.
 *   2. A separate flusher thread periodically checks globalConvTable and, for records on a given chromosome
 *      whose positions are much lower than the chromosome's current maximum position (minus a window value),
 *      assumes that no new data will be added, outputs them to the target file, and then removes them from memory to reduce memory usage.
 */

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

extern "C" {
#include "htslib/faidx.h"
#include "htslib/sam.h"
}

// ----------------------- Configuration Parameters Structure -----------------------
struct Config {
  std::string alignmentsFile;
  std::string refFile;
  std::string outputFile; // If empty, output to stdout
  char from_base;         // e.g., 'C'
  char to_base;           // e.g., 'T'
  bool unique_only = false;
  bool multiple_only = false;
  bool cg_only = false;
  bool added_chrname = false;
  bool removed_chrname = false;
  int threads = 1;
} config;

// ----------------------- Key Structure for Conversion Statistics -----------------------
struct Key {
  std::string chr; // Chromosome name
  int pos;         // 1-based reference position
  char strand;     // '+' or '-'
  bool operator==(const Key &other) const {
    return chr == other.chr && pos == other.pos && strand == other.strand;
  }
};

struct KeyHash {
  std::size_t operator()(const Key &k) const {
    std::size_t h1 = std::hash<std::string>()(k.chr);
    std::size_t h2 = std::hash<int>()(k.pos);
    std::size_t h3 = std::hash<char>()(k.strand);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

struct Value {
  std::string convQual; // Quality string for converted bases
  int convCount = 0;
  std::string unconvQual; // Quality string for unconverted bases
  int unconvCount = 0;
};

using ConvTable = std::unordered_map<Key, Value, KeyHash>;
using Entry = std::pair<Key, Value>;

// ----------------------- Global Statistical Data and Locks -----------------------
std::mutex globalMutex;
ConvTable globalConvTable;
// Records the current maximum reference position (1-based) for each chromosome
std::unordered_map<std::string, int> globalMaxPos;

// Output file stream (used exclusively by the flusher thread)
std::ofstream ofs;

// Define a "flush window": for a given chromosome, if a record's position is much lower than the chromosome's current maximum position (minus the window value),
// it is assumed that no new data will be added for that position.
const int FLUSH_WINDOW = 1000;

// ----------------------- Thread-Safe Queue (TSQueue) Template Class -----------------------
template <typename T> class TSQueue {
public:
  TSQueue(size_t capacity) : cap(capacity) {}
  void push(const T &item) {
    std::unique_lock<std::mutex> lock(mtx);
    cond_full.wait(lock, [&] { return queue.size() < cap; });
    queue.push(item);
    lock.unlock();
    cond_empty.notify_one();
  }
  bool pop(T &item) {
    std::unique_lock<std::mutex> lock(mtx);
    cond_empty.wait(lock, [&] { return !queue.empty(); });
    item = queue.front();
    queue.pop();
    lock.unlock();
    cond_full.notify_one();
    return true;
  }

private:
  std::queue<T> queue;
  std::mutex mtx;
  std::condition_variable cond_empty;
  std::condition_variable cond_full;
  size_t cap;
};

// ----------------------- Helper Function: Compute Complementary Base -----------------------
char complement(char base) {
  base = std::toupper(base);
  switch (base) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  default:
    return base;
  }
}

// ----------------------- Help Information -----------------------
void print_usage(std::ostream &out) {
  out << "Usage: hisat-3n-table [options] --alignments <alignmentFile> --ref <refFile> --base-change <char1,char2>\n"
      << "Options:\n"
      << "  --alignments <alignmentFile>   Sorted SAM file (use '-' for stdin)\n"
      << "  --ref <refFile>                Reference genome (FASTA format)\n"
      << "  --output-name <outputFile>     Output TSV file (default: stdout)\n"
      << "  --base-change <char1,char2>    Base-change rule (e.g. C,T)\n"
      << "  -u, --unique-only              Count only uniquely aligned reads\n"
      << "  -m, --multiple-only            Count only multiply aligned reads\n"
      << "  -c, --CG-only                  Count only CpG sites\n"
      << "      --added-chrname            SAM chromosome names have 'chr' prefix\n"
      << "      --removed-chrname          SAM chromosome names lack 'chr' prefix\n"
      << "  -p, --threads <int>            Number of worker threads (default: 1)\n"
      << "  -h, --help                   Print this help message and exit\n";
}

// ----------------------- Argument Parsing -----------------------
void parse_arguments(int argc, char *argv[], Config &config) {
  static struct option long_options[] = {
      {"alignments", required_argument, 0, 0},
      {"ref", required_argument, 0, 0},
      {"output-name", required_argument, 0, 0},
      {"base-change", required_argument, 0, 0},
      {"unique-only", no_argument, 0, 'u'},
      {"multiple-only", no_argument, 0, 'm'},
      {"CG-only", no_argument, 0, 'c'},
      {"added-chrname", no_argument, 0, 0},
      {"removed-chrname", no_argument, 0, 0},
      {"threads", required_argument, 0, 'p'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};
  int option_index = 0;
  int c;
  while ((c = getopt_long(argc, argv, "ump:hc", long_options, &option_index)) != -1) {
    switch (c) {
    case 0:
      if (strcmp(long_options[option_index].name, "alignments") == 0) {
        config.alignmentsFile = optarg;
      } else if (strcmp(long_options[option_index].name, "ref") == 0) {
        config.refFile = optarg;
      } else if (strcmp(long_options[option_index].name, "output-name") == 0) {
        config.outputFile = optarg;
      } else if (strcmp(long_options[option_index].name, "base-change") == 0) {
        char *comma = strchr(optarg, ',');
        if (!comma || (comma - optarg) != 1 || strlen(comma) < 2) {
          std::cerr << "Invalid --base-change format. Expected format: X,Y" << std::endl;
          exit(EXIT_FAILURE);
        }
        config.from_base = optarg[0];
        config.to_base = comma[1];
      } else if (strcmp(long_options[option_index].name, "added-chrname") == 0) {
        config.added_chrname = true;
      } else if (strcmp(long_options[option_index].name, "removed-chrname") == 0) {
        config.removed_chrname = true;
      }
      break;
    case 'u':
      config.unique_only = true;
      break;
    case 'm':
      config.multiple_only = true;
      break;
    case 'c':
      config.cg_only = true;
      break;
    case 'p':
      config.threads = std::atoi(optarg);
      if (config.threads < 1)
        config.threads = 1;
      break;
    case 'h':
      print_usage(std::cout);
      exit(EXIT_SUCCESS);
      break;
    case '?':
    default:
      print_usage(std::cerr);
      exit(EXIT_FAILURE);
    }
  }
  if (config.alignmentsFile.empty() || config.refFile.empty() ||
      config.from_base == '\0' || config.to_base == '\0') {
    std::cerr << "Missing required parameters." << std::endl;
    print_usage(std::cerr);
    exit(EXIT_FAILURE);
  }
}

// ----------------------- Reader Thread (Producer) -----------------------
// This thread sequentially reads alignment records from the SAM file and pushes each record into the thread-safe queue.
// After reading is complete, it pushes a nullptr for each worker thread as a termination signal.
void reader_thread(samFile *in, bam_hdr_t *hdr, TSQueue<bam1_t *> &queue) {
  bam1_t *aln = bam_init1();
  while (sam_read1(in, hdr, aln) >= 0) {
    bam1_t *dup = bam_dup1(aln);
    queue.push(dup);
  }
  bam_destroy1(aln);
  // Push a termination signal for each worker thread
  for (int i = 0; i < config.threads; i++) {
    queue.push(nullptr);
  }
}

// ----------------------- Worker Thread (Consumer): Process Alignment Records and Update Global Conversion Table -----------------------
void worker_thread(TSQueue<bam1_t *> &queue, const bam_hdr_t *hdr,
                   const std::vector<std::string> &refSeqs,
                   const std::vector<int> &refLens, const Config &config,
                   std::atomic<bool> &workersDone) {
  while (true) {
    bam1_t *aln = nullptr;
    queue.pop(aln);
    if (aln == nullptr) { // Termination signal encountered, exit thread
      break;
    }
    // Determine the uniqueness of the alignment (NH tag, default is unique)
    int nh = 1;
    uint8_t *nh_ptr = bam_aux_get(aln, "NH");
    if (nh_ptr) {
      nh = bam_aux2i(nh_ptr);
    }
    if (config.unique_only && nh != 1) {
      bam_destroy1(aln);
      continue;
    }
    if (config.multiple_only && nh == 1) {
      bam_destroy1(aln);
      continue;
    }

    int tid = aln->core.tid;
    if (tid < 0 || tid >= (int)hdr->n_targets) {
      bam_destroy1(aln);
      continue;
    }
    std::string chr = hdr->target_name[tid]; // Use chromosome name from SAM
    char strand = (aln->core.flag & BAM_FREVERSE) ? '-' : '+';

    char from_base, to_base;
    if (strand == '-') {
      from_base = std::toupper(complement(config.from_base));
      to_base = std::toupper(complement(config.to_base));
    } else {
      from_base = std::toupper(config.from_base);
      to_base = std::toupper(config.to_base);
    }

    // Parse CIGAR string to map alignment to the reference
    uint32_t *cigar = bam_get_cigar(aln);
    int32_t pos = aln->core.pos; // 0-based reference start position
    int readPos = 0;             // Current position in the read
    uint8_t *seq = bam_get_seq(aln);
    uint8_t *qual = bam_get_qual(aln);

    for (int j = 0; j < aln->core.n_cigar; j++) {
      int op = bam_cigar_op(cigar[j]);
      int oplen = bam_cigar_oplen(cigar[j]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        for (int k = 0; k < oplen; k++) {
          int refPos = pos + k; // 0-based reference position
          int currentReadPos = readPos + k;
          // Get the base and quality from the read in the BAM record
          uint8_t base = bam_seqi(seq, currentReadPos);
          char baseChar = std::toupper((int)seq_nt16_str[base]);
          char qchar = static_cast<char>(qual[currentReadPos] + 33);
          if (refPos < 0 || refPos >= refLens[tid])
            continue;
          char refBase = std::toupper(refSeqs[tid][refPos]);

          // If CG-only is required, check the CpG context (determine separately for forward and reverse strands)
          bool isCpG = true;
          if (config.cg_only) {
            if (strand == '+') {
              if (refBase != 'C' || refPos + 1 >= refLens[tid] ||
                  std::toupper(refSeqs[tid][refPos + 1]) != 'G')
                isCpG = false;
            } else { // For '-' strand: check the current base and the previous one
              if (refBase != 'G' || refPos - 1 < 0 ||
                  std::toupper(refSeqs[tid][refPos - 1]) != 'C')
                isCpG = false;
            }
          }
          if (config.cg_only && !isCpG)
            continue;

          // When the reference base equals from_base, determine conversion based on the read base
          if (refBase == from_base) {
            Key key{chr, refPos + 1, strand}; // Convert reference position to 1-based for output
            std::lock_guard<std::mutex> lock(globalMutex);
            // Update the maximum reference position for this chromosome
            if (globalMaxPos.find(chr) == globalMaxPos.end() || globalMaxPos[chr] < refPos + 1)
                globalMaxPos[chr] = refPos + 1;
            // Determine conversion: if the read base equals to_base, it's converted; if it equals from_base, it's unconverted
            if (baseChar == to_base) {
              globalConvTable[key].convQual.push_back(qchar);
              globalConvTable[key].convCount++;
            } else if (baseChar == from_base) {
              globalConvTable[key].unconvQual.push_back(qchar);
              globalConvTable[key].unconvCount++;
            }
          }
        }
        readPos += oplen;
        pos += oplen;
      } else if (op == BAM_CINS) {
        readPos += oplen;
      } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
        pos += oplen;
      } else if (op == BAM_CSOFT_CLIP) {
        readPos += oplen;
      }
      // BAM_CHARD_CLIP is not processed
    }
    bam_destroy1(aln);
  }
  workersDone = true;
}

// ----------------------- Flusher Thread -----------------------
// This thread periodically flushes stable entries from the global conversion table and outputs them.
void flusher_thread(const Config &config, std::atomic<bool> &allWorkersDone,
                    std::atomic<bool> &flusherDone) {
  while (!allWorkersDone) {
    // Flush every 1 second
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::vector<Entry> flushEntries;
    {
      std::lock_guard<std::mutex> lock(globalMutex);
      // Traverse globalConvTable and collect records that meet the condition into flushEntries
      for (auto it = globalConvTable.begin(); it != globalConvTable.end();) {
        const std::string &chr = it->first.chr;
        int pos = it->first.pos;
        if (globalMaxPos.find(chr) != globalMaxPos.end() &&
            pos < globalMaxPos[chr] - FLUSH_WINDOW) {
          flushEntries.push_back(*it);
          it = globalConvTable.erase(it);
        } else {
          ++it;
        }
      }
    }
    // Sort the collected records by (chr, pos, strand)
    std::sort(flushEntries.begin(), flushEntries.end(),
              [](const Entry &a, const Entry &b) {
                if (a.first.chr != b.first.chr)
                  return a.first.chr < b.first.chr;
                if (a.first.pos != b.first.pos)
                  return a.first.pos < b.first.pos;
                return a.first.strand < b.first.strand;
              });
    // Output the sorted records
    for (const auto &entry : flushEntries) {
      if (ofs.is_open())
        ofs << entry.first.chr << "\t" << entry.first.pos << "\t"
            << entry.first.strand << "\t" << entry.second.convQual << "\t"
            << entry.second.convCount << "\t" << entry.second.unconvQual
            << "\t" << entry.second.unconvCount << "\n";
      else
        std::cout << entry.first.chr << "\t" << entry.first.pos << "\t"
                  << entry.first.strand << "\t" << entry.second.convQual
                  << "\t" << entry.second.convCount << "\t"
                  << entry.second.unconvQual << "\t"
                  << entry.second.unconvCount << "\n";
    }
  }
  // After all worker threads have finished, output all remaining entries in globalConvTable
  std::vector<Entry> remainingEntries;
  {
    std::lock_guard<std::mutex> lock(globalMutex);
    for (const auto &entry : globalConvTable)
      remainingEntries.push_back(entry);
    globalConvTable.clear();
  }
  std::sort(remainingEntries.begin(), remainingEntries.end(),
            [](const Entry &a, const Entry &b) {
              if (a.first.chr != b.first.chr)
                return a.first.chr < b.first.chr;
              if (a.first.pos != b.first.pos)
                return a.first.pos < b.first.pos;
              return a.first.strand < b.first.strand;
            });
  
  for (const auto &entry : remainingEntries) {
    if (ofs.is_open())
      ofs << entry.first.chr << "\t" << entry.first.pos << "\t"
          << entry.first.strand << "\t" << entry.second.convQual << "\t"
          << entry.second.convCount << "\t" << entry.second.unconvQual << "\t"
          << entry.second.unconvCount << "\n";
    else
      std::cout << entry.first.chr << "\t" << entry.first.pos << "\t"
                << entry.first.strand << "\t" << entry.second.convQual << "\t"
                << entry.second.convCount << "\t" << entry.second.unconvQual << "\t"
                << entry.second.unconvCount << "\n";
  }
  
  flusherDone = true;
}

// ----------------------- Main Function -----------------------
int main(int argc, char *argv[]) {
  parse_arguments(argc, argv, config);

  // Open the reference genome (FASTA) using HTSlib faidx
  faidx_t *fai = fai_load(config.refFile.c_str());
  if (!fai) {
    std::cerr << "Failed to load reference file: " << config.refFile << std::endl;
    exit(EXIT_FAILURE);
  }

  // Open the SAM file
  samFile *in = sam_open(config.alignmentsFile.c_str(), "r");
  if (!in) {
    std::cerr << "Failed to open alignment file: " << config.alignmentsFile << std::endl;
    exit(EXIT_FAILURE);
  }
  bam_hdr_t *hdr = sam_hdr_read(in);
  if (!hdr) {
    std::cerr << "Failed to read SAM header" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Pre-load the reference sequence for each chromosome
  int n_targets = hdr->n_targets;
  std::vector<std::string> refSeqs(n_targets);
  std::vector<int> refLens(n_targets, 0);
  for (int i = 0; i < n_targets; i++) {
    std::string chrName = hdr->target_name[i];
    std::string fetchName = chrName;
    if (config.added_chrname) {
      if (chrName.compare(0, 3, "chr") == 0)
        fetchName = chrName.substr(3);
    } else if (config.removed_chrname) {
      if (chrName.compare(0, 3, "chr") != 0)
        fetchName = "chr" + chrName;
    }
    int len = 0;
    char *seq = fai_fetch(fai, fetchName.c_str(), &len);
    if (seq) {
      refSeqs[i] = std::string(seq);
      refLens[i] = len;
      free(seq);
    } else {
      std::cerr << "Warning: Failed to fetch reference for " << fetchName << std::endl;
      refSeqs[i] = "";
      refLens[i] = 0;
    }
  }

  // If an output file is specified, open it
  if (!config.outputFile.empty()) {
    ofs.open(config.outputFile);
    if (!ofs) {
      std::cerr << "Cannot open output file: " << config.outputFile << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // Output header
  if (ofs.is_open())
    ofs << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount\n";
  else
    std::cout << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount\n";


  // Construct a thread-safe queue (e.g., with capacity 1000)
  const size_t QUEUE_CAPACITY = 1000;
  TSQueue<bam1_t *> queue(QUEUE_CAPACITY);

  // Start the producer thread to read SAM records and push them into the queue
  std::thread prod(reader_thread, in, hdr, std::ref(queue));

  std::atomic<bool> workersDone(false);
  std::vector<std::thread> workers;
  for (int i = 0; i < config.threads; i++) {
    workers.emplace_back(worker_thread, std::ref(queue), hdr, std::ref(refSeqs),
                         std::ref(refLens), std::ref(config),
                         std::ref(workersDone));
  }

  // Start the flusher thread to periodically output stable entries from the conversion table
  std::atomic<bool> allWorkersDone(false);
  std::atomic<bool> flusherDone(false);
  std::thread flusher(flusher_thread, std::ref(config),
                      std::ref(allWorkersDone), std::ref(flusherDone));

  // Wait for the producer and worker threads to finish
  prod.join();
  for (auto &th : workers) {
    th.join();
  }
  allWorkersDone = true;
  // Wait for the flusher thread to finish flushing
  flusher.join();

  // Close the SAM file and related resources
  sam_close(in);
  bam_hdr_destroy(hdr);
  fai_destroy(fai);
  if (ofs.is_open())
    ofs.close();

  return 0;
}
