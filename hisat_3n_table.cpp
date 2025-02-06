/*
 * hisat_3n_table.cpp -- 流式生成 3N-conversion-table
 *
 * 支持参数：
 *   --alignments <alignmentFile>   （必选）排序好的 SAM 文件（若为 "-"
 * 则从标准输入读取）
 *   --ref <refFile>                （必选）参考基因组 FASTA 文件
 *   --output-name <outputFile>     （可选）输出文件名（默认输出到 stdout）
 *   --base-change <char1,char2>    （必选）碱基转换规则，如 "C,T"
 *   -u, --unique-only              仅统计唯一比对的 reads
 *   -m, --multiple-only            仅统计多重比对的 reads
 *   -c, --CG-only                  仅统计参考中 CpG
 * 位点（例如用于双硫酸盐测序数据）
 *       --added-chrname            SAM 中染色体名称有 "chr"
 * 前缀（参考中无此前缀）
 *       --removed-chrname          SAM 中染色体名称无 "chr"
 * 前缀（参考中有此前缀） -p, --threads <int>            使用的线程数（默认 1）
 *   -h, --help                   打印帮助信息并退出
 *
 * 输出格式为 TSV，共 7 列（带表头）：
 *   ref, pos, strand, convertedBaseQualities, convertedBaseCount,
 * unconvertedBaseQualities, unconvertedBaseCount
 *
 * 说明：
 *   1. 程序采用一个读取线程将比对记录流式推送到容量有限的队列中，
 *      多个工作线程并发消费队列中的记录，处理后更新各自的局部转换统计表，
 *      最后在主线程中合并各线程统计结果并输出。
 *   2. 这样设计可以避免将所有比对记录一次性加载到内存中，适用于大数据量的情况。
 */

#include <algorithm>
#include <cctype>
#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

extern "C" {
#include "htslib/faidx.h"
#include "htslib/sam.h"
}

// ----------------------- 辅助函数 -----------------------
// 返回碱基互补
char complement(char base) {
  switch (std::toupper(base)) {
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

// ----------------------- 程序配置结构 -----------------------
struct Config {
  std::string alignmentsFile;
  std::string refFile;
  std::string outputFile; // 若为空则输出到 stdout
  char from_base;         // 如 'C'
  char to_base;           // 如 'T'
  bool unique_only = false;
  bool multiple_only = false;
  bool cg_only = false;
  bool added_chrname = false;
  bool removed_chrname = false;
  int threads = 1;
} config;

// ----------------------- 转换统计的键值结构 -----------------------
struct Key {
  std::string chr; // 染色体名称
  int pos;         // 1-based参考位置
  char strand;     // '+' 或 '-'
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
  std::string convQual; // 转换碱基的质量字符串
  int convCount = 0;
  std::string unconvQual; // 未转换碱基的质量字符串
  int unconvCount = 0;
};

using ConvTable = std::unordered_map<Key, Value, KeyHash>;

// ----------------------- 线程安全队列 -----------------------
template <typename T> class TSQueue {
public:
  TSQueue(size_t capacity) : cap(capacity) {}
  // push 新数据，如果队列满则等待
  void push(const T &item) {
    std::unique_lock<std::mutex> lock(mtx);
    cond_full.wait(lock, [&] { return queue.size() < cap; });
    queue.push(item);
    lock.unlock();
    cond_empty.notify_one();
  }
  // pop 数据，如果队列为空则等待
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

// ----------------------- 帮助信息 -----------------------
void print_usage(std::ostream &out) {
  out << "Usage: hisat-3n-table [options] --alignments <alignmentFile> --ref "
         "<refFile> --base-change <char1,char2>\n"
      << "Options:\n"
      << "  --alignments <alignmentFile>   Sorted SAM file (use '-' for "
         "stdin)\n"
      << "  --ref <refFile>                Reference genome (FASTA format)\n"
      << "  --output-name <outputFile>     Output TSV file (default: stdout)\n"
      << "  --base-change <char1,char2>    Base-change rule (e.g. C,T)\n"
      << "  -u, --unique-only              Count only uniquely aligned reads\n"
      << "  -m, --multiple-only            Count only multiply aligned reads\n"
      << "  -c, --CG-only                  Count only CpG sites\n"
      << "      --added-chrname            SAM has added 'chr' prefix\n"
      << "      --removed-chrname          SAM has removed 'chr' prefix\n"
      << "  -p, --threads <int>            Number of threads (default: 1)\n"
      << "  -h, --help                   Print this help message and exit\n";
}

// ----------------------- 参数解析 -----------------------
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
  while ((c = getopt_long(argc, argv, "ump:hc", long_options, &option_index)) !=
         -1) {
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
          std::cerr << "Invalid --base-change format. Expected format: X,Y"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
        config.from_base = optarg[0];
        config.to_base = comma[1];
      } else if (strcmp(long_options[option_index].name, "added-chrname") ==
                 0) {
        config.added_chrname = true;
      } else if (strcmp(long_options[option_index].name, "removed-chrname") ==
                 0) {
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

// ----------------------- 读取线程（生产者） -----------------------
// 该线程从 SAM 文件中顺序读取比对记录，并将每个记录 push 到线程安全队列中。
// 读取结束后，为每个工作线程压入一个 nullptr 作为终止标志。
void reader_thread(samFile *in, bam_hdr_t *hdr, TSQueue<bam1_t *> &queue) {
  bam1_t *aln = bam_init1();
  while (sam_read1(in, hdr, aln) >= 0) {
    bam1_t *dup = bam_dup1(aln);
    queue.push(dup);
  }
  bam_destroy1(aln);
  // 读取完毕后，向队列中压入与工作线程数相同的终止标志（nullptr）
  // 注意：这里假设工作线程数量不会变化
  // 为保证所有消费者退出，必须推送足够的终止标志。
  // 此处在 main 中创建队列时，队列容量有限，不会无限增长。
  for(int numWorkers = config.threads; numWorkers > 0; numWorkers--) {
    queue.push(nullptr);
  }
}

// ----------------------- 工作线程函数（消费者） -----------------------
// 每个工作线程不断从队列中取出比对记录，处理后更新本地转换统计表，并释放比对记录内存。
void worker_thread(TSQueue<bam1_t *> &queue, const bam_hdr_t *hdr,
                   const std::vector<std::string> &refSeqs,
                   const std::vector<int> &refLens, const Config &config,
                   ConvTable &localTable) {
  while (true) {
    bam1_t *aln = nullptr;
    queue.pop(aln);
    if (aln == nullptr) {
      // 遇到终止标志则退出
      break;
    }
    // 判断比对的唯一性（NH 标签，默认唯一）
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
    std::string chr = hdr->target_name[tid]; // 使用 SAM 中的染色体名称
    char strand = (aln->core.flag & BAM_FREVERSE) ? '-' : '+';

    char from_base, to_base;
    if (strand == '-') {
      from_base = std::toupper(complement(config.from_base));
      to_base = std::toupper(complement(config.to_base));
    } else {
      from_base = std::toupper(config.from_base);
      to_base = std::toupper(config.to_base);
    }

    // 解析 CIGAR，将比对映射到参考上
    uint32_t *cigar = bam_get_cigar(aln);
    int32_t pos = aln->core.pos; // 0-based参考起始位置
    int readPos = 0;             // read 中当前位置
    uint8_t *seq = bam_get_seq(aln);
    uint8_t *qual = bam_get_qual(aln);

    for (int j = 0; j < aln->core.n_cigar; j++) {
      int op = bam_cigar_op(cigar[j]);
      int oplen = bam_cigar_oplen(cigar[j]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        for (int k = 0; k < oplen; k++) {
          int refPos = pos + k; // 0-based参考位置
          int currentReadPos = readPos + k;
          uint8_t base = bam_seqi(seq, currentReadPos);
          char baseChar = std::toupper((int)seq_nt16_str[base]);
          char qchar = static_cast<char>(qual[currentReadPos] + 33);
          if (refPos < 0 || refPos >= refLens[tid])
            continue;
          char refBase = std::toupper(refSeqs[tid][refPos]);

          // CG-only过滤：检查 CpG 环境
          bool isCpG = true;
          if (config.cg_only) {
            if (strand == '+') {
              if (refBase != 'C' || refPos + 1 >= refLens[tid] ||
                  std::toupper(refSeqs[tid][refPos + 1]) != 'G')
                isCpG = false;
            } else { // '-' 链
              if (refBase != 'G' || refPos - 1 < 0 ||
                  std::toupper(refSeqs[tid][refPos - 1]) != 'C')
                isCpG = false;
            }
          }
          if (config.cg_only && !isCpG)
            continue;

          // 当参考碱基等于 from_base 时，根据 read 上的碱基判断转换情况
          if (refBase == from_base) {
            Key key{chr, refPos + 1, strand}; // 输出 1-based位置
            if (baseChar == to_base) {
              localTable[key].convQual.push_back(qchar);
              localTable[key].convCount++;
            } else if (baseChar == from_base) {
              localTable[key].unconvQual.push_back(qchar);
              localTable[key].unconvCount++;
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
      // BAM_CHARD_CLIP 不处理
    }
    bam_destroy1(aln);
  }
}

// ----------------------- 主函数 -----------------------
int main(int argc, char *argv[]) {
  parse_arguments(argc, argv, config);

  // 加载参考基因组（FASTA）——使用 HTSlib faidx
  faidx_t *fai = fai_load(config.refFile.c_str());
  if (!fai) {
    std::cerr << "Failed to load reference file: " << config.refFile
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // 打开 SAM 文件
  samFile *in = sam_open(config.alignmentsFile.c_str(), "r");
  if (!in) {
    std::cerr << "Failed to open alignment file: " << config.alignmentsFile
              << std::endl;
    exit(EXIT_FAILURE);
  }
  bam_hdr_t *hdr = sam_hdr_read(in);
  if (!hdr) {
    std::cerr << "Failed to read SAM header" << std::endl;
    exit(EXIT_FAILURE);
  }

  // 预先加载每个染色体的参考序列
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
      std::cerr << "Warning: Failed to fetch reference for " << fetchName
                << std::endl;
      refSeqs[i] = "";
      refLens[i] = 0;
    }
  }

  // 构造线程安全队列（设定一个合适的容量，例如 1000 条记录）
  const size_t QUEUE_CAPACITY = 1000;
  TSQueue<bam1_t *> queue(QUEUE_CAPACITY);

  // 启动生产者线程，负责从 SAM 文件中读取比对记录
  std::thread prod([&]() {
    bam1_t *aln = bam_init1();
    while (sam_read1(in, hdr, aln) >= 0) {
      bam1_t *dup = bam_dup1(aln);
      queue.push(dup);
    }
    bam_destroy1(aln);
    // 读取完毕后，向队列中推送终止标志，每个工作线程需要收到一个 nullptr
    for (int i = 0; i < config.threads; i++) {
      queue.push(nullptr);
    }
  });

  // 启动多个工作线程（消费者）
  std::vector<ConvTable> localTables(config.threads);
  std::vector<std::thread> workers;
  for (int i = 0; i < config.threads; i++) {
    workers.emplace_back(worker_thread, std::ref(queue), hdr, std::ref(refSeqs),
                         std::ref(refLens), std::ref(config),
                         std::ref(localTables[i]));
  }

  // 等待生产者线程结束，并关闭 SAM 文件（后续只使用已提交到队列的记录）
  prod.join();
  sam_close(in);

  // 等待所有工作线程完成
  for (auto &th : workers) {
    th.join();
  }

  // 合并各线程局部统计表到全局统计表
  ConvTable globalTable;
  for (const auto &localMap : localTables) {
    for (const auto &pair : localMap) {
      const Key &key = pair.first;
      const Value &val = pair.second;
      Value &gVal = globalTable[key];
      gVal.convQual.append(val.convQual);
      gVal.convCount += val.convCount;
      gVal.unconvQual.append(val.unconvQual);
      gVal.unconvCount += val.unconvCount;
    }
  }

  // 将 globalTable 转换为 vector 并排序输出
  std::vector<std::pair<Key, Value>> entries(globalTable.begin(),
                                             globalTable.end());
  std::sort(entries.begin(), entries.end(),
            [](const std::pair<Key, Value> &a, const std::pair<Key, Value> &b) {
              if (a.first.chr != b.first.chr)
                return a.first.chr < b.first.chr;
              if (a.first.pos != b.first.pos)
                return a.first.pos < b.first.pos;
              return a.first.strand < b.first.strand;
            });

  // 输出结果
  std::ostream *outStream = &std::cout;
  std::ofstream ofs;
  if (!config.outputFile.empty()) {
    ofs.open(config.outputFile);
    if (!ofs) {
      std::cerr << "Cannot open output file: " << config.outputFile
                << std::endl;
      exit(EXIT_FAILURE);
    }
    outStream = &ofs;
  }
  (*outStream) << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount"
                  "\tunconvertedBaseQualities\tunconvertedBaseCount\n";
  for (const auto &pair : entries) {
    const Key &key = pair.first;
    const Value &val = pair.second;
    (*outStream) << key.chr << "\t" << key.pos << "\t" << key.strand << "\t"
                 << val.convQual << "\t" << val.convCount << "\t"
                 << val.unconvQual << "\t" << val.unconvCount << "\n";
  }
  if (ofs.is_open())
    ofs.close();

  // 清理：释放 header 和参考索引
  bam_hdr_destroy(hdr);
  fai_destroy(fai);

  return 0;
}
