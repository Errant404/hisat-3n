/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <atomic>
#include <iostream>
#include <getopt.h>
#include "position_3n_table.h"

using namespace std;

string alignmentFileName;
bool standardInMode = false;
string refFileName;
string outputFileName;
bool uniqueOnly = false;
bool multipleOnly = false;
bool CG_only = false;
int nThreads = 1;
long long int loadingBlockSize = 1000000;
char convertFrom = '0';
char convertTo = '0';
char convertFromComplement;
char convertToComplement;
bool addedChrName = false;
bool removedChrName = false;

bool fileExist (string& filename) {
    ifstream file(filename);
    return file.good();
}

enum {
    ARG_ADDED_CHRNAME = 256,
    ARG_REMOVED_CHRNAME
};

static const char *short_options = "s:r:t:b:umcp:h";
static struct option long_options[] {
                {"alignments",  required_argument, 0, 'a'},
                {"ref",  required_argument, 0, 'r'},
                {"output-name", required_argument, 0, 'o'},
                {"base-change", required_argument, 0, 'b'},
                {"unique-only", no_argument, 0, 'u'},
                {"multiple-only", no_argument, 0, 'm'},
                {"CG-only", no_argument, 0, 'c'},
                {"threads", required_argument, 0, 'p'},
                {"added-chrname", no_argument, 0, ARG_ADDED_CHRNAME },
                {"removed-chrname", no_argument, 0, ARG_REMOVED_CHRNAME },
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
        };

static void printHelp(ostream& out) {
    out << "hisat-3n-table developed by Yun (Leo) Zhang" << endl;
    out << "Usage:" << endl
        << "hisat-3n-table [options]* --alignments <alignmentFile> --ref <refFile> --output-name <outputFile> --base-change <char1,char2>" << endl
        << "  <alignmentFile>           SORTED SAM filename. Please enter '-' for standard input." << endl
        << "  <refFile>                 reference file (should be FASTA format)." << endl
        << "  <outputFile>              file name to save the 3n table (tsv format). By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console)." << endl
        << "  <chr1,chr2>               the char1 is the nucleotide converted from, the char2 is the nucleotide converted to." << endl;
    out << "Options (defaults in parentheses):" << endl
        << " Input:" << endl
        << "  -u/--unique-only          only count the base which is in unique mapped reads." << endl
        << "  -m/--multiple-only        only count the base which is in multiple mapped reads." << endl
        << "  -c/--CG-only              only count CG and ignore CH in reference." << endl
        << "  --added-chrname           please add this option if you use --add-chrname during HISAT-3N alignment." << endl
        << "  --removed-chrname         please add this option if you use --remove-chrname during HISAT-3N alignment." << endl
        << "  -p/--threads <int>        number of threads to launch (1)." << endl
        << "  -h/--help                 print this usage message." << endl;
}

static void parseOption(int next_option, const char *optarg) {
    switch (next_option) {
        case 'a': {
            alignmentFileName = optarg;
            if (alignmentFileName == "-") {
                standardInMode = true;
                break;
            }
            if (!fileExist(alignmentFileName)) {
                cerr << "The alignment file is not exist." << endl;
                throw (1);
            }
            break;
        }
        case 'r': {
            refFileName = optarg;
            if (!fileExist(refFileName)) {
                cerr << "reference (FASTA) file is not exist." << endl;
                throw (1);
            }
            break;
        }
        case 'o':
            outputFileName = optarg;
            break;
        case 'b': {
            string arg = optarg;
            if (arg.size() != 3 || arg[1] != ',') {
                cerr << "Error: expected 2 comma-separated "
                     << "arguments to --base-change option (e.g. C,T), got " << arg << endl;
                throw 1;
            }
            convertFrom = toupper(arg.front());
            convertTo = toupper(arg.back());
            break;
        }
        case 'u':{
            uniqueOnly = true;
            break;
        }
        case 'm': {
            multipleOnly = true;
            break;
        }
        case 'c': {
            CG_only = true;
            break;
        }
        case 'h': {
            printHelp(cerr);
            throw 0;
        }
        case 'p': {
            nThreads = stoi(optarg);
            if (nThreads < 1) {
                nThreads = 1;
            }
            break;
        }
        case ARG_ADDED_CHRNAME: {
            addedChrName = true;
            break;
        }
        case ARG_REMOVED_CHRNAME: {
            removedChrName = true;
            break;
        }
        default:
            printHelp(cerr);
            throw 1;
    }
}

static void parseOptions(int argc, const char **argv) {
    int option_index = 0;
    int next_option;
    while (true) {
        next_option = getopt_long(argc, const_cast<char **>(argv), short_options,
                                  long_options, &option_index);
        if (next_option == -1)
            break;
        parseOption(next_option, optarg);
    }

    // check filenames
    if (refFileName.empty() || alignmentFileName.empty()) {
        cerr << "No reference or SAM file specified!" << endl;
        printHelp(cerr);
        throw 1;
    }

    // give a warning for CG-only
    if (CG_only) {
        if (convertFrom != 'C' || convertTo != 'T') {
            cerr << "Warning! You are using CG-only mode. The the --base-change option is set to: C,T" << endl;
            convertFrom = 'C';
            convertTo = 'T';
        }
    }

    // check if --base-change is empty
    if (convertFrom == '0' || convertTo == '0') {
        cerr << "the --base-change argument is required." << endl;
        throw 1;
    }

    if(removedChrName && addedChrName) {
        cerr << "Error: --removed-chrname and --added-chrname cannot be used at the same time" << endl;
        throw 1;
    }

    // set complements
    convertFromComplement = asc2dnacomp[convertFrom];
    convertToComplement = asc2dnacomp[convertTo];
}

/**
 * give a SAM line, extract the chromosome and position information.
 * return true if the SAM line is mapped. return false if SAM line is not maped.
 */
bool getSAMChromosomePos(string* line, string& chr, long long int& pos) {
    int startPosition = 0;
    int endPosition = 0;
    int count = 0;

    while ((endPosition = line->find("\t", startPosition)) != string::npos) {
        if (count == 2) {
            chr = line->substr(startPosition, endPosition - startPosition);
        } else if (count == 3) {
            pos = stoll(line->substr(startPosition, endPosition - startPosition));
            if (chr == "*") {
                return false;
            } else {
                return true;
            }
        }
        startPosition = endPosition + 1;
        count++;
    }
    return false;
}

/**
 * store all reference position in this class.
 */
class Positions{
public:
    vector<Position*> refPositions; // the pool of all current reference position.
    string chromosome; // current reference chromosome name.
    long long int location; // current location (position) in reference chromosome.
    char lastBase = 'X'; // the last base of reference line. this is for CG_only mode.
    SafeQueue<string*> linePool; // pool to store unprocessed SAM line.
    SafeQueue<string*> freeLinePool; // pool to store free string pointer for SAM line.
    SafeQueue<Position*> freePositionPool; // pool to store free position pointer for reference position.
    SafeQueue<Position*> outputPositionPool; // pool to store the reference position which is loaded and ready to output.
    bool working;
    long long int refCoveredPosition; // this is the last position in reference chromosome we loaded in refPositions.
    ifstream refFile;
    int nThreads = 1;
    ChromosomeFilePositions chromosomePos; // store the chromosome name and it's streamPos. To quickly find new chromosome in file.
    bool addedChrName = false;
    bool removedChrName = false;
    atomic_size_t numTasks{0};

    Positions(string inputRefFileName, int inputNThreads, bool inputAddedChrName, bool inputRemovedChrName) {
        working = true;
        nThreads = inputNThreads;
        addedChrName = inputAddedChrName;
        removedChrName = inputRemovedChrName;
        refFile.open(inputRefFileName, ios_base::in);
        LoadChromosomeNamesPos();
    }

    ~Positions() {
        Position* pos;
        while(freePositionPool.popFront(pos)) {
            delete pos;
        }
    }

    /**
     * given the target Position output the corresponding position index in refPositions.
     */
    int getIndex(long long int &targetPos) {
        int firstPos = refPositions[0]->location;
        return targetPos - firstPos;
    }

    /**
     * given reference line (start with '>'), extract the chromosome information.
     * this is important when there is space in chromosome name. the SAM information only contain the first word.
     */
    string getChrName(string& inputLine) {
        string name;
        for (int i = 1; i < inputLine.size(); i++)
        {
            char c = inputLine[i];
            if (isspace(c)){
                break;
            }
            name += c;
        }

        if(removedChrName) {
            if(name.find("chr") == 0) {
                name = name.substr(3);
            }
        } else if(addedChrName) {
            if(name.find("chr") != 0) {
                name = string("chr") + name;
            }
        }
        return name;
    }


    /**
     * Scan the reference file. Record each chromosome and its position in file.
     */
    void LoadChromosomeNamesPos() {
        string line;
        while (getline(refFile, line)) {
            if (line.front() == '>') { // this line is chromosome name
                chromosome = getChrName(line);
                streampos currentPos = refFile.tellg();
                chromosomePos.append(chromosome, currentPos);
            }
        }
        chromosomePos.sort();
        chromosome.clear();
    }

    /**
     * get a fasta line (not header), append the bases to positions.
     */
    void appendRefPosition(string& line) {
        Position* newPos;
        // check the base one by one
        char* b;
        for (int i = 0; i < line.size(); i++) {
            getFreePosition(newPos);
            newPos->set(chromosome, location+i);
            b = &line[i];
            if (CG_only) {
                if (lastBase == 'C' && *b == 'G') {
                    refPositions.back()->set('+');
                    newPos->set('-');
                }
            } else {
                if (*b == convertFrom) {
                    newPos->set('+');
                } else if (*b == convertFromComplement) {
                    newPos->set('-');
                }
            }
            refPositions.push_back(newPos);
            lastBase = *b;
        }
        location += line.size();
    }

    /**
     * if we can go through all the workerLock, that means no worker is appending new position.
     */
    void appendingFinished() {
        while (numTasks) {}
    }

    /**
     * the output function for output thread.
     */
    void outputFunction(string outputFileName) {
        ostream* out_ = &cout;
        out_ = &cout;
        ofstream tableFile;
        if (!outputFileName.empty()) {
            tableFile.open(outputFileName, ios_base::out);
            out_ = &tableFile;
        }

        *out_ << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount\n";
        Position* pos;
        while (true) {
            // if (outputPositionPool.popFront(pos)) 
            {
                unique_lock<std::mutex> lock(outputPositionPool.mutex_);
                outputPositionPool.cv.wait(lock, [this] { return !outputPositionPool.queue_.empty() || !working; });
                if (!working) { break;}
                pos = outputPositionPool.queue_.front();
                outputPositionPool.queue_.pop();
                lock.unlock();
                outputPositionPool.cv.notify_all();
            }
            *out_ << pos->chromosome << '\t'
                        << to_string(pos->location) << '\t'
                        << pos->strand << '\t'
                        << pos->convertedQualities << '\t'
                        << to_string(pos->convertedQualities.size()) << '\t'
                        << pos->unconvertedQualities << '\t'
                        << to_string(pos->unconvertedQualities.size()) << '\n';
            returnPosition(pos);
        }
            // this_thread::sleep_for (std::chrono::microseconds(1));
        tableFile.close();
    }

    /**
     * move the position which position smaller than refCoveredPosition - loadingBlockSize, output it.
     */
    void moveBlockToOutput() {
        if (refPositions.empty()) {
            return;
        }
        int index;
        for (index = 0; index < refPositions.size(); index++) {
            if (refPositions[index]->location < refCoveredPosition - loadingBlockSize) {
                if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                    returnPosition(refPositions[index]);
                } else {
                    outputPositionPool.push(refPositions[index]);
                }
            } else {
                break;
            }
        }
        if (index != 0) {
            refPositions.erase(refPositions.begin(), refPositions.begin()+index);
        }
        outputPositionPool.cv.notify_all();
    }

    /**
     * move all the refPosition into output pool.
     */
    void moveAllToOutput() {
        if (refPositions.empty()) {
            return;
        }
        for (int index = 0; index < refPositions.size(); index++) {
            if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                returnPosition(refPositions[index]);
            } else {
                vector<uniqueID>().swap(refPositions[index]->uniqueIDs);
                outputPositionPool.push(refPositions[index]);
            }
        }
        refPositions.clear();
        outputPositionPool.cv.notify_all();
    }

    /**
     * initially load reference sequence for 2 million bp
     */
    void loadNewChromosome(string targetChromosome) {
        refFile.clear();
        // find the start position in file based on chromosome name.
        streampos startPos = chromosomePos.getChromosomePosInRefFile(targetChromosome);
        chromosome = targetChromosome;
        refFile.seekg(startPos, ios::beg);
        refCoveredPosition = 2 * loadingBlockSize;
        string line;
        lastBase = 'X';
        location = 0;
        while (getline(refFile, line)) {
            if (line.front() == '>') { // this line is chromosome name
                return; // meet next chromosome, return it.
            } else {
                if (line.empty()) { continue; }
                // change all base to upper case
                for (int i = 0; i < line.size(); i++) {
                    line[i] = toupper(line[i]);
                }
                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return;
                }
            }
        }
        freePositionPool.cv.notify_all();
    }

    /**
     * load more Position (loadingBlockSize bp) to positions
     * if we meet next chromosome, return false. Else, return ture.
     */
    void loadMore() {
        refCoveredPosition += loadingBlockSize;
        string line;
        while (getline(refFile, line)) {
            if (line.front() == '>') { // meet next chromosome, return.
                return ;
            } else {
                if (line.empty()) { continue; }

                // change all base to upper case
                for (int i = 0; i < line.size(); i++) {
                    line[i] = toupper(line[i]);
                }

                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return ;
                }
            }
        }
        freePositionPool.cv.notify_all();
    }


    /**
     * add position information from Alignment into ref position.
     */
    void appendPositions(Alignment& newAlignment) {
        if (!newAlignment.mapped || newAlignment.bases.empty()) {
            return;
        }
        long long int startPos = newAlignment.location; // 1-based position
        // find the first reference position in pool.
        int index = getIndex(newAlignment.location);

        for (int i = 0; i < newAlignment.sequence.size(); i++) {
            PosQuality* b = &newAlignment.bases[i];
            if (b->remove) {
                continue;
            }

            Position* pos = refPositions[index+b->refPos];
            assert (pos->location == startPos + b->refPos);

            if (pos->strand == '?') {
                // this is for CG-only mode. read has a 'C' or 'G' but not 'CG'.
                continue;
            }
            pos->appendBase(newAlignment.bases[i], newAlignment);
        }
    }

    /**
     * get a string pointer from freeLinePool, if freeLinePool is empty, make a new string pointer.
     */
    void getFreeStringPointer(string*& newLine) {
        if (freeLinePool.popFront(newLine)) {
            return;
        } else {
            newLine = new string();
        }
    }

    /**
     * get a Position pointer from freePositionPool, if freePositionPool is empty, make a new Position pointer.
     */
    void getFreePosition(Position*& newPosition) {
        // while (outputPositionPool.size() >= 10000) {
        //     this_thread::sleep_for (std::chrono::microseconds(1));
        // }
        {
            unique_lock<std::mutex> lock(outputPositionPool.mutex_);
            outputPositionPool.cv.wait(lock, [this] { return outputPositionPool.queue_.size() < 10000; });
        }
        if (freePositionPool.popFront(newPosition)) {
            return;
        } else {
            newPosition = new Position();
        }
    }

    /**
     * return the line to freeLinePool
     */
    void returnLine(string* line) {
        line->clear();
        freeLinePool.push(line);
    }

    /**
     * return the position to freePositionPool.
     */
    void returnPosition(Position* pos) {
        pos->initialize();
        freePositionPool.push(pos);
    }

    /**
     * this is the working function.
     * it take the SAM line from linePool, parse it.
     */
    void append(int threadID) {
        string* line;
        Alignment newAlignment;

        while (true) {
            // workerLock[threadID]->lock();
            // if(!linePool.popFront(line)) {
            //     workerLock[threadID]->unlock();
            //     this_thread::sleep_for (std::chrono::nanoseconds(1));
            //     continue;
            // }
            {
                unique_lock<std::mutex> lock(linePool.mutex_);
                linePool.cv.wait(lock, [this] { return !linePool.queue_.empty() || !working; });
                if (!working) { break;}
                line = linePool.queue_.front();
                linePool.queue_.pop();
                lock.unlock();
                linePool.cv.notify_one();
            }
            // while (refPositions.empty()) {
            //     this_thread::sleep_for (std::chrono::microseconds(1));
            // }
            {
                unique_lock<std::mutex> lock(freePositionPool.mutex_);
                freePositionPool.cv.wait(lock, [this] { return !refPositions.empty(); });
            }
            newAlignment.parse(line);
            returnLine(line);
            appendPositions(newAlignment);
            numTasks--;
        }
    }
};

// main function, initially 2 load loadingBlockSize (2,000,000) bp of reference, set reloadPos to 1 loadingBlockSize, then load SAM data.
// when the samPos larger than the reloadPos load 1 loadingBlockSize bp of reference.
// when the samChromosome is different to current chromosome, finish all sam position and output all.
int hisat_3n_table()
{
    Positions* positions;

    positions = new Positions(refFileName, nThreads, addedChrName, removedChrName);

    // open #nThreads workers
    vector<thread*> workers;
    for (int i = 0; i < nThreads; i++) {
        workers.push_back(new thread(&Positions::append, positions, i));
    }

    // open a output thread
    thread outputThread;
    outputThread = thread(&Positions::outputFunction, positions, outputFileName);

    ifstream inputFile;
    istream *alignmentFile = &cin;
    if (!standardInMode) {
        inputFile.open(alignmentFileName, ios_base::in);
        alignmentFile = &inputFile;
    }

    string* line; // temporary string to get SAM line.
    string samChromosome; // the chromosome name of current SAM line.
    long long int samPos; // the position of current SAM line.
    long long int reloadPos; // the position in reference that we need to reload.
    long long int lastPos = 0; // the position on last SAM line. compare lastPos with samPos to make sure the SAM is sorted.

    while (true) {
        positions->getFreeStringPointer(line);
        if (!getline(*alignmentFile, *line)) {
            positions->returnLine(line);
            break;
        }

        if (line->empty() || line->front() == '@') {
            positions->returnLine(line);
            continue;
        }
        // limit the linePool size to save memory
        // while(positions->linePool.size() > 1000 * nThreads) {
        //     this_thread::sleep_for (std::chrono::microseconds(1));
        // }
        {
            unique_lock<std::mutex> lock(positions->linePool.mutex_);
            positions->linePool.cv.wait(lock, [positions] { return positions->linePool.queue_.size() < 1000 * nThreads; });
        }
        // if the SAM line is empty or unmapped, get the next SAM line.
        if (!getSAMChromosomePos(line, samChromosome, samPos)) {
            positions->returnLine(line);
            continue;
        }
        // if the samChromosome is different than current positions' chromosome, finish all SAM line.
        // then load a new reference chromosome.
        if (samChromosome != positions->chromosome) {
            // wait all line is processed
            // while (!positions->linePool.empty() || positions->outputPositionPool.size() > 100000) {
            //     this_thread::sleep_for (std::chrono::microseconds(1));
            // }
            {
                unique_lock<std::mutex> lock(positions->linePool.mutex_);
                positions->linePool.cv.wait(lock, [positions] { return positions->linePool.queue_.empty(); });
            }
            {
                unique_lock<std::mutex> lock(positions->outputPositionPool.mutex_);
                positions->outputPositionPool.cv.wait(lock, [positions] { return positions->outputPositionPool.queue_.size() < 100000; });
            }
            positions->appendingFinished();
            positions->moveAllToOutput();
            positions->loadNewChromosome(samChromosome);
            reloadPos = loadingBlockSize;
            lastPos = 0;
        }
        // if the samPos is larger than reloadPos, load 1 loadingBlockSize bp in from reference.
        while (samPos > reloadPos) {
            // while (!positions->linePool.empty() || positions->outputPositionPool.size() > 100000) {
            //     this_thread::sleep_for (std::chrono::microseconds(1));
            // }
            {
                unique_lock<std::mutex> lock(positions->linePool.mutex_);
                positions->linePool.cv.wait(lock, [positions] { return positions->linePool.queue_.empty(); });
            }
            {
                unique_lock<std::mutex> lock(positions->outputPositionPool.mutex_);
                positions->outputPositionPool.cv.wait(lock, [positions] { return positions->outputPositionPool.queue_.size() < 100000; });
            }
            positions->appendingFinished();
            positions->moveBlockToOutput();
            positions->loadMore();
            reloadPos += loadingBlockSize;
        }
        if (lastPos > samPos) {
            cerr << "The input alignment file is not sorted. Please use sorted SAM file as alignment file." << endl;
            throw 1;
        }
        positions->linePool.push(line);
        positions->linePool.cv.notify_one();
        positions->numTasks++;
        lastPos = samPos;
    }
    if (!standardInMode) {
        inputFile.close();
    }

    // prepare to close everything.

    // make sure linePool is empty
    // while (!positions->linePool.empty()) {
    //     this_thread::sleep_for (std::chrono::microseconds(100));
    // }
    {
        unique_lock<std::mutex> lock(positions->linePool.mutex_);
        positions->linePool.cv.wait(lock, [positions] { return positions->linePool.queue_.empty(); });
    }
    // make sure all workers finished their appending work.
    positions->appendingFinished();
    // move all position to outputPool
    positions->moveAllToOutput();
    // wait until outputPool is empty
    // while (!positions->outputPositionPool.empty()) {
    //     this_thread::sleep_for (std::chrono::microseconds(100));
    // }
    {
        unique_lock<std::mutex> lock(positions->outputPositionPool.mutex_);
        positions->outputPositionPool.cv.wait(lock, [positions] { return positions->outputPositionPool.queue_.empty(); });
    }
    // stop all thread and clean
    while(positions->freeLinePool.popFront(line)) {
        delete line;
    }
    positions->working = false;
    positions->linePool.cv.notify_all();
    positions->outputPositionPool.cv.notify_all();
    for (int i = 0; i < nThreads; i++){
        workers[i]->join();
        delete workers[i];
    }
    outputThread.join();
    delete positions;
    return 0;
}

int main(int argc, const char** argv)
{
    int ret = 0;

    try {
        parseOptions(argc, argv);
        ret = hisat_3n_table();
    } catch(std::exception& e) {
        cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
        cerr << "Command: ";
        for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
        cerr << endl;
        return 1;
    } catch(int e) {
        if (e != 0) {
            cerr << "Error: Encountered internal HISAT-3N exception (#" << e << ")" << endl;
            cerr << "Command: ";
            for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
            cerr << endl;
        }
        return e;
    }

    return ret;
}
