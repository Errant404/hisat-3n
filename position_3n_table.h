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

#ifndef POSITION_3N_TABLE_H
#define POSITION_3N_TABLE_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <thread>
#include <cassert>
#include <tbb/concurrent_hash_map.h>
#include <tbb/global_control.h>
#include <tbb/task_group.h>
#include "alignment_3n_table.h"

using namespace std;

extern bool CG_only;
extern long long int loadingBlockSize;
extern tbb::task_group tg;

/**
 * store unique information for one base information with readID, and the quality.
 */
class uniqueID
{
public:
    unsigned long long readNameID;
    bool isConverted;
    char quality;
    bool removed;

    uniqueID() : readNameID(0), isConverted(false), quality(0), removed(false) {}

    uniqueID(unsigned long long InReadNameID,
             bool InIsConverted,
             char& InQual){
        readNameID = InReadNameID;
        isConverted = InIsConverted;
        quality = InQual;
        removed = false;
    }
};

/**
 * basic class to store reference position information
 */
class Position{
private:
    mutex converted_mutex_;
    mutex unconverted_mutex_;

    struct HashCompare {
        static size_t hash(const unsigned long long& key) {
            return key;
        }
        static bool equal(const unsigned long long& key1, const unsigned long long& key2) {
            return key1 == key2;
        }
    };

public:
    string chromosome;      // reference chromosome name
    long long int location; // 1-based position
    char strand;           // +(REF) or -(REF-RC)
    string convertedQualities;    // each char is a mapping quality on this position for converted base.
    string unconvertedQualities;  // each char is a mapping quality on this position for unconverted base.
    
    typedef tbb::concurrent_hash_map<unsigned long long, uniqueID, HashCompare> IdMap;
    IdMap uniqueIDs; // thread-safe hash map

    void initialize() {
        chromosome.clear();
        location = -1;
        strand = '?';
        convertedQualities.clear();
        unconvertedQualities.clear();
        uniqueIDs.clear();
    }

    Position(){
        initialize();
    };

    /**
     * return true if there is mapping information in this reference position.
     */
    bool empty() {
        return convertedQualities.empty() && unconvertedQualities.empty();
    }

    /**
     * set the chromosome, location (position), and strand information.
     */

    void set (string& inputChr, long long int inputLoc) {
        chromosome = inputChr;
        location = inputLoc + 1;
    }

    void set(char inputStrand) {
        strand = inputStrand;
    }

    /**
     * 使用并发哈希表和细粒度锁重新实现appendReadNameID
     * 返回true表示新增了readNameID,false表示已存在或被移除
     */
    bool appendReadNameID(PosQuality& InBase, Alignment& InAlignment) {
        IdMap::accessor acc;
        
        // 尝试插入,如果key不存在则插入并返回true
        if(uniqueIDs.insert(acc, InAlignment.readNameID)) {
            acc->second = uniqueID(InAlignment.readNameID, InBase.converted, InBase.qual);
            return true;
        }
        
        // key已存在,获取现有value
        auto& existingID = acc->second;
        
        // 如果该ID已被移除
        if(existingID.removed) {
            return false; 
        }

        // 如果转换状态不一致,则标记为移除并从相应的qualities中删除
        if(existingID.isConverted != InBase.converted) {
            existingID.removed = true;
            
            if(existingID.isConverted) {
                lock_guard<mutex> lock(converted_mutex_);
                for(int i = 0; i < convertedQualities.size(); i++) {
                    if(convertedQualities[i] == InBase.qual) {
                        convertedQualities.erase(convertedQualities.begin()+i);
                        return false;
                    }
                }
            } else {
                lock_guard<mutex> lock(unconverted_mutex_);
                for(int i = 0; i < unconvertedQualities.size(); i++) {
                    if(unconvertedQualities[i] == InBase.qual) {
                        unconvertedQualities.erase(unconvertedQualities.begin()+i);
                        return false;
                    }
                }
            }
        }
        return false;
    }

    /**
     * 使用细粒度锁重新实现appendBase
     */
    void appendBase(PosQuality& input, Alignment& a) {
        if(appendReadNameID(input, a)) {
            if(input.converted) {
                lock_guard<mutex> lock(converted_mutex_);
                convertedQualities += input.qual;
            } else {
                lock_guard<mutex> lock(unconverted_mutex_);
                unconvertedQualities += input.qual;
            }
        }
    }
};

/**
 * store all reference position in this class.
 */
class Positions{
public:
    vector<Position*> refPositions; // the pool of all current reference position.
    string chromosome; // current reference chromosome name.
    long long int location; // current location (position) in reference chromosome.
    char lastBase = 'X'; // the last base of reference line. this is for CG_only mode.
    SafeQueue<Position*> freePositionPool; // pool to store free position pointer for reference position.
    SafeQueue<Position*> outputPositionPool; // pool to store the reference position which is loaded and ready to output.
    bool working;
    long long int refCoveredPosition; // this is the last position in reference chromosome we loaded in refPositions.
    ifstream refFile;
    ChromosomeFilePositions chromosomePos; // store the chromosome name and it's streamPos. To quickly find new chromosome in file.
    bool addedChrName = false;
    bool removedChrName = false;

    Positions(string inputRefFileName, bool inputAddedChrName, bool inputRemovedChrName) {
        working = true;
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
        while (refFile.good()) {
            if (!getline(refFile, line)) {
                break;
            }
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

    void appendingFinished() {
        tg.wait();
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
        while (working) {
            if (outputPositionPool.popFront(pos)) {
                *out_ << pos->chromosome << '\t'
                          << to_string(pos->location) << '\t'
                          << pos->strand << '\t'
                          << pos->convertedQualities << '\t'
                          << to_string(pos->convertedQualities.size()) << '\t'
                          << pos->unconvertedQualities << '\t'
                          << to_string(pos->unconvertedQualities.size()) << '\n';
                returnPosition(pos);
            } else {
                this_thread::sleep_for (std::chrono::microseconds(1));
            }
        }
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
                refPositions[index]->uniqueIDs.clear();
                outputPositionPool.push(refPositions[index]);
            }
        }
        refPositions.clear();
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
        while (refFile.good()) {
            getline(refFile, line);
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
    }

    /**
     * load more Position (loadingBlockSize bp) to positions
     * if we meet next chromosome, return false. Else, return ture.
     */
    void loadMore() {
        refCoveredPosition += loadingBlockSize;
        string line;
        while (refFile.good()) {
            getline(refFile, line);
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
     * get a Position pointer from freePositionPool, if freePositionPool is empty, make a new Position pointer.
     */
    void getFreePosition(Position*& newPosition) {
        while (outputPositionPool.size() >= 10000) {
            this_thread::sleep_for (std::chrono::microseconds(1));
        }
        if (freePositionPool.popFront(newPosition)) {
            return;
        } else {
            newPosition = new Position();
        }
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
    void append(string* line) {
        Alignment newAlignment;
        while (refPositions.empty()) {
            this_thread::sleep_for (std::chrono::microseconds(1));
        }
        newAlignment.parse(line);
        delete line;
        appendPositions(newAlignment);
    }
};

#endif //POSITION_3N_TABLE_H
