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
#include "alignment_3n_table.h"

using namespace std;

extern bool CG_only;
extern long long int loadingBlockSize;

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
    mutex mutex_;
public:
    string chromosome; // reference chromosome name
    long long int location; // 1-based position
    char strand; // +(REF) or -(REF-RC)
    string convertedQualities; // each char is a mapping quality on this position for converted base.
    string unconvertedQualities; // each char is a mapping quality on this position for unconverted base.
    vector<uniqueID> uniqueIDs; // each value represent a readName which contributed the base information.
                              // readNameIDs is to make sure no read contribute 2 times in same position.

    void initialize() {
        chromosome.clear();
        location = -1;
        strand = '?';
        convertedQualities.clear();
        unconvertedQualities.clear();
        vector<uniqueID>().swap(uniqueIDs);
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
     * binary search of readNameID in readNameIDs.
     * always return a index.
     * if cannot find, return the index which has bigger value than input readNameID.
     */
    int searchReadNameID (unsigned long long&readNameID, int start, int end) {
        if (uniqueIDs.empty()) {
            return 0;
        }
        if (start <= end) {
            int middle = (start + end) / 2;
            if (uniqueIDs[middle].readNameID == readNameID) {
                return middle;
            }
            if (uniqueIDs[middle].readNameID > readNameID) {
                return searchReadNameID(readNameID, start, middle-1);
            }
            return searchReadNameID(readNameID, middle+1, end);
        }
        return start; // return the bigger one
    }


    /**
     * with a input readNameID, add it into readNameIDs.
     * if the input readNameID already exist in readNameIDs, return false.
     */
    bool appendReadNameID(PosQuality& InBase, Alignment& InAlignment) {
        int idCount = uniqueIDs.size();
        if (idCount == 0 || InAlignment.readNameID > uniqueIDs.back().readNameID) {
            uniqueIDs.emplace_back(InAlignment.readNameID, InBase.converted, InBase.qual);
            return true;
        }
        int index = searchReadNameID(InAlignment.readNameID, 0, idCount);
        if (uniqueIDs[index].readNameID == InAlignment.readNameID) {
            // if the new base is consistent with exist base's conversion status, ignore
            // otherwise, delete the exist conversion status
            if (uniqueIDs[index].removed) {
                return false;
            }
            if (uniqueIDs[index].isConverted != InBase.converted) {
                uniqueIDs[index].removed = true;
                if (uniqueIDs[index].isConverted) {
                    for (int i = 0; i < convertedQualities.size(); i++) {
                        if (convertedQualities[i] == InBase.qual) {
                            convertedQualities.erase(convertedQualities.begin()+i);
                            return false;
                        }
                    }
                } else {
                    for (int i = 0; i < unconvertedQualities.size(); i++) {
                        if (unconvertedQualities[i] == InBase.qual) {
                            unconvertedQualities.erase(unconvertedQualities.begin()+i);
                            return false;
                        }
                    }
                }
            }
            return false;
        } else {
            uniqueIDs.emplace(uniqueIDs.begin()+index, InAlignment.readNameID, InBase.converted, InBase.qual);
            return true;
        }
    }

    /**
     * append the SAM information into this position.
     */
    void appendBase (PosQuality& input, Alignment& a) {
        mutex_.lock();
        if (appendReadNameID(input,a)) {
            if (input.converted) {
                convertedQualities += input.qual;
            } else {
                unconvertedQualities += input.qual;
            }
        }
        mutex_.unlock();
    }
};

#endif //POSITION_3N_TABLE_H
