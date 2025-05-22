#ifndef SRC_COMPARECOUNTS_HPP_
#define SRC_COMPARECOUNTS_HPP_

#include <algorithm>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "KDTreeUtil.h"
#include "src/Options.h"
#include "vendor/kfunc.c"
#include "vendor/nanoflann.hpp"
#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#include <fstream>
#include <sstream>

using namespace std;

class CompareCounts {
public:
    CompareCounts(const vector<string> &filenames)
        : m_filenames(filenames),
          m_sumlogPSingle(vector<double>(filenames.size())),
          m_rawTotalCounts(vector<uint64_t>(filenames.size(), 0)),
          m_kmerSize(vector<unsigned>(filenames.size(), 0)),
          m_totalCounts(vector<uint64_t>(filenames.size(), 0)),
          m_counts(filenames.size()),
          m_counts1(filenames.size()),
          m_sum(filenames.size()),
          m_distinct(filenames.size()) {
        if (opt::verbose > 0) {
            cerr << "Reading count files" << endl;
        }

        // 并行读取每个文件，独立存储数据
#pragma omp parallel for
        for (unsigned i = 0; i < m_filenames.size(); ++i) {
            if (opt::verbose > 1) {
#pragma omp critical(cerr)
                cerr << "Reading: " << m_filenames.at(i) << endl;
            }
            ifstream fh(m_filenames.at(i));
#pragma omp critical(m_filenameToID)
            { m_filenameToID[m_filenames.at(i)] = i; }
            string line;
            if (fh.is_open()) {
                while (getline(fh, line)) {
                    if (line.length() > 0) {
                        stringstream ss(line);
                        string item;
                        getline(ss, item, '\t'); // locusID，忽略
                        if (line.at(0) == '#') {
                            // if (item == "#@TK") {
                            //     getline(ss, item, '\t');
                            //     m_rawTotalCounts[i] = std::stoull(item);
                            // } else if (item == "#@KS") {
                            //     getline(ss, item, '\t');
                            //     m_kmerSize[i] = std::stoull(item);
                            // }
                        } else {
                            // 直接按行存储数据，不依赖 locusID
                            getline(ss, item, '\t'); // maxAT
                            pair<unsigned, unsigned> countPair = loadPair(ss, item); // countAT, countCG
                            m_counts[i].push_back(countPair);
                            // m_counts1[i].push_back(countPair);
                            m_counts1[i].push_back(loadPair(ss, item));
                            m_totalCounts[i] += countPair.first + countPair.second;
                            m_sum[i].push_back(loadPair(ss, item)); // sumAT, sumCG
                            m_distinct[i].push_back(loadPair(ss, item)); // distinctAT, distinctCG
                        }
                    }
                }
                fh.close();
            }
        }
    }

    void computeScore() {
        cout << m_header;
        vector<GenotypeSummary> genotype(m_totalCounts.size());
        for (unsigned i = 0; i < m_totalCounts.size(); ++i) {
            genotype[i].cov = double(m_totalCounts[i]) / double(m_distinct[i].size());
        }
        string temp = "\n";
        cout << temp;
#pragma omp parallel for private(temp)
        for (size_t j = 0; j < m_counts.size(); j++) {
            vector<unsigned> validIndexes = gatherValidEntries(j);
            double score = std::numeric_limits<double>::max();
            if (validIndexes.size() > 0) {
                score = computeLogLikelihood(j, validIndexes);
                score /= double(validIndexes.size());
            }
            Relate info = calcRelatedness1(j, validIndexes);
            resultsStr1(temp, genotype, info, score, validIndexes.size(), j);
            temp += "\n";
#pragma omp critical(cout)
            { cout << temp; }
        }
    }

    ~CompareCounts() {}

private:
    struct Relate {
        double relatedness = 0;
        unsigned ibs0 = 0;
        unsigned ibs2 = 0;
        unsigned ibs1 = 0;
        double homConcord = 0;
        unsigned sharedHoms = 0;
        unsigned sharedHets = 0;
        unsigned hets1 = 0;
        unsigned homs1 = 0;
        unsigned hets2 = 0;
        unsigned homs2 = 0;
        unsigned nullpoint = 0;
        double threshold = 0.0;
        int bin = 0;
        int bin_index = 0;
    };

    struct GenotypeSummary {
        unsigned hets = 0;
        unsigned homs = 0;
        unsigned miss = 0;
        double errorRate = 0;
        double cov = 0;
        unsigned searchDegree = 0;
    };

    const vector<string> &m_filenames;
    typedef vector<vector<pair<unsigned, unsigned>>> PairedCount;
    typedef std::vector<std::vector<double>> vector_of_vectors_t;
    typedef KDTreeVectorOfVectorsAdaptor<vector_of_vectors_t, double> kd_tree_t;

    vector<double> m_sumlogPSingle;
    PairedCount m_counts;    // 每个文件的计数数据
    PairedCount m_counts1;   // 副本
    PairedCount m_sum;       // 每个文件的 sum 数据
    PairedCount m_distinct;  // 每个文件的 distinct 数据
    vector<uint64_t> m_rawTotalCounts;
    vector<unsigned> m_kmerSize;
    vector<uint64_t> m_totalCounts;
    tsl::robin_map<string, unsigned> m_filenameToID;
    const string m_header = "sample\tscore\tsame\trelate\tvalid_relate_sites\tvalid_sites\tsites\tcov\t";

    GenotypeSummary calcHomHetMiss(const vector<pair<unsigned, unsigned>> &counts) {
        GenotypeSummary count = {};
        vector<double> hetCount;
        for (size_t i = 0; i < counts.size(); ++i) {
            if (counts[i].first > opt::minCov) {
                if (counts[i].second > opt::minCov) {
                    ++count.hets;
                    hetCount.push_back(double(counts[i].first) / double(counts[i].first + counts[i].second));
                } else {
                    ++count.homs;
                }
            } else if (counts[i].second > opt::minCov) {
                ++count.homs;
            } else {
                ++count.miss;
            }
        }
        return count;
    }

    void resultsStr1(string &temp, const vector<GenotypeSummary> &genotype,
                     const Relate &info, double score, unsigned i, unsigned j) {
        temp.clear();
        temp += m_filenames[j].substr(m_filenames[j].find_last_of('/') + 1);
        temp += "\t";
        temp += to_string(score);
        if (score < opt::scoreThresh) {
            temp += "\t1\t";
        } else {
            temp += "\t0\t";
        }
        temp += to_string(info.relatedness);
        temp += "\t";
        temp += to_string((info.hets2 + info.homs2));
        temp += "\t";
        temp += to_string(i);
        temp += "\t";
        temp += to_string(m_distinct[j].size());
        temp += "\t";
        temp += to_string(genotype[j].cov);
        temp += "\t";
    }

    pair<unsigned, unsigned> loadPair(stringstream &ss, string &item) {
        unsigned count1 = std::stoul(item);
        getline(ss, item, '\t');
        unsigned count2 = std::stoul(item);
        getline(ss, item, '\t');
        return std::make_pair(count1, count2);
    }

    double computeSumLogPSingle(unsigned index, const vector<unsigned> &pos) const {
        double sumLogP = 0;
        for (const auto& i : pos) {
            double freqAT = 0;
            double freqCG = 0;
            if (m_counts1[index][i].first > 0) {
                freqAT = double(m_counts1[index][i].first) /
                         double(m_counts1[index][i].first + m_counts1[index][i].second);
            }
            if (m_counts1[index][i].second > 0) {
                freqCG = double(m_counts1[index][i].second) /
                         double(m_counts1[index][i].first + m_counts1[index][i].second);
            }
            sumLogP += m_counts1[index][i].first * freqAT +
                       m_counts1[index][i].second * freqCG + 1;
        }
        return sumLogP;
    }

    double computeSumLogPJoint(unsigned index, const vector<unsigned> &pos, unsigned covThresh = 0) const {
        double sumLogP = 0;
        for (const auto& i : pos) {
            double freqAT = 0;
            double freqCG = 0;
            unsigned countAT = m_counts1[index][i].first + 1;
            unsigned countCG = m_counts1[index][i].second + 1;
            if (countAT > covThresh) {
                freqAT = double(countAT) / double(countAT + countCG);
            }
            if (countCG > covThresh) {
                freqCG = double(countCG) / double(countAT + countCG);
            }
            sumLogP += countAT * freqAT + countCG * freqCG;
        }
        return sumLogP;
    }

    vector<unsigned> gatherValidEntries(unsigned index) {
        vector<unsigned> valid;
        for (unsigned j = 0; j < m_distinct[index].size(); ++j) {
            if (m_counts1[index][j].first + m_counts1[index][j].second >= 2) {
                valid.push_back(j);
            }
        }
        return valid;
    }

    double computeLogLikelihood(unsigned index, const vector<unsigned> &validIndexes) {
        return -2.0 * (computeSumLogPJoint(index, validIndexes) -
                       computeSumLogPSingle(index, validIndexes));
    }

    Relate calcRelatedness1(unsigned index, const vector<unsigned> &validIndexes) {
        Relate info = {};
        vector<double> ratios;

        for (const auto& i : validIndexes) {
            if (m_counts1[index][i].first > 2 || m_counts1[index][i].second > 2) {
                double ratio = m_counts1[index][i].first * 1.0 /
                               (m_counts1[index][i].first + m_counts1[index][i].second);
                ratios.push_back(ratio);
            }
        }

        const int num_bins = 15;
        vector<int> bin_count(num_bins, 0);
        double bin_size = 0.5 / num_bins * 1.0;

        for (double ratio : ratios) {
            if (ratio < 0.5) {
                for (int i = 0; i < num_bins; i++) {
                    if (ratio < (i + 1.0) * bin_size) {
                        ++bin_count[i];
                        break;
                    }
                }
            }
        }

        info.bin_index = 0;
        info.bin = bin_count[0];
        for (int i = 0; i < num_bins; ++i) {
            if (bin_count[i] < info.bin) {
                info.bin = bin_count[i];
                info.bin_index = i;
            }
        }

        info.threshold = (info.bin_index + 0.5) * bin_size * 1.0;
        double upper_threshold = 1.0 - info.threshold;
        for (const auto& i : validIndexes) {
            if (m_counts1[index][i].first > 2 || m_counts1[index][i].second > 2) {
                double ratio = m_counts1[index][i].first * 1.0 /
                               (m_counts1[index][i].first + m_counts1[index][i].second);
                if (ratio < info.threshold) {
                    ++info.homs2;
                }
                if (m_counts1[index][i].second * 1.0 /
                    (m_counts1[index][i].first + m_counts1[index][i].second) < info.threshold) {
                    ++info.homs2;
                }
                if (ratio > info.threshold && ratio < upper_threshold) {
                    ++info.hets2;
                }
            } else {
                info.nullpoint++;
            }
        }
        info.relatedness = (double(info.hets2 + info.homs2) != 0) ? double(info.hets2) / double(info.hets2 + info.homs2) : -1;
        return info;
    }

    double computeErrorRate(unsigned index) const {
        if (m_rawTotalCounts[index] > 0 && m_kmerSize[index] > 0) {
            uint64_t sum = 0;
            uint64_t distinctKmers = 0;
            for (unsigned i = 0; i < m_distinct[index].size(); ++i) {
                sum += m_sum[index][i].first + m_sum[index][i].second;
                distinctKmers += m_distinct[index][i].first + m_distinct[index][i].second;
            }
            double expected = double(m_rawTotalCounts[index]) *
                              double(distinctKmers) / double(opt::genomeSize);
            return (1.0 - pow(double(sum) / expected, 1.0 / double(m_kmerSize[index])));
        } else {
            return -1.0;
        }
    }
};

#endif /* SRC_COMPARECOUNTS_HPP_ */