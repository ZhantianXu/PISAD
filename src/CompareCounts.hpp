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
        m_totalCounts(vector<uint64_t>(filenames.size(), 0)) {
    if (opt::verbose > 0) {
      cerr << "Reading count files" << endl;
    }
    // read first file twice to init vectors
    {
      ifstream fh(m_filenames.at(0));
      string line;
      if (fh.is_open()) {
        while (getline(fh, line)) {
          if (line.length() > 0) {
            stringstream ss;
            ss.str(line);
            string item;
            getline(ss, item, '\t');
            if (line.at(0) != '#') {
              string locusID = item;
              m_locusIDToIndex[locusID] = m_locusIDs.size();
              m_locusIDs.emplace_back(locusID);
              getline(ss, item, '\t');
              getline(ss, item, '\t');

              getline(ss, item, '\t');
              getline(ss, item, '\t');

              getline(ss, item, '\t');
              getline(ss, item, '\t');

              getline(ss, item, '\t');
              m_distinct.push_back(loadPair(ss, item));
            }
          }
        }
      }
    }
    m_counts = PairedCount(filenames.size(),vector<pair<unsigned, unsigned>>(m_distinct.size()));
    m_counts1 = PairedCount(filenames.size(), vector<pair<unsigned, unsigned>>(m_distinct.size()));
    m_sum = PairedCount(filenames.size(),vector<pair<unsigned, unsigned>>(m_distinct.size()));

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
            stringstream ss;
            ss.str(line);
            string item;
            getline(ss, item, '\t');
            if (line.at(0) == '#') {
              if (item == "#@TK") {
                getline(ss, item, '\t');
                m_rawTotalCounts[i] = std::stoull(item);
              } else if (item == "#@KS") {
                getline(ss, item, '\t');
                m_kmerSize[i] = std::stoull(item);
              }
            } else {
              string locusID = item;
              // locusID\t maxAT\tmaxCG\t countAT\tcountCG\t sumAT\tsumCG\t
              // distinctAT\tdistinctCG\n
              getline(ss, item, '\t');
              m_counts[i][m_locusIDToIndex.at(locusID)] = loadPair(ss, item);
              m_counts1[i][m_locusIDToIndex.at(locusID)] = loadPair(ss, item);
              m_totalCounts[i] +=
                  m_counts[i][m_locusIDToIndex.at(locusID)].first +
                  m_counts[i][m_locusIDToIndex.at(locusID)].second;
              m_sum[i][m_locusIDToIndex.at(locusID)] = loadPair(ss, item);
            }
          }
        }
      }
    }
  }

  /*
   * Compute comparisons between all combinations
   * Slower than using PCA to lower comparison numbers
   */
  void computeScore() {
    cout << m_header;
    vector<GenotypeSummary> genotype(m_totalCounts.size());
    for (unsigned i = 0; i < m_totalCounts.size(); ++i) {
      genotype[i].cov = double(m_totalCounts[i]) / double(m_distinct.size());
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

  ~CompareCounts() {
    // TODO Auto-generated destructor stub
  }

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
    //		double madFreq = 0;
    //		double median = 0;
    //		double mean = 0;
    //		double var50 = 0;
    double errorRate = 0;
    double cov = 0;
    unsigned searchDegree = 0;
  };

  const vector<string> &m_filenames;
  typedef vector<vector<pair<unsigned, unsigned>>> PairedCount;
  typedef std::vector<std::vector<double>> vector_of_vectors_t;
  typedef KDTreeVectorOfVectorsAdaptor<vector_of_vectors_t, double> kd_tree_t;

  vector<double> m_sumlogPSingle;
  PairedCount m_counts;
  PairedCount m_counts1;
  PairedCount m_sum;
  vector<pair<unsigned, unsigned>> m_distinct;
  vector<uint64_t> m_rawTotalCounts;
  vector<unsigned> m_kmerSize;
  vector<uint64_t> m_totalCounts;
  vector<string> m_locusIDs;
  tsl::robin_map<string, unsigned> m_locusIDToIndex;
  tsl::robin_map<string, unsigned> m_filenameToID;
  const string m_header = "sample\tscore\tsame\trelate\tvalid_relate_sites\tvalid_sites\tsites\tcov\t";

  struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &pair) const {
      return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
  };

  GenotypeSummary
  calcHomHetMiss(const vector<pair<unsigned, unsigned>> &counts) {
    GenotypeSummary count = {};
    vector<double> hetCount; // heterozygous frequency
    for (size_t i = 0; i < counts.size(); ++i) {
      if (counts.at(i).first > opt::minCov) {
        if (counts.at(i).second > opt::minCov) {
          ++count.hets;
          hetCount.push_back(double(counts.at(i).first) /double(counts.at(i).first + counts.at(i).second));
        } else {
          ++count.homs;
        }
      } else if (counts.at(i).second > opt::minCov) {
        ++count.homs;
      } else {
        ++count.miss;
      }
    }
    //		count.median = 0.5;
    //		count.madFreq = calculateMAD(hetCount, count.median);
    //		count.mean = 0.5;
    //		count.var50 = calculateVar50(hetCount, count.mean);
    return (count);
  }

  /*
   * prepare results string
   */
  void resultsStr1(string &temp, const vector<GenotypeSummary> &genotype,
                   const Relate &info, double score, unsigned i, unsigned j) {
    temp.clear();
    temp += m_filenames[j];
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
    temp += to_string(m_distinct.size());
    temp += "\t";
    temp += to_string(genotype.at(j).cov);
    temp += "\t";
  }

  pair<unsigned, unsigned> loadPair(stringstream &ss, string &item) {
    unsigned count1 = std::stoul(item);
    getline(ss, item, '\t');
    unsigned count2 = std::stoul(item);
    getline(ss, item, '\t');
    return (std::make_pair(count1, count2));
  }

  double computeSumLogPSingle(unsigned index,const vector<unsigned> &pos) const {
    double sumLogP = 0;
    for (vector<unsigned>::const_iterator i = pos.begin(); i != pos.end();
         ++i) {
      double freqAT = 0;
      double freqCG = 0;
      if (m_counts1.at(index).at(*i).first > 0) // opt::minCov
      {
        freqAT = double(m_counts1.at(index).at(*i).first) /
                 double(m_counts1.at(index).at(*i).first + m_counts1.at(index).at(*i).second);
      }
      if (m_counts1.at(index).at(*i).second > 0) {
        freqCG = double(m_counts1.at(index).at(*i).second) /
                 double(m_counts1.at(index).at(*i).first + m_counts1.at(index).at(*i).second);
      }
      sumLogP += m_counts1.at(index).at(*i).first * freqAT +
                 m_counts1.at(index).at(*i).second * freqCG + 1;
    }
    return (sumLogP);
  }

  // opt::minCov
  double computeSumLogPJoint(unsigned index1, const vector<unsigned> &pos,
                             unsigned covThresh = 0) const {
    double sumLogP = 0;
    for (vector<unsigned>::const_iterator i = pos.begin(); i != pos.end();
         ++i) {
      double freqAT = 0;
      double freqCG = 0;
      unsigned countAT = m_counts1.at(index1).at(*i).first + 1;
      unsigned countCG = m_counts1.at(index1).at(*i).second + 1;
      if (countAT > covThresh) {
        freqAT = double(countAT) / double(countAT + countCG);
      }
      if (countCG > covThresh) {
        freqCG = double(countCG) / double(countAT + countCG);
      }
      sumLogP += countAT * freqAT + countCG * freqCG;
    }
    return (sumLogP);
  }

  // TODO: Easily parelizable
  vector<unsigned> gatherValidEntries(unsigned index1) {
    vector<unsigned> valid;
    vector<bool> binValid(m_distinct.size(), true);
    for (unsigned j = 0; j < m_distinct.size(); ++j) {
      if (m_counts1.at(index1).at(j).first + m_counts1.at(index1).at(j).second <
          2) {
        binValid[j] = false;
      }
    }
    for (unsigned i = 0; i < binValid.size(); ++i) {
      if (binValid[i]) {
        valid.push_back(i);
      }
    }
    return (valid);
  }

  // compute only sites that aren't missing
  double computeLogLikelihood(unsigned index1,const vector<unsigned> &validIndexes) {
    return -2.0 * (computeSumLogPJoint(index1, validIndexes) -
                   computeSumLogPSingle(index1, validIndexes));
  }

  void initLogPSum(const vector<unsigned> &pos) {
    for (unsigned i = 0; i < m_counts.size(); ++i) {
      m_sumlogPSingle[i] = computeSumLogPSingle(i, pos);
    }
  }

  Relate calcRelatedness1(unsigned index2,const vector<unsigned> &validIndexes) {
    Relate info = {};
    vector<double> ratios;

    for (vector<unsigned>::const_iterator i = validIndexes.begin();
         i != validIndexes.end(); ++i) {
      if (m_counts1.at(index2).at(*i).first > 2 or
          m_counts1.at(index2).at(*i).second > 2) {
        double ratio = m_counts1.at(index2).at(*i).first * 1.0 /
                       (m_counts1.at(index2).at(*i).first + m_counts1.at(index2).at(*i).second);
        ratios.push_back(ratio);
      }
    }

    const int num_bins = 15;
    vector<int> bin_count(num_bins, 0);
    double bin_size = 0.5 / num_bins * 1.0;

    for (double ratio : ratios) {
      if (ratio < 0.5) {
        for (int i = 0; i < num_bins; i++) {
          if ((ratio < (i + 1.0) * bin_size)) {
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
    for (vector<unsigned>::const_iterator i = validIndexes.begin();
         i != validIndexes.end(); ++i) {
      if (m_counts1.at(index2).at(*i).first > 2 or m_counts1.at(index2).at(*i).second > 2) {

        if (m_counts1.at(index2).at(*i).first * 1.0 /
                (m_counts1.at(index2).at(*i).first + m_counts1.at(index2).at(*i).second) < info.threshold) {
          ++info.homs2;
        }
        if (m_counts1.at(index2).at(*i).second * 1.0 /
                (m_counts1.at(index2).at(*i).first +m_counts1.at(index2).at(*i).second) < info.threshold) {
          ++info.homs2;
        }

        if (m_counts1.at(index2).at(*i).first * 1.0 /
                    (m_counts1.at(index2).at(*i).first + m_counts1.at(index2).at(*i).second) > info.threshold &&
            m_counts1.at(index2).at(*i).first * 1.0 /
                    (m_counts1.at(index2).at(*i).first + m_counts1.at(index2).at(*i).second) < upper_threshold) {
          ++info.hets2;
        }
      } else {
        info.nullpoint++;
      }
    }
    info.relatedness = double(info.hets2) / double(info.hets2 + info.homs2);
    return info;
  }

  double computeErrorRate(unsigned index) const {
    if (m_rawTotalCounts.at(index) > 0 && m_kmerSize.at(index) > 0) {
      uint64_t sum = 0;
      uint64_t distinctKmers = 0;
      for (unsigned i = 0; i < m_distinct.size(); ++i) {
        sum += m_sum.at(index).at(i).first + m_sum.at(index).at(i).second;
        distinctKmers += m_distinct.at(i).first + m_distinct.at(i).second;
      }
      double expected = double(m_rawTotalCounts.at(index)) *
                        double(distinctKmers) / double(opt::genomeSize);
      return (1.0 -
              pow(double(sum) / expected, 1.0 / double(m_kmerSize.at(index))));
    } else {
      return -1.0;
    }
  }
};

#endif /* SRC_COMPARECOUNTS_HPP_ */
