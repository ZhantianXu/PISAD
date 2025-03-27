

#ifndef SRC_FINGERPRINT_HPP_
#define SRC_FINGERPRINT_HPP_
#include "Options.h"
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#include "vendor/KseqHashIterator.hpp"
#include "vendor/concurrentqueue.h"
#include "vendor/kseq_util.h"
#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class FingerPrint {
public:
  typedef uint16_t AlleleID;
  typedef uint64_t HashedKmer;

  FingerPrint()
      : m_totalCounts(0), m_maxCounts(0), m_totalBases(0), m_totalReads(0) {
    // read in fasta files
    // generate hash table
    initCountsHash();
    if (opt::covThresh != 0) {
      m_maxCounts = (m_counts.size() * opt::covThresh) / 2;
    }
  }

  void computeCounts(const vector<string> &filenames) {
#pragma omp parallel for
    for (unsigned i = 0; i < filenames.size(); ++i) {
      gzFile fp;
      fp = gzopen(filenames[i].c_str(), "r");
      if (fp == Z_NULL) {
#pragma omp critical(stderr)
        {
          std::cerr << "file " << filenames[i] << " cannot be opened"
                    << std::endl;
        }
        exit(1);
      } else if (opt::verbose) {
#pragma omp critical(stderr)
        { std::cerr << "Opening " << filenames[i] << std::endl; }
      }
      // read in seq
      kseq_t *seq = kseq_init(fp);
      int l = kseq_read(seq);
      // cerr << l << m_earlyTerm << endl;
      while (l >= 0 && !m_earlyTerm) {
        processSingleRead(seq);
        l = kseq_read(seq);
        if (opt::verbose > 2) {
#pragma omp critical(stderr)
          if ((m_totalReads % 1000000) == 0) {
            cerr << "Current Total: " << m_totalReads << " reads, "
                 << m_totalKmers << " k-mers, " << m_totalCounts
                 << " total counts, and " << m_totalBases << " total bases "
                 << endl;
          }
        }
      }
      kseq_destroy(seq);
      gzclose(fp);
    }
  }

  void insertCount(const char *seqs, uint64_t seql, unsigned multiplier = 1) {
    for (KseqHashIterator itr(seqs, seql, opt::k); itr != itr.end(); ++itr) {
      if (m_counts.find(*itr) != m_counts.end()) {
#pragma omp atomic update
        // #pragma omp critical
        m_counts[*itr] += multiplier;
        // m_counts.modify_if(*itr, [multiplier](std::pair<const uint64_t,
        // size_t> &item) 				   { item.second +=
        // multiplier;
        // });

#pragma omp atomic update
        m_totalCounts += multiplier;
      }
#pragma omp atomic update
      ++m_totalKmers;
    }
#pragma omp atomic update
    m_totalBases += seql;
  }

  // use only if threads > number of files
  void computeCountsProducerConsumer(const vector<string> &filenames) {
    if (opt::threads <= filenames.size()) {
      // not enough threads to saturate
      computeCounts(filenames);
    } else {
      uint64_t numReads = 0, processedCount = 0;

      moodycamel::ConcurrentQueue<kseq_t> workQueue(opt::threads * s_bulkSize);
      moodycamel::ConcurrentQueue<kseq_t> recycleQueue(opt::threads * s_bulkSize * 2);
      bool good = true;
      typedef std::vector<kseq_t>::iterator iter_t;

      // fill recycleQueue with empty objects
      {
        std::vector<kseq_t> buffer(opt::threads * s_bulkSize * 2, kseq_t());
        recycleQueue.enqueue_bulk(std::move_iterator<iter_t>(buffer.begin()),buffer.size());
      }

#pragma omp parallel
      {
        std::vector<kseq_t> readBuffer(s_bulkSize);
        string outBuffer;
        if (unsigned(omp_get_thread_num()) < filenames.size()) {
          // file reading init
          gzFile fp;
          fp = gzopen(filenames.at(omp_get_thread_num()).c_str(), "r");
          std::cerr << "Opening " << filenames.at(omp_get_thread_num()) << std::endl;

          kseq_t *seq = kseq_init(fp);

          // per thread token
          moodycamel::ProducerToken ptok(workQueue);

          // tokens for recycle queue
          moodycamel::ConsumerToken rctok(recycleQueue);
          moodycamel::ProducerToken rptok(recycleQueue);

          unsigned dequeueSize = recycleQueue.try_dequeue_bulk(rctok, std::move_iterator<iter_t>(readBuffer.begin()),s_bulkSize);
          while (dequeueSize == 0) {
            dequeueSize = recycleQueue.try_dequeue_bulk(rctok, std::move_iterator<iter_t>(readBuffer.begin()),s_bulkSize);
          }

          unsigned size = 0;
          while (kseq_read(seq) >= 0 && !m_earlyTerm) {
            cpy_kseq(&readBuffer[size++], seq);
            if (dequeueSize == size) {
              // try to insert, if cannot queue is full
              while (!workQueue.try_enqueue_bulk(ptok, std::move_iterator<iter_t>(readBuffer.begin()), size)) {
                // try to work
                if (kseq_read(seq) >= 0) {
                  //------------------------WORK CODE
                  // START---------------------------------------
                  processSingleRead(seq);
                  //------------------------WORK CODE
                  // END-----------------------------------------
                } else {
                  goto fileEmpty;
                }
              }
              // reset buffer
              dequeueSize = recycleQueue.try_dequeue_bulk(rctok, std::move_iterator<iter_t>(readBuffer.begin()),s_bulkSize);
              while (dequeueSize == 0) {
                // try to work
                if (kseq_read(seq) >= 0) {
                  //------------------------WORK CODE
                  // START---------------------------------------
                  processSingleRead(seq);
                  //------------------------WORK CODE
                  // END-----------------------------------------
                } else {
                  goto fileEmpty;
                }
                dequeueSize = recycleQueue.try_dequeue_bulk(rctok, std::move_iterator<iter_t>(readBuffer.begin()),s_bulkSize);
              }
              size = 0;
            }
          }
        fileEmpty:
          // finish off remaining work
          for (unsigned i = 0; i < size; ++i) {
            //------------------------WORK CODE
            // START---------------------------------------
            processSingleRead(seq);
            //------------------------WORK CODE
            // END-----------------------------------------
          }
          assert(recycleQueue.enqueue_bulk(rptok, std::move_iterator<iter_t>(readBuffer.begin()), size));
          if (processedCount < numReads) {
            moodycamel::ConsumerToken ctok(workQueue);
            // join in if others are still not finished
            while (processedCount < numReads) {
              size_t num = workQueue.try_dequeue_bulk(ctok, std::move_iterator<iter_t>(readBuffer.begin()),s_bulkSize);
              if (num) {
                for (unsigned i = 0; i < num; ++i) {
                  //------------------------WORK CODE
                  // START---------------------------------------
                  processSingleRead(seq);
                  //------------------------WORK CODE
                  // END-----------------------------------------
                }
                assert(recycleQueue.enqueue_bulk(rptok, std::move_iterator<iter_t>(readBuffer.begin()),num));
              }
            }
          }
#pragma omp atomic update
          good &= false;
          kseq_destroy(seq);
          gzclose(fp);
        } else {
          moodycamel::ConsumerToken ctok(workQueue);
          moodycamel::ProducerToken rptok(recycleQueue);
          while (good) {
            if (workQueue.size_approx() >= s_bulkSize) {
              size_t num = workQueue.try_dequeue_bulk(ctok, std::move_iterator<iter_t>(readBuffer.begin()),s_bulkSize);
              if (num) {
                for (unsigned i = 0; i < num; ++i) {
                  //------------------------WORK CODE
                  // START---------------------------------------
                  processSingleRead(&readBuffer[i]);
                  //------------------------WORK CODE
                  // END-----------------------------------------
                }
                assert(recycleQueue.enqueue_bulk(rptok, std::move_iterator<iter_t>(readBuffer.begin()),num));
              }
            }
          }
        }
      }
    }
  }

  /*
   * Prints summary info needed for computing error rate
   */
  void printOptionalHeader() const {
    string outStr = "";
    outStr += "#@TK\t";
    outStr += std::to_string(m_totalKmers);
    outStr += "\n#@KS\t";
    outStr += std::to_string(opt::k);
    cout << outStr;
  }

  void printCountsMode(std::string &extracted) const {
    for (unsigned y = 0; y < opt::snp.size(); ++y) {
      std::string name1;
      size_t pos = opt::snp[y].find_last_of("/");
      size_t pos1 = opt::snp[y].find_last_of(".");
      if (pos != std::string::npos && pos1 != std::string::npos && pos1 > pos) {
        name1 = opt::snp[y].substr(pos + 1, pos1 - pos - 1);
      } else {
        std::cerr << "file " << opt::snp[y] << " can not be calculated" << std::endl;
        exit(1);
      }
      std::ofstream file(opt::name + extracted + "_" + name1 + ".txt");

      if (file.is_open()) {
        std::cerr << "file " << opt::name + extracted + "_" + name1 + ".txt"
                  << " is opened" << std::endl;

        string outStr = "";
        outStr += "#@TK\t";
        outStr += std::to_string(m_totalKmers);
        outStr += "\n#@KS\t";
        outStr += std::to_string(opt::k);
        file << outStr;

        if (opt::information) {
          file << "\n#locusID\tmaxAT\tmaxCG\tcountAT\tcountCG\tsumAT\tsumCG\tdistinctAT\tdistinctCG\tref\tval\n";
        } else {
          file << "\n#locusID\tmaxAT\tmaxCG\tcountAT\tcountCG\tsumAT\tsumCG\tdistinctAT\tdistinctCG\n";
        }
        string tempStr_ref, tempStr_val;
        for (size_t i = 0; i < m_alleleIDs[y].size(); ++i) {
          outStr.clear();
          const vector<uint64_t> &allele1 = *m_alleleIDToKmerRef[y][i];
          const vector<uint64_t> &allele2 = *m_alleleIDToKmerVar[y][i];
          unsigned countSumAT = 0;
          unsigned countSumCG = 0;

          tempStr_ref.clear();
          tempStr_val.clear();
          for (size_t j = 0; j < allele1.size(); ++j) {
            if (m_counts.find(allele1.at(j)) != m_counts.end()) {
              tempStr_ref += std::to_string(m_counts.at(allele1.at(j)));
              tempStr_ref += ",";
            }
          }

          for (size_t j = 0; j < allele2.size(); ++j) {
            if (m_counts.find(allele2.at(j)) != m_counts.end()) {
              tempStr_val += std::to_string(m_counts.at(allele2.at(j)));
              tempStr_val += ",";
            }
          }

          auto calculateModeAdjusted = [&](const vector<uint64_t> &allele,unsigned &countSum) -> int {
            unordered_map<uint64_t, unsigned> freqMap;
            for (auto id : allele) {
              unsigned freq = m_counts.at(id);
              ++freqMap[freq];  // freqMap计算频数  数字：频数
              countSum += freq; // Accumulate counts for sum
            }

            unsigned maxCount = 0;
            uint64_t mode = 0;
            for (const auto &p : freqMap) {
              if (p.second > maxCount) {
                maxCount = p.second;
                mode = p.first; // 取最高频数的数字
              }
            }

            size_t countInRange = 0;
            uint64_t sumInRange = 0;
            for (auto id : allele) {
              unsigned freq = m_counts.at(id);
              if (mode <= 2) // 防止为负数
              {
                if (freq <= mode + 2) {
                  ++countInRange;
                  sumInRange += freq;
                }
              } else {
                if (freq >= mode - 2 &&
                    freq <= mode + 2) // 最高频数的数字±2，的个数之和
                {
                  ++countInRange;
                  sumInRange += freq;
                }
              }
            }

            if (countInRange > allele.size() / 2) // 如果个数大于一半
            {
              return sumInRange / countInRange; // Floor division
            } else {
              std::vector<int> copy;
              for (auto id : allele) {
                unsigned freq = m_counts.at(id);
                copy.push_back(freq);
              }
              std::sort(copy.begin(), copy.end());
              size_t size = copy.size();

              if (size % 2 == 0) {
                return (copy[size / 2 - 1] + copy[size / 2]) / 2;
              } else {
                return copy[size / 2];
              }
            }
          };

          auto modelmax = [&](const vector<uint64_t> &allele) -> int {
            unsigned covmax = 0;
            for (auto id : allele) {
              unsigned freq = m_counts.at(id);
              if (freq > covmax) {
                covmax = freq;
              }
            }
            return covmax;
          };

          if (allele1.size() >= 4 && allele2.size() >= 4) {
            int modeAdjustedAT = calculateModeAdjusted(allele1, countSumAT);
            int modeAdjustedCG = calculateModeAdjusted(allele2, countSumCG);
            int maxAT = modelmax(allele1);
            int maxCG = modelmax(allele2);

            outStr += m_alleleIDs[y].at(i);
            outStr += "\t";
            outStr += std::to_string(maxAT);
            outStr += "\t";
            outStr += std::to_string(maxCG);
            outStr += "\t";
            outStr += std::to_string(modeAdjustedAT);
            outStr += "\t";
            outStr += std::to_string(modeAdjustedCG);
            outStr += "\t";
            outStr += std::to_string(countSumAT);
            outStr += "\t";
            outStr += std::to_string(countSumCG);
            outStr += "\t";
            outStr += std::to_string(allele1.size());
            outStr += "\t";
            outStr += std::to_string(allele2.size());
            if (opt::information) {
              outStr += "\t";
              outStr += tempStr_ref;
              outStr += "\t";
              outStr += tempStr_val;
              outStr += "\n";
            } else {
              outStr += "\n";
            }
            file << outStr;
          }
        }
      } else {
        std::cerr << "file " << opt::snp[y] << " can not be opened"
                  << std::endl;
        exit(1);
      }
    }
  }
  uint64_t getTotalKmerCounts() { return m_totalCounts; }

  uint64_t getTotalCounts() { return m_totalBases; }

  string printInfoSummary() {

    string outStr = "";
    outStr += "Total Bases Considered: ";
    outStr += std::to_string(getTotalCounts());
    outStr += "\n";
    outStr += "Total k-mers Considered: ";
    outStr += std::to_string(m_totalKmers);
    outStr += "\n";
    outStr += "Total k-mers Recorded: ";
    outStr += std::to_string(getTotalKmerCounts());
    outStr += "\n";
    outStr += "Distinct k-mers in initial set: ";
    outStr += std::to_string(m_counts.size());
    outStr += "\n";
    outStr += "Total Sites: ";
    outStr += std::to_string(m_alleleIDToKmerRef.size());
    outStr += "\n";
    return (outStr);
  }

private:
  const static size_t s_bulkSize = 1024;
  uint64_t m_totalCounts;
  uint64_t m_totalKmers;
  uint64_t m_maxCounts;
  // using Map1 = phmap::parallel_flat_hash_map<uint64_t, size_t,
  // 										   std::hash<uint64_t>,
  // 										   std::equal_to<uint64_t>,
  // 										   std::allocator<std::pair<const
  // uint64_t,
  // size_t>>, 6, std::mutex>;
  // // Map1 m_counts;
  // tsl::robin_map<uint64_t, size_t, robin_hood::hash<uint64_t>> m_counts; //
  // k-mer to count
  tsl::robin_map<uint64_t, size_t> m_counts; // k-mer to count
  // vector<shared_ptr<vector<HashedKmer>>> m_alleleIDToKmerRef;
  // vector<shared_ptr<vector<HashedKmer>>> m_alleleIDToKmerVar;
  vector<vector<shared_ptr<vector<HashedKmer>>>> m_alleleIDToKmerRef;
  vector<vector<shared_ptr<vector<HashedKmer>>>> m_alleleIDToKmerVar;
  vector<vector<string>> m_alleleIDs;
  uint64_t m_totalBases;
  uint64_t m_totalReads; // for debugging purposes
  bool m_earlyTerm = false;
  //	unsigned m_maxSiteSize;
  //	static const unsigned interval = 65536;

  void processSingleRead(kseq_t *seq) {
    // k-merize and insert
    insertCount(seq->seq.s, seq->seq.l);
    // cerr << m_maxCounts << m_totalCounts << m_earlyTerm << endl;
    m_totalReads++;
    if (m_maxCounts != 0 && m_totalCounts > m_maxCounts) {
      if (opt::verbose > 0) {
#pragma omp critical(stderr)
        {
          cerr << "max count reached at " << m_totalReads << " reads, "
               << m_totalKmers << " k-mers, " << m_totalCounts
               << " total counts, and " << m_totalBases << " total bases "
               << endl;
        }
      }
      m_earlyTerm = true;
    }
  }

  void initCountsHash() {
    m_alleleIDToKmerRef.resize(opt::snp.size());
    m_alleleIDToKmerVar.resize(opt::snp.size());
    m_alleleIDs.resize(opt::snp.size());
    // #pragma omp parallel for
    for (unsigned i = 0; i < opt::snp.size(); ++i) {
      gzFile fp;
      fp = gzopen(opt::snp[i].c_str(), "r");
      if (fp == Z_NULL) {
#pragma omp critical(stderr)
        {
          std::cerr << "file " << opt::snp[i] << " cannot be opened"
                    << std::endl;
        }
        exit(1);
      } else if (opt::verbose) {
#pragma omp critical(stderr)
        { std::cerr << "Opening " << opt::snp[i] << std::endl; }
      }

      kseq_t *seq = kseq_init(fp);
      int l = kseq_read(seq);
      size_t entryNum = 0;

      // m_alleleIDToKmerRef.emplace_back(vector<shared_ptr<vector<HashedKmer>>>());
      // m_alleleIDToKmerVar.emplace_back(vector<shared_ptr<vector<HashedKmer>>>());

      while (l >= 0) {
        if (entryNum % 2 == 0) {
          unsigned index = entryNum / 2;
          m_alleleIDToKmerRef[i].emplace_back(
              shared_ptr<vector<uint64_t>>(new vector<uint64_t>()));
          // k-merize and
          for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
               itr != itr.end(); ++itr) {
            uint64_t hv = *itr;
            // check for duplicates
            m_alleleIDToKmerRef[i][index]->emplace_back(hv);
            m_counts[hv] = 0;
          }
          m_alleleIDs[i].emplace_back(seq->name.s);
        } else {
          unsigned index = entryNum / 2;
          m_alleleIDToKmerVar[i].emplace_back(
              shared_ptr<vector<uint64_t>>(new vector<uint64_t>()));
          // k-merize and insert
          for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
               itr != itr.end(); ++itr) {
            uint64_t hv = *itr;
            // check for duplicates
            m_alleleIDToKmerVar[i][index]->emplace_back(hv);
            m_counts[hv] = 0;
          }
        }
        l = kseq_read(seq);
        entryNum++;
      }
      kseq_destroy(seq);
      gzclose(fp);

      cerr << "m_alleleIDs " << m_alleleIDs[i].size() << endl;
      cerr << "m_alleleIDToKmerRef " << m_alleleIDToKmerRef[i].size() << endl;
      cerr << "m_alleleIDToKmerVar " << m_alleleIDToKmerVar[i].size() << endl;
      cerr << "mcounts " << m_counts.size() << endl;
    }
  }
};
#endif /* SRC_FINGERPRINT_HPP_ */
