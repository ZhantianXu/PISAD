

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
#include "vendor/FastxParser.hpp"
#include <thread>
#include <vector>

#include "parallel_hashmap/btree.h"
#include "parallel_hashmap/meminfo.h"
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_base.h"
#include "parallel_hashmap/phmap_bits.h"
#include "parallel_hashmap/phmap_config.h"
#include "parallel_hashmap/phmap_dump.h"
#include "parallel_hashmap/phmap_fwd_decl.h"
#include "parallel_hashmap/phmap_utils.h"

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
    if(opt::threads<=filenames.size() && opt::threads!=0){
      computeCountsSingle(filenames);
    }
    else{
      processWithFastxParseropenmp(filenames);
    }
  }

  void computeCountsSingle(const vector<string> &filenames) {
#pragma omp parallel for
    for (unsigned i = 0; i < filenames.size(); ++i) {
      gzFile fp;
      fp = gzopen(filenames[i].c_str(), "r");
      if (fp == Z_NULL) {
#pragma omp critical(stderr)
        {
          std::cerr << "file " << filenames[i] << " cannot be opened"<< std::endl;
        }
        exit(1);
      } else if (opt::verbose) {
#pragma omp critical(stderr)
        { std::cerr << "Opening " << filenames[i] << std::endl; }
      }
      kseq_t *seq = kseq_init(fp);
      int l = kseq_read(seq);
      while (l >= 0 && !m_earlyTerm) {
        insertCount(seq->seq.s, seq->seq.l);
        l = kseq_read(seq);
      }
      kseq_destroy(seq);
      gzclose(fp);
    }
  }

  void insertCount(const char *seqs, uint64_t seql, unsigned multiplier = 1) {
    for (KseqHashIterator itr(seqs, seql, opt::k); itr != itr.end(); ++itr) {
      if (m_counts.find(*itr) != m_counts.end()) {
  #pragma omp critical
        m_counts[*itr] += multiplier;

#pragma omp atomic update
        m_totalCounts += multiplier;
      }
    }
  }

  void processWithFastxParseropenmp(const vector<string> &filenames) {
    size_t numProducers,numConsumers;
    if(opt::threads==0){
      unsigned int cores = std::thread::hardware_concurrency();
      if(cores>=filenames.size()*6){
        numProducers = filenames.size();
        numConsumers = filenames.size()*5;
      }else{
        numProducers = cores/6;
        numConsumers = numProducers*5;
      }
    }else{
      numProducers = (opt::threads / 6 < filenames.size()) ? (opt::threads / 6) : filenames.size();
      numConsumers = numProducers*5;
    }
    size_t batchSize = 1024;
    
    
    // 使用FastxParser库解析FASTQ文件
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(filenames, numConsumers,numProducers, batchSize);
    parser.start();
    
    #pragma omp parallel num_threads(numConsumers)
    {
      // 每个线程获取一个读取组
      auto rg = parser.getReadGroup();
      
      while (true) {

        // 尝试填充读取组
        if (parser.refill(rg)) {
          // 处理读取组中的所有序列
          for (auto& read : rg) {
            insertCount(read.seq.c_str(), read.seq.length());
          }
          
          // 通知解析器我们已经完成了这个读取组的处理
          parser.finishedWithGroup(rg);
  
          // 在处理完一组后检查是否达到最大计数阈值
          if (m_maxCounts != 0 && m_totalCounts > m_maxCounts ) {
            parser.requestStop();
          }
          
        } else {
          // 没有更多数据可处理，退出循环
          break;
        }
      }
    } 
    // 停止解析器
    parser.stop();
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
        // outStr += "#@TK\t";
        // outStr += std::to_string(m_totalKmers);
        // outStr += "\n#@KS\t";
        // outStr += std::to_string(opt::k);
        // file << outStr;

        if (opt::information) {
          file << "#locusID\tmaxAT\tmaxCG\tcountAT\tcountCG\tsumAT\tsumCG\tdistinctAT\tdistinctCG\tref\tval\n";
        } else {
          file << "#locusID\tmaxAT\tmaxCG\tcountAT\tcountCG\tsumAT\tsumCG\tdistinctAT\tdistinctCG\n";
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
  // size_t>>, 10, std::mutex>;
  // Map1 m_counts;
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


  void initCountsHash() {
    m_alleleIDToKmerRef.resize(opt::snp.size());
    m_alleleIDToKmerVar.resize(opt::snp.size());
    m_alleleIDs.resize(opt::snp.size());
    for (unsigned i = 0; i < opt::snp.size(); ++i) {
      gzFile fp;
      fp = gzopen(opt::snp[i].c_str(), "r");
      if (fp == Z_NULL) {
          std::cerr << "file " << opt::snp[i] << " cannot be opened" << std::endl;
          exit(1);
      } else{
        std::cerr << "Reading " << opt::snp[i] << std::endl; 
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

      cerr << "alleleIDs_nums " << m_alleleIDs[i].size() << endl;
    }
  }
};
#endif /* SRC_FINGERPRINT_HPP_ */
