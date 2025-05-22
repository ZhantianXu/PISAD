#include <H5Cpp.h>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <omp.h>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "parallel_hashmap/phmap.h"

#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

// Hash function for pairs
struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, p.second);
        return seed;
    }
};

class spinlock_mutex {
  // https://rigtorp.se/spinlock/
  // https://vorbrodt.blog/2019/02/12/fast-mutex/
  public:
     // Assignment is disabled.
     spinlock_mutex& operator=(const spinlock_mutex& rhs) = delete;

     void lock() noexcept {
        for (;;) {
           if (!lock_.exchange(true, std::memory_order_acquire))
              break;
           while (lock_.load(std::memory_order_relaxed))
              __builtin_ia32_pause();
        }
     }

     void unlock() noexcept {
        lock_.store(false, std::memory_order_release);
     }

  private:
     alignas(4 * sizeof(std::max_align_t)) std::atomic_bool lock_ = false;
  };

// Type aliases for clarity
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::no_property,
                              boost::property<boost::edge_weight_t, int>> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;
typedef std::pair<int, int> EdgePair;

// Map type definitions using parallel hash maps
using KmerCharMap = phmap::parallel_flat_hash_map<
    unsigned long long, char, std::hash<unsigned long long>,
    std::equal_to<unsigned long long>,
    std::allocator<std::pair<const unsigned long long, char>>, 10, spinlock_mutex>;

using KmerCountMap = phmap::parallel_flat_hash_map<
    unsigned long long, unsigned short, std::hash<unsigned long long>,
    std::equal_to<unsigned long long>,
    std::allocator<std::pair<const unsigned long long, unsigned short>>, 10,
    spinlock_mutex>;

using KmerPairStringMap = phmap::parallel_flat_hash_map<
    std::pair<unsigned long long, unsigned long long>,
    std::pair<std::string, std::string>, PairHash,
    std::equal_to<std::pair<unsigned long long, unsigned long long>>,
    std::allocator<std::pair<const std::pair<unsigned long long, unsigned long long>,
    std::pair<std::string, std::string>>>, 10, spinlock_mutex>;

// Set type definition
using KmerSet = phmap::parallel_flat_hash_set<
    unsigned long long, std::hash<unsigned long long>,
    std::equal_to<unsigned long long>, std::allocator<unsigned long long>, 10,
    spinlock_mutex>;

// Utility class for sequence manipulation
class SequenceUtils {
public:
    // Static member initialization
    static constexpr std::array<const char*, 4> RbinaryEle = {
      "10", "11", "00", "01"
  };

  static std::string reverseComplement(const std::string& s) {
      int lenS = s.length();
      if (lenS % 2 != 0) {
          std::cerr << "Binary string " << s << " has incorrect length: " << lenS << std::endl;
          std::exit(1);
      }

      std::string result(lenS, '0');

      for (int i = 0; i < lenS; i += 2) {
          // 二进制字符串转换为整数索引
          int index = ((s[i] - '0') << 1) | (s[i + 1] - '0');
          const char* complement = RbinaryEle[index];

          result[lenS - i - 2] = complement[0];
          result[lenS - i - 1] = complement[1];
      }

      return result;
  }
    
    static unsigned long long kmerToInt(const std::string& s) {
        unsigned long long result = 0;
        for (char c : s) {
            result <<= 2;
            switch (c) {
                case 'A': result |= 0b00; break;
                case 'C': result |= 0b01; break;
                case 'G': result |= 0b11; break;
                case 'T': result |= 0b10; break;
            }
        }
        return result;
    }
    
    static std::string intToKmer(unsigned long long val, int k) {
        std::string s = std::bitset<64>(val).to_string().substr(64 - k * 2);
        return binaryStringToKmer(s);
    }
    
    static std::string binaryStringToKmer(const std::string& s) {
        std::string kmer;
        for (size_t i = 0; i < s.size(); i += 2) {
            std::string sub = s.substr(i, 2);
            if (sub == "00")
                kmer += 'A';
            else if (sub == "01")
                kmer += 'C';
            else if (sub == "11")
                kmer += 'G';
            else if (sub == "10")
                kmer += 'T';
        }
        return kmer;
    }
    
    static unsigned long long reverseInt(unsigned long long var, int k) {
        unsigned long long mask = (1ULL << k) - 1;
        unsigned long long result = var & mask;
        
        unsigned long long reversed = 0;
        int groupCount = k / 2;
        for (int i = 0; i < groupCount; ++i) {
            unsigned long long twoBits = (result >> (2 * i)) & 0x3;
            
            switch (twoBits) {
                case 0b00: twoBits = 0b10; break;
                case 0b01: twoBits = 0b11; break;
                case 0b10: twoBits = 0b00; break;
                case 0b11: twoBits = 0b01; break;
            }
            
            reversed |= twoBits << (2 * (groupCount - 1 - i));
        }
        
        return reversed;
    }
    
    static std::tuple<unsigned long long, char> removeBits(unsigned long long val, int k, int n) {
        unsigned long long mask = (1ULL << k) - 1;
        unsigned long long binarykmer = val & mask;
        int x1, x2;
        
        if (n == 0) {
            x1 = k / 2 - 1;
            x2 = k / 2;
        } else if (n == 1) {
            x1 = 0;
            x2 = 1;
        } else if (n == 2) {
            x1 = k - 2;
            x2 = k - 1;
        } else {
            std::cerr << "Error: Invalid n value: " << n << std::endl;
            exit(1);
        }
        
        char removedBits = (binarykmer >> (x1)) & 0b11;
        
        unsigned long long newVal = 0;
        int bitPos = 0;
        for (int i = 0; i < k; ++i) {
            if (i != (x1) && i != (x2)) {
                newVal |= ((binarykmer >> i) & 0b1) << bitPos;
                ++bitPos;
            }
        }
        
        return {newVal, removedBits};
    }
    
    static int getHighestBitPosition(unsigned short value) {
        if (value == 0)
            return -1;
        int pos = 0;
        while (value) {
            pos++;
            value >>= 1;
        }
        return pos - 1;
    }
    
    static unsigned long long computeHighestOneBitPosition(unsigned long long value,
                                                unsigned long long value1, size_t x, int k) {
        unsigned long long lowBits =
            ((value1 >> (2 * x)) & 1) | (((value1 >> (2 * x + 1)) & 1) << 1);
        unsigned long long highBits = (value >> (k - 1)) << 2;
        unsigned long long result =
            ((highBits | lowBits) << (k - 1)) | (value & ((1ULL << (k - 1)) - 1));
        return result;
    }
    
    static std::string getFilename(const std::string& path) {
        size_t pos = path.find_last_of('/');
        return (pos == std::string::npos) ? path : path.substr(pos + 1);
    }
    
    static std::string ensureTrailingSlash(std::string path) {
        if (path.empty() || path.back() != '/') {
            path += '/';
        }
        return path;
    }
    
    static double getPeakRSSInGB() {
        std::ifstream file("/proc/self/status");
        std::string line;
        size_t peakRss = 0;
        while (std::getline(file, line)) {
            if (line.find("VmPeak:") == 0) {
                std::istringstream iss(line);
                std::string key;
                size_t value;
                std::string unit;
                iss >> key >> value >> unit;
                if (unit == "kB") {
                    peakRss = value;
                }
                break;
            }
        }
        return peakRss / (1024.0 * 1024.0); // Convert kB to GB
    }
};

// Class for handling k-mer processing
class KmerProcessor {
private:
    int k;
    std::string inputFilename;  
    int lowCov;                 
    int highCov;                
    std::vector<std::tuple<unsigned long long, unsigned long long, char>> edges;  
    KmerCharMap leftIndex, rightIndex;
    KmerCountMap mIndex;
    KmerPairStringMap extendKmers;

public:
    KmerProcessor(int k, const std::string& inputFilename, int lowCov, int highCov) 
        : k(k), inputFilename(inputFilename), lowCov(lowCov), highCov(highCov) {}
    
    void readKmersFromHDF5() {
        // Read histogram to estimate size
        std::ifstream infile(inputFilename + ".histo");
        if (!infile) {
            throw std::runtime_error("Failed to open histogram file");
        }
        
        std::string line;
        int index;
        long long sum = 0;
        
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (iss >> index) {
                long long value;
                if (iss >> value && index >= lowCov && index <= highCov) {
                    sum += value;
                }
            }
        }
        
        std::cout << "Reserve number: " << sum << std::endl;
        leftIndex.reserve(sum);
        rightIndex.reserve(sum);
        mIndex.reserve(sum);
        
        const int k2 = k * 2;
        KmerSet filterLKey, filterRKey;
        
        try {
            // Open HDF5 file
            H5::H5File file(inputFilename + ".h5", H5F_ACC_RDONLY);
            H5::Group root = file.openGroup("/dsk/solid");
            hsize_t numObjs = root.getNumObjs();
            
            #pragma omp parallel for schedule(dynamic)
            for (hsize_t i = 0; i < numObjs; ++i) {
                std::string datasetName = root.getObjnameByIdx(i);
                H5::DataSet dataset = root.openDataSet(datasetName);
                H5::DataSpace dataspace = dataset.getSpace();
                
                hsize_t dataSize;
                dataspace.getSimpleExtentDims(&dataSize, nullptr);
                
                struct KmerData {
                    unsigned long long value;
                    uint32_t abundance;
                };
                std::vector<KmerData> dataBuffer(dataSize);
                
                H5::CompType mtype(sizeof(KmerData));
                mtype.insertMember("value", HOFFSET(KmerData, value), H5::PredType::NATIVE_UINT64);
                mtype.insertMember("abundance", HOFFSET(KmerData, abundance), H5::PredType::NATIVE_UINT32);
                
                dataset.read(dataBuffer.data(), mtype);
                
                for (const auto& kmer : dataBuffer) {
                    unsigned long long val = kmer.value;
                    uint32_t coverage = kmer.abundance;
                    
                    if (coverage < lowCov || coverage > highCov)
                        continue;
                    
                    unsigned long long newval = SequenceUtils::reverseInt(val, k2);
                    if (val > newval)
                        val = newval;
                    
                    // 使用SequenceUtils中的函数
                    auto [key, keyC] = SequenceUtils::removeBits(val, k2, 0);
                    
                    mIndex.try_emplace_l(
                      key,
                      [&](auto& pair) {
                          if (SequenceUtils::getHighestBitPosition(pair.second) < 14) {
                              pair.second = (pair.second << 2) + keyC;
                          }
                      },
                      static_cast<unsigned short>(4 + keyC)
                  );
                  
                    
                    auto [lkey, lkeyC] = SequenceUtils::removeBits(val, k2, 1);
                    if (!leftIndex.try_emplace(lkey, lkeyC).second)
                        filterLKey.emplace(lkey);
                    
                    auto [rkey, rkeyC] = SequenceUtils::removeBits(val, k2, 2);
                    if (!rightIndex.try_emplace(rkey, rkeyC).second)
                        filterRKey.emplace(rkey);
                }
            }
            
            std::cout << "Number of elements in filterLKey: " << filterLKey.size() << std::endl;
            std::cout << "Number of elements in filterRKey: " << filterRKey.size() << std::endl;
            
            for (const auto& key : filterLKey)
                leftIndex.erase(key);
            for (const auto& key : filterRKey)
                rightIndex.erase(key);
                
            std::cout << "left_index count: " << leftIndex.size() << std::endl;
            std::cout << "right_index count: " << rightIndex.size() << std::endl;
            std::cout << "m_index count: " << mIndex.size() << std::endl;
            
        } catch (const H5::Exception& e) {
            std::cerr << "HDF5 Error: " << e.getDetailMsg() << std::endl;
            throw;
        }
    }
    
    std::pair<std::string, std::string> extendToLeft(unsigned long long h1) {
        int mid = (k - 1) / 2;
        std::string h1Binary = std::bitset<64>(h1).to_string().substr(64 - k * 2);
        std::string temp = h1Binary, Rtemp = SequenceUtils::reverseComplement(h1Binary);
        
        unsigned long long key = std::stoull(temp.substr(0, 2 * k - 2), nullptr, 2);
        unsigned long long Rkey = std::stoull(Rtemp.substr(2), nullptr, 2);
        std::string add = "", Radd = "";
        
        for (int i = 0; i < mid; ++i) {
            bool flag = false, flagR = false;
            if (rightIndex.find(key) != rightIndex.end()) {
                temp = std::bitset<2>(rightIndex.at(key)).to_string() + temp;
                Rtemp = SequenceUtils::reverseComplement(temp);
                flag = true;
            }
            if (leftIndex.find(Rkey) != leftIndex.end()) {
                Rtemp = Rtemp + std::bitset<2>(leftIndex.at(Rkey)).to_string();
                temp = SequenceUtils::reverseComplement(Rtemp);
                flagR = true;
            }
            if (flag == flagR) {
                if (flag) {
                    temp = temp.substr(2);
                    Rtemp = Rtemp.substr(0, Rtemp.size() - 2);
                }
                break;
            }
            if (flag && !flagR) {
                add = std::bitset<2>(rightIndex.at(key)).to_string() + add;
                Radd = SequenceUtils::reverseComplement(add);
            } else if (!flag && flagR) {
                Radd = Radd + std::bitset<2>(leftIndex.at(Rkey)).to_string();
                add = SequenceUtils::reverseComplement(Radd);
            }
            key = std::stoull(temp.substr(0, 2 * k - 2), nullptr, 2);
            Rkey = std::stoull(Rtemp.substr(Rtemp.size() - (2 * k - 2)), nullptr, 2);
        }
        return {temp, add};
    }
    
    std::pair<std::string, std::string> extendToRight(const std::string& h1Binary) {
        int mid = k / 2;
        std::string temp = h1Binary, Rtemp = SequenceUtils::reverseComplement(h1Binary);
        unsigned long long key = std::stoull(temp.substr(temp.size() - (2 * k - 2)), nullptr, 2);
        unsigned long long Rkey = std::stoull(Rtemp.substr(0, 2 * k - 2), nullptr, 2);
        std::string add = "", Radd = "";
        
        for (int i = 0; i < mid; ++i) {
            bool flag = false, flagR = false;
            if (leftIndex.find(key) != leftIndex.end()) {
                temp = temp + std::bitset<2>(leftIndex.at(key)).to_string();
                Rtemp = SequenceUtils::reverseComplement(temp);
                flag = true;
            }
            if (rightIndex.find(Rkey) != rightIndex.end()) {
                Rtemp = std::bitset<2>(rightIndex.at(Rkey)).to_string() + Rtemp;
                temp = SequenceUtils::reverseComplement(Rtemp);
                flagR = true;
            }
            if (flag == flagR) {
                if (flag) {
                    temp = temp.substr(0, temp.size() - 2);
                    Rtemp = Rtemp.substr(2);
                }
                break;
            }
            if (flag && !flagR) {
                add = add + std::bitset<2>(leftIndex.at(key)).to_string();
                Radd = SequenceUtils::reverseComplement(add);
            } else if (!flag && flagR) {
                Radd = std::bitset<2>(rightIndex.at(Rkey)).to_string() + Radd;
                add = SequenceUtils::reverseComplement(Radd);
            }
            key = std::stoull(temp.substr(temp.size() - (2 * k - 2)), nullptr, 2);
            Rkey = std::stoull(Rtemp.substr(0, 2 * k - 2), nullptr, 2);
        }
        return {temp, add};
    }
    
    std::tuple<std::string, std::string, char, bool> extendOnePair(unsigned long long h1, unsigned long long h2) {
        bool flag = true;
        char supportPairL = 0, supportPairR = 0;
        
        auto [temp1, add1] = extendToLeft(h1);
        auto [temp2, add2] = extendToLeft(h2);
        int minL = std::min(add1.size(), add2.size());
        add1 = add1.substr(add1.size() - minL);
        add2 = add2.substr(add2.size() - minL);
        
        for (int i = minL; i > 0; i -= 2) {
            if (add1.substr(i - 2, 2) == add2.substr(i - 2, 2))
                supportPairL++;
            else
                break;
        }
        
        temp1 = temp1.substr(temp1.size() - (supportPairL + k) * 2);
        temp2 = temp2.substr(temp2.size() - (supportPairL + k) * 2);
        auto [ekmer1, add1R] = extendToRight(temp1);
        auto [ekmer2, add2R] = extendToRight(temp2);
        int minR = std::min(add1R.size(), add2R.size());
        
        for (int i = 0; i < minR; i += 2) {
            if (add1R.substr(i, 2) == add2R.substr(i, 2))
                supportPairR++;
            else
                break;
        }
        
        int len = (supportPairL + k + supportPairR) * 2;
        ekmer1 = ekmer1.substr(0, len);
        ekmer2 = ekmer2.substr(0, len);
        ekmer1 = SequenceUtils::binaryStringToKmer(ekmer1);
        ekmer2 = SequenceUtils::binaryStringToKmer(ekmer2);
        
        return {ekmer1, ekmer2, supportPairL + supportPairR, flag};
    }
    
    void findSnpEdges() {
        int numThreads = omp_get_max_threads();
        
        std::vector<std::vector<std::tuple<unsigned long long, unsigned long long, char>>> threadLocalEdges(numThreads);
        
        size_t totalKeys = mIndex.size();
        size_t chunkSize = (totalKeys + numThreads - 1) / numThreads;
        
        #pragma omp parallel num_threads(numThreads)
        {
            int threadId = omp_get_thread_num();
            auto& localEdges = threadLocalEdges[threadId];
            
            auto itStart = mIndex.begin();
            std::advance(itStart, threadId * chunkSize);
            
            auto itEnd = itStart;
            std::advance(itEnd, std::min(chunkSize, totalKeys - threadId * chunkSize));
            
            for (auto it = itStart; it != itEnd; ++it) {
                const auto& keyVal = *it;
                const auto& key = keyVal.first;
                // 使用SequenceUtils中的函数
                const auto& mKeyLen = SequenceUtils::getHighestBitPosition(keyVal.second);
                
                if (mKeyLen == 2)
                    continue;
                
                for (size_t i = 0; i < mKeyLen / 2 - 1; ++i) {
                    for (size_t j = i + 1; j < mKeyLen / 2; ++j) {
                        // 使用SequenceUtils中的函数
                        unsigned long long k1 = SequenceUtils::computeHighestOneBitPosition(key, keyVal.second, i, k);
                        unsigned long long k2 = SequenceUtils::computeHighestOneBitPosition(key, keyVal.second, j, k);
                        
                        auto [ek1, ek2, supportPair, flag] = extendOnePair(k1, k2);
                        
                        if (flag) {
                            if (k1 < k2) {
                                localEdges.emplace_back(k1, k2, supportPair);
                                extendKmers[{k1, k2}] = {ek1, ek2};
                            } else {
                                localEdges.emplace_back(k2, k1, supportPair);
                                extendKmers[{k2, k1}] = {ek2, ek1};
                            }
                        }
                    }
                }
            }
        }
        
        edges.clear();
        for (auto& localEdges : threadLocalEdges) {
            edges.insert(edges.end(), std::make_move_iterator(localEdges.begin()),
                        std::make_move_iterator(localEdges.end()));
        }
    }
    
    int calculateWeightThreshold() {
        return (k - 3) / 2;
    }
    
    void buildGraphAndFindComponents(const std::string& outputPrefix) {
        
        // Create output files
        std::string snpFile = outputPrefix + "_" + std::to_string(k) + "_" +
                              std::to_string(lowCov) + "_" + std::to_string(highCov) + "_pair.snp";
        std::ofstream snpOut(snpFile);
        if (!snpOut) {
            throw std::runtime_error("Failed to open output file: " + snpFile);
        }
        
        std::string exsnpFile = outputPrefix + "_" + std::to_string(k) + "_" +
                              std::to_string(lowCov) + "_" + std::to_string(highCov) + "_pairex.snp";
        std::ofstream exsnpOut(exsnpFile);
        if (!exsnpOut) {
            throw std::runtime_error("Failed to open output file: " + exsnpFile);
        }
        
        // Build the graph
        Graph G;
        std::unordered_map<unsigned long long, Vertex> vertexMap;
        std::unordered_map<Vertex, unsigned long long> reverseVertexMap;
        
        for (const auto& [u, v, weight] : edges) {
            if (vertexMap.find(u) == vertexMap.end()) {
                vertexMap[u] = boost::add_vertex(G);
                reverseVertexMap[vertexMap[u]] = u;
            }
            if (vertexMap.find(v) == vertexMap.end()) {
                vertexMap[v] = boost::add_vertex(G);
                reverseVertexMap[vertexMap[v]] = v;
            }
            boost::add_edge(vertexMap[u], vertexMap[v], weight, G);
        }
        
        // Find connected components
        std::vector<int> component(boost::num_vertices(G));
        int numComponents = boost::connected_components(G, &component[0]);
        std::cout << "Number of components: " << numComponents << std::endl;
        
        std::vector<std::vector<Vertex>> components(numComponents);
        for (auto v : boost::make_iterator_range(vertices(G))) {
            components[component[v]].push_back(v);
        }
        
        // Maximum weighted matching
        WeightMap weightMap = get(boost::edge_weight, G);
        std::vector<Vertex> mate(boost::num_vertices(G));
        bool success = boost::checked_edmonds_maximum_cardinality_matching(G, &mate[0]);
        
        if (!success) {
            throw std::runtime_error("Error in matching algorithm.");
        }
        
        std::cout << "Size of matching: " << matching_size(G, &mate[0]) << std::endl;
        
        // Select mates and calculate weight distribution
        std::unordered_map<std::pair<unsigned long long, unsigned long long>, char,
            boost::hash<std::pair<unsigned long long, unsigned long long>>> selectedMates;
        
        std::unordered_map<int, int> weightDis;
        
        auto edgePair = boost::edges(G);
        for (auto edgeIter = edgePair.first; edgeIter != edgePair.second; ++edgeIter) {
            auto e = *edgeIter;
            auto u = boost::source(e, G);
            auto v = boost::target(e, G);
            
            if (mate[u] == v && mate[v] == u) {
                int w = boost::get(weightMap, e);
                selectedMates[std::make_pair(reverseVertexMap[u], reverseVertexMap[v])] = w;
                
                weightDis[w] = weightDis.find(w) == weightDis.end() ? 1 : weightDis[w] + 1;
            }
        }
        
        // Sort weight distribution
        std::vector<std::pair<int, int>> sortedWeightDis(weightDis.begin(), weightDis.end());
        std::sort(sortedWeightDis.begin(), sortedWeightDis.end());
        
        int weightThreshold = calculateWeightThreshold();
        std::cout << "Selected kmer pairs: " << selectedMates.size() << std::endl;
        std::cout << "Weight threshold: " << weightThreshold << std::endl;
        
        // Process selected mates and write SNP pairs to files
        int numDeleted = 0;
        for (const auto& [e1_e2, w] : selectedMates) {
            auto e1 = e1_e2.first;
            auto e2 = e1_e2.second;
            
            std::string ele1 = SequenceUtils::intToKmer(e1, k);
            std::string ele2 = SequenceUtils::intToKmer(e2, k);
            
            if (w < weightThreshold) {
                numDeleted++;
                continue;
            }
            
            snpOut << ele1 << " " << ele2 << " " << static_cast<int>(w) << std::endl;
            
            auto [ek1, ek2] = extendKmers[{e1, e2}];
            exsnpOut << ek1 << " " << ek2 << " " << static_cast<int>(w) << std::endl;
        }
        
        std::cout << "Deleted " << numDeleted << " kmer pairs" << std::endl;
        
        snpOut.close();
        exsnpOut.close();
    }
    
    KmerCharMap& getLeftIndex() { return leftIndex; }
    KmerCharMap& getRightIndex() { return rightIndex; }
    KmerCountMap& getMIndex() { return mIndex; }
    KmerPairStringMap& getExtendKmers() { return extendKmers; }
    int getK() const { return k; }
};

// Main function
int main(int argc, char **argv) {
    // Default parameter values
    std::string inputFile;        // Input file (required)
    int lowCov = 0;               // Low coverage (required)
    int highCov = 0;              // High coverage (required)
    int k = 21;                   // Kmer size (optional, default: 21)
    int threads = 1;              // Number of threads (optional, default: 1)
    std::string outputPath = "./"; // Output path (optional, default: ./)
    
    // Define command-line options
    static struct option longOptions[] = {
        {"input", required_argument, 0, 'i'},
        {"left", required_argument, 0, 'l'},
        {"right", required_argument, 0, 'r'},
        {"kmer", required_argument, 0, 'k'},
        {"thread", required_argument, 0, 't'},
        {"output", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    // Parse command-line arguments
    int opt;
    int optionIndex = 0;
    
    while ((opt = getopt_long(argc, argv, "i:l:r:k:t:o:h", longOptions, &optionIndex)) != -1) {
        try {
            switch (opt) {
                case 'i':
                    inputFile = optarg;
                    break;
                case 'l':
                    lowCov = std::stoi(optarg);
                    break;
                case 'r':
                    highCov = std::stoi(optarg);
                    break;
                case 'k':
                    k = std::stoi(optarg);
                    break;
                case 't':
                    threads = std::stoi(optarg);
                    break;
                case 'o':
                    outputPath = SequenceUtils::ensureTrailingSlash(optarg);
                    break;
                case 'h':
                    std::cerr << "Usage: " << argv[0] << " [options]\n"
                              << "Required options:\n"
                              << "  -i, --input <file>    Kmer frequency file\n"
                              << "  -l, --left <int>      Low coverage threshold\n"
                              << "  -r, --right <int>     High coverage threshold\n"
                              << "Optional options:\n"
                              << "  -k, --kmer <int>      Kmer size (default: 21)\n"
                              << "  -t, --thread <int>    Number of threads (default: 1)\n"
                              << "  -o, --output <path>   Output path (default: ./)\n"
                              << "  -h, --help            Show this help message\n"
                              << std::endl;
                    return 1;
                default:
                    return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error parsing argument -" << (char)opt << ": " << e.what() << std::endl;
            return 1;
        }
    }
    
    // Validate required arguments
    if (inputFile.empty() || lowCov == 0 || highCov == 0) {
        std::cerr << "Error: Missing required arguments (-i, -l, -r)\n"
                  << "Use " << argv[0] << " --help for usage information" << std::endl;
        return 1;
    }
    
    try {
        std::string filename = SequenceUtils::getFilename(inputFile);
        omp_set_num_threads(threads);
        
        auto startTime = std::chrono::high_resolution_clock::now();
        std::cout << "Heterozygous kmer coverage range: " << lowCov << " " << highCov << std::endl;
        
        // Create KmerProcessor instance and process data
        KmerProcessor processor(k, inputFile, lowCov, highCov);
        
        // Read k-mers from HDF5 file
        processor.readKmersFromHDF5();
        
        auto readTime = std::chrono::high_resolution_clock::now();
        std::cout << "Finished reading (kmer cov) file, time cost: "
                  << std::chrono::duration_cast<std::chrono::duration<double>>(readTime - startTime).count()
                  << " seconds" << std::endl;
        
        // Find SNP edges
        processor.findSnpEdges();
        
        auto edgesTime = std::chrono::high_resolution_clock::now();
        std::cout << "Added SNP edges, time cost: "
                  << std::chrono::duration_cast<std::chrono::duration<double>>(edgesTime - readTime).count()
                  << " seconds" << std::endl;
        
        // Build graph and find components
        std::string outputFilename = outputPath + filename;
        processor.buildGraphAndFindComponents(outputFilename);
        
        auto graphTime = std::chrono::high_resolution_clock::now();
        std::cout << "Built graph, time cost: "
                  << std::chrono::duration_cast<std::chrono::duration<double>>(graphTime - edgesTime).count()
                  << " seconds" << std::endl;
        
        // Report memory usage
        double peakRssGB = SequenceUtils::getPeakRSSInGB();
        std::cout << "Peak Memory Usage: " << peakRssGB << " GB" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Error: Unknown exception" << std::endl;
        return 1;
    }
    
    return 0;
}