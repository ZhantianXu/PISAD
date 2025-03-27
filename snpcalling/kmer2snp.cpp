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

#include "parallel_hashmap/btree.h"
#include "parallel_hashmap/meminfo.h"
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_base.h"
#include "parallel_hashmap/phmap_bits.h"
#include "parallel_hashmap/phmap_config.h"
#include "parallel_hashmap/phmap_dump.h"
#include "parallel_hashmap/phmap_fwd_decl.h"
#include "parallel_hashmap/phmap_utils.h"

#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::no_property,
                              boost::property<boost::edge_weight_t, int>>
    Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;
typedef std::pair<int, int> EdgePair;

struct pair_hash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &p) const {
    std::size_t seed = 0;
    boost::hash_combine(seed, p.first);
    boost::hash_combine(seed, p.second);
    return seed;
  }
};

using Map1 = phmap::parallel_flat_hash_map<
    unsigned long long, char, std::hash<unsigned long long>,
    std::equal_to<unsigned long long>,
    std::allocator<std::pair<const unsigned long long, char>>, 6, std::mutex>;

using Map2 = phmap::parallel_flat_hash_map<
    unsigned long long, unsigned short, std::hash<unsigned long long>,
    std::equal_to<unsigned long long>,
    std::allocator<std::pair<const unsigned long long, unsigned short>>, 6,
    std::mutex>;

using Map3 = phmap::parallel_flat_hash_map<
    std::pair<unsigned long long, unsigned long long>,
    std::pair<std::string, std::string>, pair_hash,
    std::equal_to<std::pair<unsigned long long, unsigned long long>>,
    std::allocator<std::pair<const std::pair<unsigned long long, unsigned long long>,
                  std::pair<std::string, std::string>>>,
    6, std::mutex>;

// Function prototypes
int calc_weightThreshold(int k);
std::vector<std::tuple<unsigned long long, unsigned long long, char>>
snp_edges(const Map2 &m_index, int k, const Map1 &left_index,
          const Map1 &right_index, Map3 &extendKmers);
std::pair<std::string, std::string> extend_to_left(unsigned long long h1,
                                                   const Map1 &left_index,
                                                   const Map1 &right_index,
                                                   int k);
std::pair<std::string, std::string>
extend_to_right(const std::string &h1_binary, const Map1 &left_index,const Map1 &right_index, int k);
std::tuple<std::string, std::string, char, bool>
extend_one_pair(unsigned long long h1, unsigned long long h2,const Map1 &left_index, const Map1 &right_index, int k);
std::string reverse_binary_string(const std::string &s);
unsigned long long transfer_kmer_int(const std::string &s);
std::string transfer_int_kmer(unsigned long long val, int k);
std::string transfer_binary_string_2_kmer(const std::string &s);
unsigned long long reverse_int(unsigned long long v, int k);
std::string reverse_binary_string(const std::string &s);
void print_list(const std::vector<std::string> &l);
int file_lines(const std::string &filename);
void build_graph_and_find_components(const std::vector<std::tuple<unsigned long long, unsigned long long, char>>
        &edges,int k, std::string &name, Map3 &extendKmers, int lowCov, int highCov);
std::tuple<Map1, Map1, Map2> read(const std::string &input_filename, int low,int high, int k);

void print_resident_memory() {
  std::ifstream statm("/proc/self/statm");
  if (!statm.is_open()) {
    std::cerr << "Failed to open /proc/self/statm" << std::endl;
    return;
  }

  unsigned long size, resident, shared, text, lib, data, dt;
  statm >> size >> resident >> shared >> text >> lib >> data >> dt;
  statm.close();

  long page_size = sysconf(_SC_PAGE_SIZE);
  double resident_gb = static_cast<double>(resident) * page_size /(1024 * 1024 * 1024);
  std::cout << "Resident memory: " << resident_gb << " GB" << std::endl;
}
double getPeakRSSInGB() {
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
  // 将 kB 转换为 GB
  double peakRssGB = peakRss / (1024.0 * 1024.0);
  return peakRssGB;
}

std::tuple<unsigned long long, char> remove_bits(unsigned long long val, int k,int n) {
  // Calculate the mask to extract the lower k bits
  unsigned long long mask = (1ULL << k) - 1;
  unsigned long long binarykmer = val & mask;
  int x1;
  int x2;
  if (n == 0)
  {
    x1 = k / 2 - 1;
    x2 = k / 2;
  } else if (n == 1)
  {
    x1 = 0;
    x2 = 1;
  } else if (n == 2)
  {
    x1 = k - 2;
    x2 = k - 1;
  } else {
    std::cout << "error" << n << std::endl;
    exit(1);
  }

  // Extract the bits to be removed
  char removed_bits = (binarykmer >> (x1)) & 0b11;

  // Construct the new value by skipping the n-th and (n+1)-th bits
  unsigned long long new_val = 0;
  int bit_pos = 0;
  for (int i = 0; i < k; ++i) {
    if (i != (x1) && i != (x2)) {
      new_val |= ((binarykmer >> i) & 0b1) << bit_pos;
      ++bit_pos;
    }
  }

  return {new_val, removed_bits};
}
std::string get_filename(const std::string &path) {
  size_t pos = path.find_last_of('/');
  if (pos == std::string::npos) {
    return path;
  }
  return path.substr(pos + 1);
}

std::string ensure_trailing_slash(std::string path) {
  if (path.empty() || path.back() != '/') {
    path += '/';
  }
  return path;
}

int main(int argc, char **argv) {
  // Define default values
  std::string t1;            // Input file (required)
  int lowCov = 0;            // Low coverage (required)
  int highCov = 0;           // High coverage (required)
  int k = 21;                // Kmer size (optional, default: 21)
  int thread = 1;            // Number of threads (optional, default: 1)
  std::string output = "./"; // Output path (optional, default: ./)
  std::string name;          // Unused, kept empty

  // Define options
  static struct option long_options[] = {{"input", required_argument, 0, 'i'},
                                         {"left", required_argument, 0, 'l'},
                                         {"right", required_argument, 0, 'r'},
                                         {"kmer", required_argument, 0, 'k'},
                                         {"thread", required_argument, 0, 't'},
                                         {"output", required_argument, 0, 'o'},
                                         {"help", no_argument, 0, 'h'},
                                         {0, 0, 0, 0}};

  int opt;
  int option_index = 0;

  // Parse arguments
  while ((opt = getopt_long(argc, argv, "i:l:r:k:t:o:h", long_options,
                            &option_index)) != -1) {
    try {
      switch (opt) {
      case 'i':
        t1 = optarg;
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
        thread = std::stoi(optarg);
        break;
      case 'o':
        output = ensure_trailing_slash(optarg);
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
    } catch (const std::exception &e) {
      std::cerr << "Error parsing argument -" << (char)opt << ": " << e.what()
                << std::endl;
      return 1;
    }
  }

  // Check only required arguments
  if (t1.empty() || lowCov == 0 || highCov == 0) {
    std::cerr << "Error: Missing required arguments (-i, -l, -r)\n"
              << "Use " << argv[0] << " --help for usage information"
              << std::endl;
    return 1;
  }

  try {
    std::string filename = get_filename(t1);
    omp_set_num_threads(thread);
    Map3 extendKmers;

    auto time1 = std::chrono::high_resolution_clock::now();

    std::cout << "Heterozygous kmer coverage range: " << lowCov << " "
              << highCov << std::endl;

    auto [left_index, right_index, m_index] = read(t1, lowCov, highCov, k);

    auto time2 = std::chrono::high_resolution_clock::now();

    std::cout << "Finished reading (kmer cov) file, time cost: "
              << std::chrono::duration_cast<std::chrono::duration<double>>(time2 - time1).count()
              << " seconds" << std::endl;

    std::cout << "left_index count: " << left_index.size() << std::endl;
    std::cout << "right_index count: " << right_index.size() << std::endl;
    std::cout << "m_index count: " << m_index.size() << std::endl;

    std::vector<std::tuple<unsigned long long, unsigned long long, char>>
        edges = snp_edges(m_index, k, left_index, right_index, extendKmers);

    auto time3 = std::chrono::high_resolution_clock::now();
    std::cout << "Added SNP edges, time cost: "
              << std::chrono::duration_cast<std::chrono::duration<double>>(time3 - time2).count()
              << " seconds" << std::endl;
    std::cout << "Edge count: " << edges.size() << std::endl;

    auto time4 = std::chrono::high_resolution_clock::now();
    std::string output_filename = output + filename;
    build_graph_and_find_components(edges, k, output_filename, extendKmers,lowCov, highCov);

    auto time6 = std::chrono::high_resolution_clock::now();
    std::cout << "Built graph, time cost: "
              << std::chrono::duration_cast<std::chrono::duration<double>>(time6 - time4).count()
              << " seconds" << std::endl;

    double peakRssGB = getPeakRSSInGB();
    std::cout << "Peak Memory Usage: " << peakRssGB << " GB" << std::endl;
  } catch (const std::out_of_range &e) {
    std::cerr << "Error: std::out_of_range exception - " << e.what()<< std::endl;
    return 1;
  } catch (const std::invalid_argument &e) {
    std::cerr << "Error: std::invalid_argument exception - " << e.what()<< std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Error: Unknown exception" << std::endl;
    return 1;
  }

  return 0;
}

// Function definitions

int calc_weightThreshold(int k) {
  return (k - 3) / 2;
  // return heteRate >= 0.005 ? 4 : (k - 3) / 2;
}

int get_m_num(unsigned short value) {
  if (value == 0)
    return -1;
  int pos = 0;
  while (value) {
    pos++;
    value >>= 1;
  }
  return pos - 1;
}
unsigned long long highestOneBitPosition(unsigned long long value,
                                         unsigned long long value1, size_t x,
                                         int k) {
  unsigned long long low_bits =
      ((value1 >> (2 * x)) & 1) | (((value1 >> (2 * x + 1)) & 1) << 1);
  unsigned long long high_bits = (value >> (k - 1)) << 2;
  unsigned long long result =
      ((high_bits | low_bits) << (k - 1)) | (value & ((1ULL << (k - 1)) - 1));
  return result;
}

std::vector<std::tuple<unsigned long long, unsigned long long, char>>
snp_edges(const Map2 &m_index, int k, const Map1 &left_index,
          const Map1 &right_index, Map3 &extendKmers) {
  int num_threads = omp_get_max_threads();

  std::vector<
      std::vector<std::tuple<unsigned long long, unsigned long long, char>>>
      thread_local_edges(num_threads);

  // Calculate work distribution manually
  size_t total_keys = m_index.size();
  size_t chunk_size =
      (total_keys + num_threads - 1) / num_threads; // Ceiling division

#pragma omp parallel num_threads(num_threads)
  {
    int thread_id = omp_get_thread_num();
    auto &local_edges = thread_local_edges[thread_id];

    // Calculate start and end iterators for this thread
    auto it_start = m_index.begin();
    std::advance(it_start, thread_id * chunk_size);

    auto it_end = it_start;
    std::advance(it_end,std::min(chunk_size, total_keys - thread_id * chunk_size));

    // Iterate over the assigned range of keys
    for (auto it = it_start; it != it_end; ++it) {
      const auto &key_val = *it;
      const auto &key = key_val.first;
      const auto &m_key_len = get_m_num(key_val.second);
      if (m_key_len == 2)
        continue;

      for (size_t i = 0; i < m_key_len / 2 - 1; ++i) {
        for (size_t j = i + 1; j < m_key_len / 2; ++j) {
          unsigned long long k1 = highestOneBitPosition(key_val.first, key_val.second, i, k);
          unsigned long long k2 = highestOneBitPosition(key_val.first, key_val.second, j, k);
          auto [ek1, ek2, supportPair, flag] =
              extend_one_pair(k1, k2, left_index, right_index, k);
          if (flag) {
            if (k1 < k2) {
              local_edges.emplace_back(k1, k2, supportPair);
              extendKmers[{k1, k2}] = {ek1, ek2};
            } else {
              local_edges.emplace_back(k2, k1, supportPair);
              extendKmers[{k2, k1}] = {ek2, ek1};
            }
          }
        }
      }
    }
  }

  // Combine all local_edges from different threads
  std::vector<std::tuple<unsigned long long, unsigned long long, char>> edges;
  for (auto &local_edges : thread_local_edges) {
    edges.insert(edges.end(), std::make_move_iterator(local_edges.begin()),
                 std::make_move_iterator(local_edges.end()));
  }

  return edges;
}

unsigned long long extractBits(unsigned long long key, int k, int n) {
  unsigned long long result;
  if (n == 1) {
    // 取 key 的低 k 位
    result = key & ((1ULL << k) - 1);
  } else if (n == 2) {
    // 取 key 的高 k 位
    result = (key >> (sizeof(key) * 8 - k)) & ((1ULL << k) - 1);
  } else {
    // 如果 n 不是 1 或 2，返回 0 或者其他错误处理
    result = 0;
  }
  return result;
}
std::pair<std::string, std::string> extend_to_left(unsigned long long h1,
                                                   const Map1 &left_index,
                                                   const Map1 &right_index,
                                                   int k) {
  int mid = (k - 1) / 2;
  std::string h1_binary = std::bitset<64>(h1).to_string().substr(64 - k * 2);
  std::string temp = h1_binary, Rtemp = reverse_binary_string(h1_binary);

  unsigned long long key = std::stoull(temp.substr(0, 2 * k - 2), nullptr, 2);
  unsigned long long Rkey = std::stoull(Rtemp.substr(2), nullptr, 2);
  std::string add = "", Radd = "";

  for (int i = 0; i < mid; ++i) {
    bool flag = false, flagR = false;
    if (right_index.find(key) != right_index.end()) {
      temp = std::bitset<2>(right_index.at(key)).to_string() + temp;
      Rtemp = reverse_binary_string(temp);
      flag = true;
    }
    if (left_index.find(Rkey) != left_index.end()) {
      Rtemp = Rtemp + std::bitset<2>(left_index.at(Rkey)).to_string();
      temp = reverse_binary_string(Rtemp);
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
      add = std::bitset<2>(right_index.at(key)).to_string() + add;
      Radd = reverse_binary_string(add);
    } else if (!flag && flagR) {
      Radd = Radd + std::bitset<2>(left_index.at(Rkey)).to_string();
      add = reverse_binary_string(Radd);
    }
    key = std::stoull(temp.substr(0, 2 * k - 2), nullptr, 2);
    Rkey = std::stoull(Rtemp.substr(Rtemp.size() - (2 * k - 2)), nullptr, 2);
  }
  return {temp, add};
}

std::pair<std::string, std::string>
extend_to_right(const std::string &h1_binary, const Map1 &left_index,
                const Map1 &right_index, int k) {
  int mid = k / 2;
  std::string temp = h1_binary, Rtemp = reverse_binary_string(h1_binary);
  unsigned long long key = std::stoull(temp.substr(temp.size() - (2 * k - 2)), nullptr, 2);
  unsigned long long Rkey = std::stoull(Rtemp.substr(0, 2 * k - 2), nullptr, 2);
  std::string add = "", Radd = "";

  for (int i = 0; i < mid; ++i) {
    bool flag = false, flagR = false;
    if (left_index.find(key) != left_index.end()) {
      temp = temp + std::bitset<2>(left_index.at(key)).to_string();
      Rtemp = reverse_binary_string(temp);
      flag = true;
    }
    if (right_index.find(Rkey) != right_index.end()) {
      Rtemp = std::bitset<2>(right_index.at(Rkey)).to_string() + Rtemp;
      temp = reverse_binary_string(Rtemp);
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
      add = add + std::bitset<2>(left_index.at(key)).to_string();
      Radd = reverse_binary_string(add);
    } else if (!flag && flagR) {
      Radd = std::bitset<2>(right_index.at(Rkey)).to_string() + Radd;
      add = reverse_binary_string(Radd);
    }
    key = std::stoull(temp.substr(temp.size() - (2 * k - 2)), nullptr, 2);
    Rkey = std::stoull(Rtemp.substr(0, 2 * k - 2), nullptr, 2);
  }
  return {temp, add};
}

std::tuple<std::string, std::string, char, bool>
extend_one_pair(unsigned long long h1, unsigned long long h2,
                const Map1 &left_index, const Map1 &right_index, int k) {
  bool flag = true;
  char supportPairL = 0, supportPairR = 0;
  auto [temp1, add1] = extend_to_left(h1, left_index, right_index, k);
  auto [temp2, add2] = extend_to_left(h2, left_index, right_index, k);
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
  auto [ekmer1, add1R] = extend_to_right(temp1, left_index, right_index, k);
  auto [ekmer2, add2R] = extend_to_right(temp2, left_index, right_index, k);
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
  ekmer1 = transfer_binary_string_2_kmer(ekmer1);
  ekmer2 = transfer_binary_string_2_kmer(ekmer2);

  return {ekmer1, ekmer2, supportPairL + supportPairR, flag};
}

unsigned long long count_elements_in_group(const std::string &file_path,
                                           const std::string &group_path) {

  // 打开HDF5文件
  H5::H5File file(file_path, H5F_ACC_RDONLY);

  // 打开指定的组
  H5::Group group = file.openGroup(group_path);

  // 获取组中的对象数量
  hsize_t num_objs = group.getNumObjs();

  // 用于存储元素总数的变量
  unsigned long long total_count = 0;

  // 遍历组中的所有对象
  for (hsize_t i = 0; i < num_objs; ++i) {
    // 获取对象名称
    std::string obj_name = group.getObjnameByIdx(i);

    // 打开dataset
    H5::DataSet dataset = group.openDataSet(obj_name);

    // 获取dataset的数据空间
    H5::DataSpace dataspace = dataset.getSpace();

    // 获取数据集的尺寸
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, nullptr);

    // 累加dataset的元素数量
    total_count += dims[0];
  }

  return total_count;
}

std::tuple<Map1, Map1, Map2> read(const std::string &input_filename, int low,
                                  int high, int k) {
  Map1 left_index, right_index;
  Map2 m_index;
  using Set1 = phmap::parallel_flat_hash_set<
      unsigned long long, std::hash<unsigned long long>,
      std::equal_to<unsigned long long>, std::allocator<unsigned long long>, 6,
      std::mutex>;
  Set1 filterLKey, filterRKey;

  // unsigned long long number = count_elements_in_group(input_filename,
  // "/dsk/solid"); unsigned long long expected_elements = number / 2;

  std::ifstream infile(input_filename + ".histo");
  if (!infile) {
    std::cerr << "open error" << std::endl;
  }

  std::string line;
  int index;
  long long sum = 0;

  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    if (iss >> index) {
      long long value;
      if (iss >> value) {
        if (index >= low && index <= high) {
          sum += value;
        }
      }
    }
  }

  std::cout << "rserve_number: " << sum << std::endl;
  left_index.reserve(sum);
  right_index.reserve(sum);
  m_index.reserve(sum);

  const int k2 = k * 2;

  // 打开HDF5文件
  H5::H5File file(input_filename + ".h5", H5F_ACC_RDONLY);

  // 打开根组
  H5::Group root = file.openGroup("/dsk/solid");

  // 获取根组下的所有成员
  hsize_t num_objs = root.getNumObjs();

#pragma omp parallel for
  for (hsize_t i = 0; i < num_objs; ++i)
  // for (hsize_t i = 0; i < 20; ++i)meige
  {
    std::string datasetName = root.getObjnameByIdx(i);

    // 打开数据集
    H5::DataSet dataset = root.openDataSet(datasetName);

    // 获取数据空间
    H5::DataSpace dataspace = dataset.getSpace();

    // 获取数据集大小
    hsize_t dataSize;
    dataspace.getSimpleExtentDims(&dataSize, nullptr);

    // 创建缓冲区来存储数据
    struct KmerData {
      unsigned long long value;
      uint32_t abundance;
    };
    std::vector<KmerData> dataBuffer(dataSize);

    // 定义内存数据类型
    H5::CompType mtype(sizeof(KmerData));
    mtype.insertMember("value", HOFFSET(KmerData, value),H5::PredType::NATIVE_UINT64);
    mtype.insertMember("abundance", HOFFSET(KmerData, abundance),H5::PredType::NATIVE_UINT32);

    // 从数据集中读取数据
    dataset.read(dataBuffer.data(), mtype);

    for (const auto &kmer : dataBuffer) {
      unsigned long long val = kmer.value;
      uint32_t coverage = kmer.abundance;

      if (coverage < low || coverage > high)
        continue;

      unsigned long long newval = reverse_int(val, k2);
      if (val > newval)
        val = newval;
      auto [key, key_c] = remove_bits(val, k2, 0);

      //   m_index[key].emplace_back(val);
      auto it = m_index.find(key);
      if (it != m_index.end()) {
        if (get_m_num(it->second) < 14) {
          unsigned short new_value = (it->second << 2) + key_c;
          m_index[key] = new_value;
        }

      } else {

        m_index[key] = 4 + key_c;
      }

      auto [lkey, lkey_c] = remove_bits(val, k2, 1);
      if (!left_index.insert({lkey, lkey_c}).second)
        filterLKey.emplace(lkey);

      auto [rkey, rkey_c] = remove_bits(val, k2, 2);
      if (!right_index.insert({rkey, rkey_c}).second)
        filterRKey.emplace(rkey);
    }
  }
  std::cout << "Number of elements in filterLKey: " << filterLKey.size()
            << std::endl;
  std::cout << "Number of elements in filterRKey: " << filterRKey.size()
            << std::endl;

  for (const auto &key : filterLKey)
    left_index.erase(key);
  for (const auto &key : filterRKey)
    right_index.erase(key);

  return {std::move(left_index), std::move(right_index), std::move(m_index)};
}

std::unordered_map<std::string, std::string> RbinaryEle = {
    {"00", "10"}, {"01", "11"}, {"11", "01"}, {"10", "00"}};

std::string reverse_binary_string(const std::string &s) {
  std::string news;
  news.reserve(s.length());
  int lenS = s.length();
  if (lenS % 2 != 0) {
    std::cout << "binary string " << s << " is not correct len: " << lenS
              << std::endl;
    exit(1);
  }
  for (int i = 0; i < lenS; i += 2) {
    std::string pair = s.substr(i, 2);
    // news += RbinaryEle[pair];
    news = RbinaryEle[pair] + news;
  }
  return news;
}
unsigned long long transfer_kmer_int(const std::string &s) {
  unsigned long long result = 0;
  for (char c : s) {
    result <<= 2;
    switch (c) {
    case 'A':
      result |= 0b00;
      break;
    case 'C':
      result |= 0b01;
      break;
    case 'G':
      result |= 0b11;
      break;
    case 'T':
      result |= 0b10;
      break;
    }
  }
  return result;
}

std::string transfer_int_kmer(unsigned long long val, int k) {
  std::string s = std::bitset<64>(val).to_string().substr(64 - k * 2);
  return transfer_binary_string_2_kmer(s);
}

std::string transfer_binary_string_2_kmer(const std::string &s) {
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

unsigned long long reverse_int(unsigned long long var, int k) {
  // 取低k位
  unsigned long long mask = (1ULL << k) - 1;
  unsigned long long result = var & mask;

  // 按两位为一组进行顺序反转
  unsigned long long reversed = 0;
  int groupCount = k / 2;
  for (int i = 0; i < groupCount; ++i) {
    unsigned long long twoBits = (result >> (2 * i)) & 0x3; // 取第i组的两位

    // 按照 {00, 10}, {01, 11} 的规则进行替换
    switch (twoBits) {
    case 0b00:
      twoBits = 0b10;
      break;
    case 0b01:
      twoBits = 0b11;
      break;
    case 0b10:
      twoBits = 0b00;
      break;
    case 0b11:
      twoBits = 0b01;
      break;
    }

    reversed |= twoBits << (2 * (groupCount - 1 - i)); // 将第i组的两位放置到倒数第i组的位置
  }

  return reversed;
}

void print_list(const std::vector<std::string> &l) {
  for (const auto &ele : l)
    std::cout << ele << " ";
  std::cout << std::endl;
}

int file_lines(const std::string &filename) {
  std::ifstream infile(filename);
  return std::count(std::istreambuf_iterator<char>(infile),std::istreambuf_iterator<char>(), '\n');
}

void build_graph_and_find_components(
    const std::vector<std::tuple<unsigned long long, unsigned long long, char>>
        &edges,
    int k, std::string &name, Map3 &extendKmers, int lowCov, int highCov) {
  // Output SNP pairs to file
  std::string snpFile = name + "_" + std::to_string(k) + "_" +std::to_string(lowCov) + "_" + std::to_string(highCov) +"_pair.snp";
  std::ofstream snpOut(snpFile);
  if (!snpOut) {
    std::cerr << "Failed to open output file: " << snpFile << std::endl;
    return;
  }
  std::string exsnpFile = name + "_" + std::to_string(k) + "_" +std::to_string(lowCov) + "_" +
                          std::to_string(highCov) + "_pairex.snp";
  std::ofstream exsnpOut(exsnpFile);
  if (!exsnpOut) {
    std::cerr << "Failed to open output file: " << exsnpFile << std::endl;
    return;
  }

  Graph G;
  std::unordered_map<unsigned long long, Vertex> vertex_map;
  std::unordered_map<Vertex, unsigned long long> reverse_vertex_map;
  for (const auto &[u, v, weight] : edges) {
    if (vertex_map.find(u) == vertex_map.end()) {
      vertex_map[u] = boost::add_vertex(G);
      reverse_vertex_map[vertex_map[u]] = u;
    }
    if (vertex_map.find(v) == vertex_map.end()) {
      vertex_map[v] = boost::add_vertex(G);
      reverse_vertex_map[vertex_map[v]] = v;
    }
    boost::add_edge(vertex_map[u], vertex_map[v], weight, G);
  }

  std::vector<int> component(boost::num_vertices(G));
  int num = boost::connected_components(G, &component[0]);
  std::cout << "number of components: " << num << std::endl;

  std::vector<std::vector<Vertex>> components(num);
  for (auto v : boost::make_iterator_range(vertices(G))) {
    components[component[v]].push_back(v);
  }

  // Maximum weighted matching
  WeightMap weight_map = get(boost::edge_weight, G);
  std::vector<Vertex> mate(boost::num_vertices(G));
  bool success =
      boost::checked_edmonds_maximum_cardinality_matching(G, &mate[0]);
  if (!success) {
    std::cerr << "Error in matching algorithm." << std::endl;
    return;
  }
  std::cout << "size of matching: " << matching_size(G, &mate[0]) << std::endl;

  // Prepare to select mates and calculate weight distribution
  std::unordered_map<
      std::pair<unsigned long long, unsigned long long>, char,
      boost::hash<std::pair<unsigned long long, unsigned long long>>>
      selectedMates;
  std::unordered_map<int, int> weightDis;

  auto edge_pair = boost::edges(G);
  for (auto edge_iter = edge_pair.first; edge_iter != edge_pair.second;
       ++edge_iter) {
    auto e = *edge_iter;
    auto u = boost::source(e, G);
    auto v = boost::target(e, G);
    if (mate[u] == v && mate[v] == u) {
      int w = boost::get(weight_map, e);
      selectedMates[std::make_pair(reverse_vertex_map[u],reverse_vertex_map[v])] = w;
      if (weightDis.find(w) == weightDis.end()) {
        weightDis[w] = 0;
      }
      weightDis[w] += 1;
    }
  }

  // Sort weight distribution
  std::vector<std::pair<int, int>> sortedWeightDis(weightDis.begin(),weightDis.end());
  std::sort(sortedWeightDis.begin(), sortedWeightDis.end());

  int weightThreshold = calc_weightThreshold(k);
  std::cout << "selected kmer pairs " << selectedMates.size() << std::endl;
  std::cout << "weightThreshold " << weightThreshold << std::endl;

  // Processing selected mates and SNP pairs
  int numid = 0;
  for (const auto &[e1_e2, w] : selectedMates) {
    auto e1 = e1_e2.first;
    auto e2 = e1_e2.second;
    std::string ele1 = transfer_int_kmer(e1, k);
    std::string ele2 = transfer_int_kmer(e2, k);
    if (w < weightThreshold) {
      numid++;
      continue;
    }

    snpOut << ele1 << " " << ele2 << " " << static_cast<int>(w) << std::endl;
    auto [ek1, ek2] = extendKmers[{e1, e2}];
    exsnpOut << ek1 << " " << ek2 << " " << static_cast<int>(w) << std::endl;
  }

  std::cout << "delete " << numid << " kmer pairs" << std::endl;

  snpOut.close();
  exsnpOut.close();
}