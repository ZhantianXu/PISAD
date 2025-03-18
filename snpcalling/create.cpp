#include "parallel_hashmap/btree.h"
#include "parallel_hashmap/meminfo.h"
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_base.h"
#include "parallel_hashmap/phmap_bits.h"
#include "parallel_hashmap/phmap_config.h"
#include "parallel_hashmap/phmap_dump.h"
#include "parallel_hashmap/phmap_fwd_decl.h"
#include "parallel_hashmap/phmap_utils.h"
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
#include <malloc.h>
#include <omp.h>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using Map1 = phmap::parallel_flat_hash_set<
    unsigned long long, std::hash<unsigned long long>,
    std::equal_to<unsigned long long>, std::allocator<unsigned long long>, 6,
    std::mutex>;
using Map2 = phmap::parallel_flat_hash_map<
    unsigned long long, char, std::hash<unsigned long long>,
    std::equal_to<unsigned long long>,
    std::allocator<std::pair<const unsigned long long, char>>, 6, std::mutex>;
Map2 kmer_count;
unsigned long long transfer_kmer_int(const std::string &s);
unsigned long long reverse_int(unsigned long long v, int k);
void read(const std::string &input_filename, int low, int high, int k);
std::vector<std::string> splitByN(const std::string &sequence);
void processFile(const std::string &inputFile, const std::string &outputFile,int k, int limlt);
std::string transfer_int_kmer(unsigned long long val, int k);
std::string transfer_binary_string_2_kmer(const std::string &s);

std::string get_filename(const std::string &path) {
  // First find the last '/'
  size_t slash_pos = path.find_last_of('/');
  size_t start_pos = (slash_pos == std::string::npos) ? 0 : slash_pos + 1;

  // Then find the last '.snp'
  size_t snp_pos = path.rfind("_pairex.snp");

  // If '.snp' not found, return from last '/' to end
  if (snp_pos == std::string::npos || snp_pos < start_pos) {
    return path.substr(start_pos);
  }

  // Return substring between last '/' and last '.snp'
  return path.substr(start_pos, snp_pos - start_pos);
}

std::string ensure_trailing_slash(std::string path) {
  if (path.empty() || path.back() != '/') {
    path += '/';
  }
  return path;
}

int main(int argc, char **argv) {
  std::string inputFile, filename;
  std::string output; // 设置默认输出文件名
  int k = 21;         // 设置默认k值
  int limlt = k;      // 设置默认limit值
  bool inputSpecified = false;

  // 选项解析
  int opt;
  while ((opt = getopt(argc, argv, "i:o:k:l:")) != -1) {
    switch (opt) {
    case 'i':
      inputFile = optarg;
      filename = get_filename(inputFile);
      inputSpecified = true;
      break;
    case 'o':
      output = ensure_trailing_slash(optarg);
      break;
    case 'k':
      k = std::stoi(optarg);
      break;
    case 'l':
      limlt = std::stoi(optarg);
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " -i <inputFile> [-o <outputFile>] [-k <k>] [-l <limit>]\n";
      return 1;
    }
  }

  // 检查是否提供了必须的输入文件参数
  if (!inputSpecified) {
    std::cerr << "Error: Input file (-i) is required\n";
    std::cerr << "Usage: " << argv[0]
              << " -i <inputFile> [-o <outputFile>] [-k <k>] [-l <limit>]\n";
    return 1;
  }

  omp_set_num_threads(1); // 设置使用指定线程数量
  std::string output_filename = output + filename + ".fa";
  processFile(inputFile, output_filename, k, limlt);
  return 0;
}

// 按照'N'分隔序列
std::vector<std::string> splitByN(const std::string &sequence) {
  std::vector<std::string> kmers;
  std::stringstream ss(sequence);
  std::string kmer;

  while (std::getline(ss, kmer, 'N')) {
    if (!kmer.empty()) {
      kmers.push_back(kmer);
    }
  }

  return kmers;
}

// 处理文件
void processFile(const std::string &inputFile, const std::string &outputFile,
                 int k, int limlt) {
  std::ifstream file(inputFile);
  std::ofstream outFile(outputFile);
  std::vector<std::string> result;
  std::string line;

  if (!file.is_open()) {
    std::cerr << "open filed: " << inputFile << std::endl;
  }
  std::cout << "first read ,no repeat and rev" << std::endl;

  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string sequence1, sequence2;
    int distance;

    // 读取每行的两个碱基串和最后的数字
    if (!(iss >> sequence1 >> sequence2 >> distance)) {
      std::cerr << "Incorrect format" << std::endl;
      continue;
    }
    // 对两个碱基串进行21-mer操作
    for (size_t i = 0; i + k <= sequence1.size(); ++i) {
      unsigned long long val = transfer_kmer_int(sequence1.substr(i, k));
      unsigned long long newval = reverse_int(val, k * 2);
      unsigned long long val1 = transfer_kmer_int(sequence2.substr(i, k));
      unsigned long long newval1 = reverse_int(val1, k * 2);
      if (!(kmer_count.find(val) != kmer_count.end())) {
        kmer_count[val] = 0;
        kmer_count[newval] = 0;
      }
      kmer_count[val] += 1;
      kmer_count[newval] += 1;

      if (!(kmer_count.find(val1) != kmer_count.end())) {
        kmer_count[val1] = 0;
        kmer_count[newval1] = 0;
      }
      kmer_count[val1] += 1;
      kmer_count[newval1] += 1;
    }
  }

  file.clear();                 // 清除eof和其他标志
  file.seekg(0, std::ios::beg); // 将文件指针移动到文件开头

  std::cout << "second read, Duplicate area and output" << std::endl;
  int number = 0;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string sequence1, sequence2;
    int distance;
    std::string modifiedSequence;
    std::string modifiedSequence1;
    int num = 0;
    int num1 = 0;

    // 读取每行的两个碱基串和最后的数字
    if (!(iss >> sequence1 >> sequence2 >> distance)) {
      std::cerr << "Incorrect format" << std::endl;
      continue;
    }
    for (size_t i = 0; i + k <= sequence1.size(); ++i) {

      unsigned long long val = transfer_kmer_int(sequence1.substr(i, k));
      unsigned long long newval = reverse_int(val, k * 2);
      unsigned long long val1 = transfer_kmer_int(sequence2.substr(i, k));
      unsigned long long newval1 = reverse_int(val1, k * 2);

      auto it1 = kmer_count.find(val);
      auto it2 = kmer_count.find(newval);
      if (it1 != kmer_count.end() && it1->second == 1 &&
          it2 != kmer_count.end() && it2->second == 1) {
        num++;
        if (!modifiedSequence.empty()) {
          modifiedSequence += 'N'; // 重新添加'N'分隔符
        }
        modifiedSequence += transfer_int_kmer(val, k);
      }
      auto it3 = kmer_count.find(val1);
      auto it4 = kmer_count.find(newval1);
      if (it3 != kmer_count.end() && it3->second == 1 &&
          it3 != kmer_count.end() && it3->second == 1) {

        num1++;
        if (!modifiedSequence1.empty()) {
          modifiedSequence1 += 'N'; // 重新添加'N'分隔符
        }
        modifiedSequence1 += transfer_int_kmer(val1, k);
      }
    }
    if (num >= limlt && num1 >= limlt) {
      outFile << (">" + std::to_string(number) + " " + "ref") << "\n"
              << modifiedSequence << "\n";
      outFile << (">" + std::to_string(number) + " " + "val") << "\n"
              << modifiedSequence1 << "\n";
      // number++;
    }
    number++;
  }

  file.close();
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

    reversed |= twoBits << (2 * (groupCount - 1 -
                                 i)); // 将第i组的两位放置到倒数第i组的位置
  }

  return reversed;
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