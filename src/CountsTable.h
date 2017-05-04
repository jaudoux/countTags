#ifndef COUNTS_TABLE_H
#define COUNTS_TABLE_H

#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "utils.h"
#include "dna.h"

class CountsTable {

private:
  std::unordered_map<uint64_t,uint32_t*> tags_counts;
  std::vector<std::string> sample_names;
  uint nb_samples;
  uint kmer_length;
  bool stranded;

public:
  CountsTable(uint nb_samples, uint kmer_length, bool stranded = false);

  ~CountsTable();

  void setSampleName(uint sample_id, const char *sample_name);

  uint32_t getCount(const char *kmer, uint sample_id);
  uint32_t getCount(uint64_t kmer, uint sample_id);

  /*
  * Set the counts value for one kmer of one sample
  */
  void setCount(const char *kmer, uint sample_id, uint value);
  void setCount(uint64_t kmer, uint sample_id, uint value);

  /*
  * Increment the count for one kmer of one sample
  */
  void incrementCount(const char *kmer, uint sample_id, uint value);
  void incrementCount(uint64_t kmer, uint sample_id, uint value);

  /*
  * Filters counts based on their recurrency accross samples
  */
  void recurrencyFilter(uint recurrency_threshold = 2);

  /*
  * Print counts on STDOUT
  */
  void printCounts(char sep = '\t');

};

#endif // COUNTS_TABLE_H
