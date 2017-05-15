#include <iostream>
#include <sstream>
#include <iterator>
#include <string.h>
#include <cstdint>
#include <vector>
#include <unordered_map>
//#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>

#include "optionparser.h"
#include "dna.h"
//#include <zlib.h>

#define MILLION 1000000

//return the minumum value of the k-mer at pos p between strand rev and stran fwd
//TODO add a function that get a DNA string and a k value, and return a array of vector values
inline uint64_t valns(uint32_t p, char *dna,uint32_t k,int64_t *last,uint64_t *valfwd,uint64_t *valrev, bool stranded = false){
  int e=p-*last;
  if(e!=1){
    *last=p;
    *valfwd=DNAtoInt(&dna[p], k, true);
    *valrev=intRevcomp(*valfwd, k);
  }else{
    // Compute the new value from the previous one.
    uint64_t m=1;
    *valfwd%=m<<(2*k-2);
    *valfwd<<=2;
    int new_nuc = convNuc(dna[p+k-1]);
    *valfwd += new_nuc;
    *last=p;
    *valrev/=1<<(2);
    *valrev+=(uint64_t)compNuc(new_nuc)<<(2*k-2);
  }
  if(stranded || *valfwd < *valrev) {
    return *valfwd;
  } else {
    return *valrev;
  }
}

std::string join( const std::vector<std::string>& elements, const char* const separator)
{
  switch (elements.size())
  {
    case 0:
      return "";
    case 1:
      return elements[0];
    default:
      std::ostringstream os;
      std::copy(elements.begin(), elements.end()-1, std::ostream_iterator<std::string>(os, separator));
      os << *elements.rbegin();
      return os.str();
  }
}

struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "ERROR: %s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }

  static option::ArgStatus Numeric(const option::Option& option, bool msg)
  {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }
};

enum  optionIndex {UNKNOWN,HELP,PROBE_LENGTH,STRANDED,MAX_READS,NB_THREADS,NORMALIZE,TAG_NAMES};
const option::Descriptor usage[] =
{
  {UNKNOWN, 0,"" , ""    ,
    option::Arg::None, "USAGE: countTags [options] tags.fa seq.fastq[.gz]\n\nOptions:" },
  {HELP,    0,"" , "help",
    option::Arg::None, "  --help  \tPrint usage and exit." },
  {PROBE_LENGTH, 0, "k","",
    Arg::Numeric,      "  -k INT      \ttags length" },
  {MAX_READS, 0, "m","",
    Arg::Numeric,      "  -m INT      \tmax number of reads" },
  //{NB_THREADS, 0, "t","",
  //  Arg::Numeric,      "  -t=INT      \tnumber of threads" },
  {STRANDED,    0,"" , "stranded",
    option::Arg::None, "  --stranded  \tstrand-specific protocol" },
  {NORMALIZE,    0,"" , "normalize",
    option::Arg::None, "  --normalize  \tnormalize counts" },
  {TAG_NAMES,    0,"" , "tag-names",
    option::Arg::None, "  --tag-names  \tprint tag names in the output" },
  {0,0,0,0,0,0}
};

class Tag {
public:
  uint64_t tag;
};

int main (int argc, char *argv[]) {

  // Config vars
  const char * tags_file;
  const char * seq_file;
  uint32_t tag_length = 22;
  bool stranded = false;
  bool normalize = false;
  bool print_tag_names = false;
  uint32_t nb_tags = 0;
  uint32_t max_reads = UINT32_MAX;
  int nb_threads = 1;
  uint32_t nb_samples;
  uint32_t i;
  uint32_t read_id;

  char * seq;
  std::string tag_name;
  uint32_t seq_length;
  uint64_t tag;
  uint64_t valrev,valfwd;
  uint64_t nb_factors;
  int64_t last;
  std::unordered_map<uint64_t,double*>  tags_counts;
  std::unordered_map<uint64_t,std::vector<std::string>> tags_names;
  std::unordered_map<uint64_t,double*>::iterator it_counts;
  std::vector<uint64_t> nb_factors_by_sample;

  // File vars
  std::string gzip_pipe = "gunzip -fc ";
  std::string cat_pipe = "cat ";
  std::string tmp;
  FILE * file;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  uint32_t line_id = 0;

  /**********************************
   *
   *        Parsing options
   *
   *********************************/
  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  option::Stats  stats(usage, argc, argv);
  option::Option options[stats.options_max], buffer[stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);

  if (parse.error())
    return 1;

  if (options[HELP] || argc == 0) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  // TODO we should check the length to choose the appropriate
  // integer size to use
  if (options[PROBE_LENGTH]) {
    tag_length = atoi(options[PROBE_LENGTH].arg);
  }

  if (options[MAX_READS]) {
    max_reads = atoi(options[MAX_READS].arg);
  }

  //if (options[NB_THREADS]) {
  //  nb_threads = atoi(options[NB_THREADS].arg);
  //}

  if (options[STRANDED]) {
    stranded = true;
  }

  if (options[NORMALIZE]) {
    normalize = true;
  }

  if (options[TAG_NAMES]) {
    print_tag_names = true;
  }

  if(parse.nonOptionsCount() < 2) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  tags_file = parse.nonOption(0);
  nb_samples = parse.nonOptionsCount() - 1;

  /**********************************
   *
   *   Create tags counting table
   *
   *********************************/
  file = popen(tmp.append(cat_pipe).append(tags_file).c_str(),"r");
  line_id = 0;

  std::cerr << "Counting k-mers" << std::endl;
  // Create hash table of k-mer counts
  line_id = 0;
  nb_tags = 0;
  while ((read = getline(&line, &len, file)) != -1) {
    if(line_id % 2 == 0) {
      tag_name = line;
      tag_name.pop_back(); // remove new_line char
      tag_name.erase(0,1); // remove the fasta ">" prefix
    }
    if(line_id % 2 == 1) {
      tag = DNAtoInt(line,tag_length,stranded);
      tags_counts[tag] = new double[nb_samples]();
      tags_names[tag].push_back(tag_name);
      nb_tags++;
    }
    line_id++;
  }
  fclose(file);
  if (line)
    free(line);

  std::cerr << "Finished indexing tags" << std::endl;


  /**********************************
   *
   *            First pass
   *
   *********************************/
//#pragma omp parallel num_threads(nb_threads)
  for (int s = 0; s < nb_samples; ++s) {
    std::cerr << "Counting tags for file: " << parse.nonOption(s+1) << "\n";

    line = NULL;
    len = 0;
    line_id = 0;
    tmp = "";
    file = popen(tmp.append(gzip_pipe).append(parse.nonOption(s+1)).c_str(),"r");
    nb_factors = 0;

    while ((read = getline(&line, &len, file)) != -1) {
      // If this line is a sequence
      if(line_id % 4 == 1) {
        read_id = ((int)((double)line_id*0.25) + 1);
        if(read_id >= max_reads) {
          break;
        }
        // Print a user-friendly output on STDERR every each XXXX reads processed
        if(read_id % MILLION == 0) {
          std::cerr << (int)((double)line_id*0.25) + 1 << " reads parsed" << std::endl;
        }
        // set seq to line
        seq = line;
        seq_length = strlen(seq) - 1; // Minus 1 because we have a new line

        // Skip the sequence if the read length is < to the tag_length
        if(seq_length < tag_length)
          continue;

        nb_tags = seq_length - tag_length + 1;
        nb_factors += nb_tags;

        //uint64_t valrev,valfwd;
        last = -3;

        for(i = 0; i < nb_tags; i++) {
          it_counts = tags_counts.find(valns(i,seq,tag_length,&last,&valfwd,&valrev,stranded));
          if(it_counts != tags_counts.end()) {
            it_counts->second[s]++;
          }
        }
      }
      line_id++;
    }

    nb_factors_by_sample.push_back(nb_factors);

    if(normalize && nb_factors > 0) {
      std::cerr << "Normalize counts" << std::endl;
      for (it_counts=tags_counts.begin(); it_counts!=tags_counts.end(); ++it_counts) {
        // TODO We should take into accout the error rate...
        if(it_counts->second[s] > 0)
          it_counts->second[s] = it_counts->second[s] * MILLION / nb_factors;
      }
    }

    // Close file and clear line buffer
    fclose(file);
    if (line)
      free(line);
  }

  /****
   * PRINT THE RESULTS
   */
  // First print headers
  std::cout << "tag";
  if (print_tag_names)
    std::cout << "\ttag_names";
  for (int s = 0; s < nb_samples; ++s) {
    std::cout << "\t" << parse.nonOption(s+1);
  }
  std::cout << "\n";
  char *tag_seq = new char[tag_length];
  for (it_counts=tags_counts.begin(); it_counts!=tags_counts.end(); ++it_counts) {
    intToDNA(it_counts->first,tag_length,tag_seq);
    std::cout << tag_seq;
    if(print_tag_names) {
      std::cout << "\t" << join(tags_names[it_counts->first],",");
    }
    for (int s = 0; s < nb_samples; ++s) {
      std::cout << "\t" << it_counts->second[s];
    }
    std::cout << std::endl;
  }

  std::cout << "total_factors";
  if (print_tag_names)
    std::cout << "\t*";
  for (int s = 0; s < nb_samples; ++s) {
    std::cout << "\t" << nb_factors_by_sample[s];
  }
  std::cout << std::endl;

}
