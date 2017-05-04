#ifndef DNA_H
#define DNA_H

#include <cstdint>

#include "utils.h"

inline uint convNuc(char nuc){
  switch (nuc){
    case 'a' : case 'A' : return 0 ;
    case 'c' : case 'C' : return 1 ;
    case 'g' : case 'G' : return 2 ;
    case 't' : case 'T' : return 3 ;
    default : return 0;
  }
  return 0;
}

inline uint compNuc(uint nuc){
  switch (nuc){
    case 0 : return 3 ;
    case 1 : return 2 ;
    case 2 : return 1 ;
    case 3 : return 0 ;
    default  : return 0;
  }
  return 0;
}

inline char intToNuc(uint c) {
  switch(c) {
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'T';
  default : return 'A';
  }
}

//int comparUint(const void *a1, const void* a2){
//  return (int) (* (uint64_t *)a1 - * (uint64_t *) a2);
//}

inline void intToDNA(uint64_t code, uint dna_length, char *dna) {
  uint64_t mask = 3;
  for (uint i=0; i < dna_length; i++) {
    dna[dna_length-i-1] = intToNuc(code & mask);
    code >>=2;
  }
}

inline uint64_t intRevcomp(uint64_t factor, uint32_t length) {
  uint64_t mask;
  if (length == 32)
    mask = ~0;
  else
    mask =  ((uint64_t) 1 << (2*length)) - 1;

  factor ^= mask;

  uint64_t mask_lsb;
  // Corresponds to the rightmost nucleotide
  mask_lsb = 3;
  uint64_t shift = 0;
  uint64_t result=0;
  for(uint32_t j(0);j<length;j++){
    result <<= 2;
    // get the leftmost nucleotide and put it at the end
    result |= (factor & mask_lsb) >> shift;
    mask_lsb <<= 2;
    shift += 2;
  }

  return result;
}

inline uint64_t DNAtoInt(const char *dna, uint32_t dna_length, bool stranded = false){
  uint64_t dna_int = 0;
  for (uint32_t i = 0; i< dna_length ; i++){
    dna_int <<= 2;
    dna_int |= convNuc(dna[i]);
  }
  // If the conversion is not "strand-specific" we calculate the reverse DNA it
  // and return the one that has the smallest value
  if(!stranded) {
    uint64_t rev_comp = intRevcomp(dna_int,dna_length);
    if(rev_comp < dna_int) {
      return rev_comp;
    }
  }
  return dna_int;
}

#endif // DNA_H_
