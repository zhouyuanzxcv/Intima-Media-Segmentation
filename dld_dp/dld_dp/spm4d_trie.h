#ifndef _SPM4D_TRIE_H
#define _SPM4D_TRIE_H

#include <string>
#include <list>
#include <math.h>
#include "carray.h"
#include "trie_map.h"

template <class T>
class spm4d_trie
{
public:
  typedef carray<unsigned char, 4> key_type;
  typedef trie_map<key_type, T, 256> mat_type;
  typedef typename mat_type::iterator element_iter;
  
  struct cell {
    size_t d1,d2,d3,d4;
    T value;
  };
  typedef std::list<cell> cell_list;
  

  spm4d_trie(size_t d1, size_t d2, size_t d3, size_t d4) 
    : d1_size(d1),d2_size(d2),d3_size(d3),d4_size(d4) {          
  }

  inline T& operator() (size_t d1, size_t d2, size_t d3, size_t d4) {
    key[0] = d1; key[1] = d2; key[2] = d3; key[3] = d4;
    return mat[key];
  }

  inline T operator() (size_t d1, size_t d2, size_t d3, size_t d4) const {
    key[0] = d1; key[1] = d2; key[2] = d3; key[3] = d4;
    return mat[key];
  }

  inline void increment_by(size_t d1, size_t d2, size_t d3, size_t d4, T value) {
    key[0] = d1; key[1] = d2; key[2] = d3; key[3] = d4;    
    T* v = &(mat[key]);
    (*v) = (*v)+value;
  }

  void convert_to_list(cell_list &cl) {
    element_iter i;    
    element_iter b = mat.begin();
    element_iter e = mat.end();
    cell new_cell;

    for (i = b; i != e; ++i) {
      new_cell.d1 = i->first[0];
      new_cell.d2 = i->first[1];
      new_cell.d3 = i->first[2];
      new_cell.d4 = i->first[3];
      new_cell.value = i->second;
      cl.push_back(new_cell);
    }
  }

  ~spm4d_trie() {    
  }

public:
  mat_type mat;  
  key_type key;
  size_t d1_size,d2_size,d3_size,d4_size;  
};

#endif