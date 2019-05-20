#ifndef _SPM4D_HASH_H
#define _SPM4D_HASH_H

#include <string>
#include <hash_map>
#include <list>
#include <math.h>

template <class T>
class spm4d_hash
{
public:
  
  struct coord {
    coord (size_t d1, size_t d2, size_t d3, size_t d4, size_t index)
      : d1(d1),d2(d2),d3(d3),d4(d4),index(index) {
    }

    bool operator== (const coord &right) const {
      return (d1==right.d1 && d2==right.d2 && d3==right.d3 && d4==right.d4);
    }

    size_t d1,d2,d3,d4,index;
  };
  
  struct hashfcn {
    size_t operator() (const coord &key) const {
      return key.index;
    }
  };

  typedef std::hash_map<coord, T, hashfcn> mat4d_type;
  typedef typename mat4d_type::const_iterator element_iter;  
  
  struct cell {
    size_t d1,d2,d3,d4;
    T value;
  };
  typedef std::list<cell> cell_list;

  spm4d_hash (size_t d1, size_t d2, size_t d3, size_t d4) 
    : d1_size(d1),d2_size(d2),d3_size(d3),d4_size(d4) {    
  }

  inline T& operator() (size_t d1, size_t d2, size_t d3, size_t d4) {
    coord key( d1,d2,d3,d4,((d1*d2_size+d2)*d3_size+d3)*d4_size+d4 );
    return mat[key];
  }

  inline T operator() (size_t d1, size_t d2, size_t d3, size_t d4) const {
    coord key( d1,d2,d3,d4,((d1*d2_size+d2)*d3_size+d3)*d4_size+d4 );
    return mat[key];
  }

  inline void increment_by(size_t d1, size_t d2, size_t d3, size_t d4, T value) {
    coord key( d1,d2,d3,d4,((d1*d2_size+d2)*d3_size+d3)*d4_size+d4 );
    T* v = &(mat[key]);
    (*v) = (*v)+value;
  }

  void convert_to_list(cell_list &cl) {
    element_iter i;
    cell new_cell;        

    for (i = mat.begin(); i != mat.end(); i++) {
      coord key = i->first;
      new_cell.value = i->second;
      
      new_cell.d1 = key.d1;
      new_cell.d2 = key.d2;
      new_cell.d3 = key.d3;
      new_cell.d4 = key.d4;

      cl.push_back(new_cell);
    }
  }

  ~spm4d_hash() {    
  }

public:
  mat4d_type mat;
  size_t d1_size,d2_size,d3_size,d4_size;
};

#endif