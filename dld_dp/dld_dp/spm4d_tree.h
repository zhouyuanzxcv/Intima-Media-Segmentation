#ifndef _SPM4D_TREE_H
#define _SPM4D_TREE_H

#include <string>
#include <map>
#include <list>


template <class T>
class spm4d_tree
{
public:
  typedef std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, T> > > > mat4d;
  typedef typename mat4d::iterator d1_iter;
  typedef typename std::map<size_t, std::map<size_t, std::map<size_t, T> > >::iterator d2_iter;
  typedef typename std::map<size_t, std::map<size_t, T> >::iterator d3_iter;
  typedef typename std::map<size_t, T>::iterator d4_iter;

  struct cell {
    size_t d1,d2,d3,d4;
    T value;
  };
  typedef std::list<cell> cell_list;

  spm4d_tree(size_t d1, size_t d2, size_t d3, size_t d4) 
    : d1(d1),d2(d2),d3(d3),d4(d4) {    
  }

  inline T& operator() (size_t d1, size_t d2, size_t d3, size_t d4) {
    return mat[d1][d2][d3][d4];
  }

  inline T operator() (size_t d1, size_t d2, size_t d3, size_t d4) const {
    return mat[d1][d2][d3][d4];
  }

  inline void increment_by(size_t d1, size_t d2, size_t d3, size_t d4, T value) {
    T* v = &(mat[d1][d2][d3][d4]);
    (*v) = (*v)+value;
  }

  void convert_to_list(cell_list &cl) {
    d1_iter i;
    d2_iter j;
    d3_iter k;
    d4_iter l;
    cell new_cell;

    for (i = this->mat.begin(); i != this->mat.end(); i++)
      for (j = i->second.begin(); j != i->second.end(); j++)
        for (k = j->second.begin(); k != j->second.end(); k++)
          for (l = k->second.begin(); l != k->second.end(); l++){
            new_cell.d1 = i->first;
            new_cell.d2 = j->first;
            new_cell.d3 = k->first;
            new_cell.d4 = l->first;   
            new_cell.value = l->second;
            cl.push_back(new_cell);
          }    
    
  }

  void max_pos(size_t &d1, size_t &d2, size_t &d3, size_t &d4) {
    d1_iter i;
    d2_iter j;
    d3_iter k;
    d4_iter l;
    
    T m = 0;
    for (i = this->mat.begin(); i != this->mat.end(); i++)
      for (j = i->second.begin(); j != i->second.end(); j++)
        for (k = j->second.begin(); k != j->second.end(); k++)
          for (l = k->second.begin(); l != k->second.end(); l++){
            if (l->second > m) {
              m = l->second;
              d1 = i->first;
              d2 = j->first;
              d3 = k->first;
              d4 = l->first;              
            }
          }    
  }

  ~spm4d_tree() {    
  }

private:
  mat4d mat;
  size_t d1,d2,d3,d4;    
};

#endif