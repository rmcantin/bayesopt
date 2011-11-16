#ifndef __INDEX_VECTOR_HPP__
#define __INDEX_VECTOR_HPP__

struct c_unique {
  int current;
  c_unique() {current=0;}
  int operator()() {return ++current;}
} UniqueNumber;

inline std::vector<int> returnIndexVector(size_t n)
{
  std::vector<int> arr(n);
  generate (arr.begin(), arr.end(), UniqueNumber);
  return arr;
}

inline void modifyIndexVector(std::vector<int>& arr)
{
  generate (arr.begin(), arr.end(), UniqueNumber);
}

#endif
