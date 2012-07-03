#ifndef __INDEX_VECTOR_HPP__
#define __INDEX_VECTOR_HPP__

struct c_unique {
  int current;
  c_unique() {current=0;}
  int operator()() {return ++current;}
};

inline std::vector<int> returnIndexVector(size_t n)
{
  c_unique UniqueNumber;
  std::vector<int> arr(n);
  generate (arr.begin(), arr.end(), UniqueNumber);
  return arr;
}

inline void modifyIndexVector(std::vector<int>& arr)
{
  c_unique UniqueNumber;
  generate (arr.begin(), arr.end(), UniqueNumber);
}

#endif
