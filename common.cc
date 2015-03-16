#include "common.hpp"

// comparator for the <int, int> key pair
bool myeq::operator() (const std::pair<int, int> & x, const std::pair<int, int> & y) const{
  return x.first == y.first && x.second == y.second;
}

// hash function
size_t myhash::operator()(const std::pair<int, int> & p) const{
  size_t seed = h_int(p.first);
  return h_int(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// split the string by the separator c
std::vector<std::string> split_string(std::string s, std::string c){
  std::vector<std::string> ret;
  for(int i = 0, n = 0; i <= s.length(); i = n + 1){
    n = s.find_first_of(c, i);
    if(n == std::string::npos) n = s.length();
    std::string tmp = s.substr(i, n-i);
    ret.push_back(tmp);
  }
  return ret;
}

// get a uniform random number bewteen [0, 1)
double uniform_rand(){
  return (double)rand() * (1.0 / (RAND_MAX + 1.0));
}
