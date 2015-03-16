#include <iostream>
#include <vector>

typedef std::pair<int, int> key;

struct myeq : std::binary_function<std::pair<int, int>, std::pair<int, int>, bool>{
  bool operator() (const std::pair<int, int> & x, const std::pair<int, int> & y) const;
};

struct myhash : std::unary_function<std::pair<int, int>, size_t>{
private:
  const std::hash<int> h_int;
public:
  myhash() : h_int() {}
  size_t operator()(const std::pair<int, int> & p) const;
};

std::vector<std::string> split_string(std::string s, std::string c);

double uniform_rand();