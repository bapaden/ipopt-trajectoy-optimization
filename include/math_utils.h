#include <vector>

std::vector<double> operator+(std::vector<double> lhs, std::vector<double> rhs){
  assert(lhs.size()==rhs.size());
  std::vector<double> sum;
  for(int i=0;i<rhs.size();i++)
    sum.push_back(lhs[i]+rhs[i]);
  return sum;
}

std::vector<double> operator-(std::vector<double> lhs,  std::vector<double> rhs){
  assert(lhs.size()==rhs.size());
  std::vector<double> difference;
  for(int i=0;i<rhs.size();i++)
    difference.push_back(lhs[i]-rhs[i]);
  return difference;
}

std::vector<double> operator*(double c, std::vector<double> rhs){
  std::vector<double> scaled;
  for(int i=0;i<rhs.size();i++)
    scaled.push_back(c*rhs[i]);
  return scaled;
  
}