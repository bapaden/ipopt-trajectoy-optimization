#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <vector>
#include <iostream>

//Vector of vectors of doubles (i.e. 2D double array). Adds a constructor for resizing matrix
class Matrix{
  std::vector< std::vector<double> > data;
public:
//   Matrix(){}
  Matrix(unsigned int rows,unsigned int columns){
    data.resize(rows);
    for(auto& row : data)
      row.resize(columns);
    
    assert(rows>0 and columns >0);
    assert(data.size()==rows);
    assert(data[0].size()==columns);
  }
 std::vector<double>& operator[](unsigned int index){return data[index];}
 unsigned int size(){return data.size();}
};

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

void printVector(  std::vector<double>& xdot){
  std::cout << "~" <<std::endl;
  for(int i=0;i<xdot.size();i++){
    std::cout << xdot[i] << std::endl;
  }
  std::cout << "~" << std::endl;
}

void printMatrix(   Matrix& matrix ){
  if(matrix.size()>0){
    for(int i=0;i<matrix.size();i++){
      std::cout << "|";
      for(int j=0;j<matrix[i].size();j++){
        std::cout << matrix[i][j] << " ";
      }
      std::cout << "|" << std::endl;
    }
  }
}

#endif