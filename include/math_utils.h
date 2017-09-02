#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <vector>
#include <iostream>
#include <cstdio>

  typedef std::vector<double> Vector ;
  
  class tester{
    int count;
  public:
    tester(){count=0;}
    void msg(){std::cout << "Test msg: " << count << std::endl;count++;}
  };
  
  class Matrix;
  void printMatrix(Matrix mat);
  
  //Vector of vectors of doubles (i.e. 2D double array). Adds a constructor for resizing matrix
  class Matrix{
    std::vector< Vector > data;
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
    Vector& operator[](unsigned int index){return data[index];}
    unsigned int size(){return data.size();}
    bool empty(){return data.empty();}
    void append_bottom(Matrix mat){
      data.insert(data.end(),mat.data.begin(),mat.data.end());
    }
    void append_right(Matrix mat){
//       std::cout << "before " << std::endl;
//       printMatrix(*this);      
      assert(mat.size()==data.size());
      for(int i=0;i<data.size();i++){
        data[i].insert(data[i].end(),mat[i].begin(),mat[i].end());
      }
//       std::cout << "after " << std::endl;
//       printMatrix(*this);
      
    }
  };
  
  Vector operator+(Vector lhs, Vector rhs){
    assert(lhs.size()==rhs.size());
    Vector sum;
    for(int i=0;i<rhs.size();i++)
      sum.push_back(lhs[i]+rhs[i]);
    return sum;
  }
  
  Vector operator-(Vector lhs,  Vector rhs){
    assert(lhs.size()==rhs.size());
    Vector difference;
    for(int i=0;i<rhs.size();i++)
      difference.push_back(lhs[i]-rhs[i]);
    return difference;
  }
  
  
  
  Vector operator*(double c, Vector rhs){
    Vector scaled;
    for(int i=0;i<rhs.size();i++)
      scaled.push_back(c*rhs[i]);
    return scaled;
    
  }
  
  void printVector(  Vector& xdot){
    std::cout << "~" <<std::endl;
    for(int i=0;i<xdot.size();i++){
      std::cout << xdot[i] << std::endl;
    }
    std::cout << "~" << std::endl;
  }
  
  void printMatrix( Matrix matrix ){
    std::cout << "size: (" << matrix.size() << "," << matrix[0].size() << std::endl;
    if(matrix.size()>0){
      for(int i=0;i<matrix.size();i++){
        std::cout << "|";
        for(int j=0;j<matrix[i].size();j++){
          std::cout << double(matrix[i][j]) << " ";
        }
        std::cout << "|" << std::endl;
      }
    }
  }
  
  
#endif