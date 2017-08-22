#include <cmath>
#include <iostream>
#include <assert.h>
#include <iterator>
#include <math_utils.h>

typedef std::vector<double> State;
typedef std::vector<double> Control;
typedef std::vector< std::vector<double> > Matrix;
typedef std::vector< State > Trajectory;
typedef std::vector< Control > Input;
typedef std::vector< std::vector < double > > Jacobian;
typedef std::vector< double > DecisionVar;
typedef DecisionVar::iterator Iter;

void printVector(  std::vector<double>& xdot){
  std::cout << "(";
  for(int i=0;i<xdot.size()-1;i++){
    std::cout << xdot[i] << "," << std::endl;
  }
  std::cout << xdot[xdot.size()-1] << ")" << std::endl;
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

std::vector<double> reshape(std::vector< std::vector<double> > matrix){
  assert(matrix.size()>0);
  
  std::vector<double> mat2vec;
  printMatrix(matrix);
  for(int i=0;i<matrix.size();i++){
    for(int j=0;j<matrix[0].size();j++){
      mat2vec.push_back(matrix[i][j]);
    }
  }
  return mat2vec;
}

class DynamicalSystem {
protected: 
  unsigned int stateDim;
  unsigned int controlDim;
  unsigned int numSteps;
  double dt;
  
public:
  //ctor
  DynamicalSystem(unsigned int stateDim_, 
                  unsigned int controlDim_,
                  unsigned int numSteps_,
                  double dt_):
                  stateDim(stateDim_), 
                  controlDim(controlDim_),
                  numSteps(numSteps_),
                  dt(dt_)
                  {}
                  
                  //evaluate the vector field at the state x_ with the input u_
                  virtual std::vector<double> vectorField(std::vector<double> x_, std::vector<double> u_)=0;
                                    
                  std::vector<double> vectorField(DecisionVar y, unsigned int index) {
                    State x(stateDim);
                    Control u(controlDim);
                    
                    for(int i=0;i<stateDim;i++){
                      x[i] = y[index];
                      index++;
                    }
                    for(int i=0;i<controlDim;i++){
                      u[i] = y[index];
                      index++;
                    }
                    
                    return vectorField(x,u);
                  }                   
                  
                  //evaluate the jacobian of the vector field at x_ with respect to x_ and u_
                  virtual Matrix systemJacobian(  std::vector<double> x_,   std::vector<double> u_)=0;
                  
                  //vector structure [x[0] u[0] x[1] u[1] ... dt]
                  std::vector<double> constraintResidual(DecisionVar& y){
                    int n = stateDim;
                    int m = controlDim;
                    int N = numSteps;
                    assert(y.size()==(N+1)*(n+m)+1);
                    double dt = y[(N+1)*(n+m)];
                    
                    std::vector<double> residual;
                    //time step j from 0 to N-2 
                    for(int j=0;j<N;j++){
                      //state i from 0 to n-1
                      State xdot = vectorField(y,j*(n+m));
                      for(int i=0;i<n;i++){
                        double error = y[(j+1)*(n+m)+i]-y[j*(n+m)+i]-y[(N+1)*(n+m)]*xdot[i];
                        residual.push_back(error);
                      }
                    }
                    return residual;
                  }
                  
                  Matrix dynamicInfeasibility(  Trajectory x,   Input u,   double dt){
                    assert(x.size()==u.size());
                    
                    //x[i+1] = x[i] + dt*f(x,u) => res[i] = x[i+1]-x[i]-dt*f(x,u)
                    Matrix residual(x.size()-1);
                    for(int i=0;i<x.size()-1;i++){
                      residual[i] = x[i+1] - x[i] - dt*vectorField(x[i],u[i]);
                    }
                    return residual;
                  }
                  
                  double timeStep(){return dt;}
                  unsigned int steps(){return numSteps;}
                  unsigned int state_dimension(){return stateDim;}
                  unsigned int control_dimension(){return controlDim;}
};

class Pendulum : public DynamicalSystem {

public:
  Pendulum(unsigned int numSteps, double dt) : DynamicalSystem(2, 1,numSteps,dt) {}
  
  std::vector<double> vectorField(std::vector<double> x, std::vector<double> u) final {
    std::vector<double> xdot(stateDim);
    xdot[0] = x[1];
    xdot[1] = -sin(x[0])+u[0];
    
    return xdot;
  }
  
  Matrix systemJacobian(State x, Control u) final {
    assert(x.size()==stateDim);
    assert(u.size()==controlDim);
    //create 2x3 matrix for jacobian
    Matrix Jac(stateDim);
    for(int i=0;i<stateDim;i++){Jac[i].resize(controlDim+stateDim);}
    
    Jac[0][0] = 0; 
    Jac[0][1] = 1; 
    Jac[0][2] = 0;
    Jac[1][0] = cos(x[0]); 
    Jac[1][1] = 0;
    Jac[1][2] = 1;
    
    return Jac;
  }
};