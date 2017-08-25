#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

#include <cmath>
#include <iostream>
#include <assert.h>
#include <iterator>
#include <math_utils.h>

typedef std::vector<double> State;
typedef std::vector<double> Control;
typedef std::vector< State > Trajectory;
typedef std::vector< Control > Input;
typedef std::vector< std::vector < double > > Jacobian;
typedef std::vector< double > DecisionVar;
typedef DecisionVar::iterator Iter;

class DynamicalSystem {
protected: 
  unsigned int stateDim,controlDim,numSteps,n,m,N;
  double dt;
  std::vector<double> initCond;
  std::vector<double> termCond;
  
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
                  {
                    n = stateDim;
                    m = controlDim;
                    N = numSteps;
                  }
                  
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
                  virtual Matrix systemJacobian(DecisionVar y, unsigned int index)=0;
                  
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
                  
                  Matrix constraintJacobian(DecisionVar y){
                    
                    int n = stateDim;
                    int m = controlDim;
                    int N = numSteps;
                    double dt = y[(N+1)*(n+m)];
                    double count = 0;
                    Matrix dgdy(n*N,(N+1)*(n+m));
                    for(int p=0;p<N;p++){
                      State xdot = vectorField(y,p*(n+m));
                      for(int q=0;q<n;q++){
                        dgdy[p*n+q][(N+1)*(n+m)-1] = xdot[q];
                        for(int j=0;j<N+1;j++){
                          Matrix dfdx = systemJacobian(y,j*(m+n));
                          for(int i=0;i<n+m;i++){
                            if(p==j){
                              if(i<n){dgdy[p*n+q][j*(n+m)+i]= -1.0 - dt*dfdx[q][i];}
                              else{dgdy[p*n+q][j*(n+m)+i]=dt*dfdx[q][i];}
                            }
                            if(j-1==p and i == q){dgdy[p*n+q][j*(n+m)+i]=1;}
                          }
                        }
                      }    
                    }  
                        
                    return dgdy;
                  }
                  
                  //vector of constraints of the form h(x(t),u(t))<=0
                  virtual std::vector<double> signalConstraint(std::vector<double> x, std::vector<double> u)=0;
                  //Apply the user specified constraint to each timestep
                  std::vector<double> constraints(DecisionVar y){
                    std::vector<double> h;
                    std::vector<double>::const_iterator x_first,x_last,u_first,u_last;
                    for(int j=0;j<N+1;j++){

                      x_first = y.begin()+j*(n+m);
                      x_last = x_first+n-1;
                      u_first = x_last+1;
                      u_last = u_first+m-1;
                      
                      std::vector<double> x(x_first,x_last+1);
                      std::vector<double> u(u_first,u_last+1);
                      std::vector<double> h_t( signalConstraint(x,u) );
                      h.insert(h.end(),h_t.begin(),h_t.end());
                    }
                    return h;
                  }
                  
                  std::vector<double> constraintJacobian(DecisionVar y){
                    
                  }
                  
                  double timeStep(){return dt;}
                  unsigned int steps(){return numSteps;}
                  unsigned int state_dimension(){return stateDim;}
                  unsigned int control_dimension(){return controlDim;}
};

#endif