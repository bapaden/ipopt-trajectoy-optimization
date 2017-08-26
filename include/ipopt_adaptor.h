#ifndef IPOPT_ADAPTOR_H
#define IPOPT_ADAPTOR_H

#include <cmath>
#include <iostream>
#include <assert.h>
#include <iterator>
#include <math_utils.h>


class TPBVP{
protected: 
  unsigned int stateDim,controlDim,numSteps,n,m,N;
  double dt;
  Vector initCond;
  Vector termCond;
  
public:
  //ctor
  TPBVP(unsigned int stateDim_, 
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
                  
                  virtual Matrix bv()=0;
                  
                  //For initializing optimization recursion
                  virtual Matrix initialTraj()=0;
                  Vector initialGuess(){
                    Matrix init = initialTraj();
                    assert(not init.empty());
                    
                    Vector init_vec;
                    for(int i=0;i<init.size();i++){
                      for(int j=0;j<init[0].size();j++){
                        init_vec.push_back(init[i][j]);
                      }
                    }
                    return init_vec;
                  }
                  
                  //Running cost for cost function
                  virtual double g(Vector x, Vector u)=0;
                  double g(Vector y, unsigned int index){
                    Vector x(stateDim);
                    Vector u(controlDim);
                    
                    for(int i=0;i<stateDim;i++){
                      x[i] = y[index];
                      index++;
                    }
                    for(int i=0;i<controlDim;i++){
                      u[i] = y[index];
                      index++;
                    }
                    return g(x,u);
                  }
                  
                  double J(Vector y){
                    double cost(0);
                    for(int j=0;j<N+1;j++){
                      cost+=g(y,j*(n+m));
                    }
                    return cost;
                  }
                  
                  
                  //evaluate the vector field at the state x_ with the input u_
                  virtual Vector f(Vector x_, Vector u_)=0;
                                    
                  Vector f(Vector y, unsigned int index) {
                    Vector x(stateDim);
                    Vector u(controlDim);
                    
                    for(int i=0;i<stateDim;i++){
                      x[i] = y[index];
                      index++;
                    }
                    for(int i=0;i<controlDim;i++){
                      u[i] = y[index];
                      index++;
                    }
                    return f(x,u);
                  } 
                  
                  //vector of constraints of the form h(x(t),u(t))<=0
                  virtual Vector h(Vector x, Vector u)=0;
                  
                  Vector h(Vector y, unsigned int index){
                    Vector x(n);
                    Vector u(m);
                    
                    for(int i=0;i<n;i++){
                      x[i] = y[index];
                      index++;
                    }
                    for(int i=0;i<n;i++){
                      u[i] = y[index];
                      index++;
                    }
                    return h(x,u);
                  }
                  
                  //vector structure [x[0] u[0] x[1] u[1] ... dt]
                  Vector phi(Vector& y){
                    int n = stateDim;
                    int m = controlDim;
                    int N = numSteps;
                    assert(y.size()==(N+1)*(n+m)+1);
                    double dt = y[(N+1)*(n+m)];
                    
                    Vector residual;
                    //time step j from 0 to N-2 
                    for(int j=0;j<N;j++){
                      //state i from 0 to n-1
                      Vector xdot = f(y,j*(n+m));
                      for(int i=0;i<n;i++){
                        double error = y[(j+1)*(n+m)+i]-y[j*(n+m)+i]-y[(N+1)*(n+m)]*xdot[i];
                        residual.push_back(error);
                      }
                    }
                    return residual;
                  }
                  
                  ///~~~Gradients~~~///
                  
                  //Gradient of the running cost
                  virtual Vector Dg(Vector x, Vector y)=0;
                  Vector Dg(Vector y, unsigned int index){
                    Vector x(n);
                    Vector u(m);
                    
                    for(int i=0;i<n;i++){
                      x[i] = y[index];
                      index++;
                    }
                    for(int i=0;i<n;i++){
                      u[i] = y[index];
                      index++;
                    }
                    return Dg(x,u);
                  }
                    
                  
                  //First variation of cost
                  Vector DJ(Vector y){
                    Vector dJdx(y.size());
                    for(int j=0;j<N+1;j++){
                      Vector gradRunningCost = Dg(y,j*(n+m));
                      std::cout << "local grad " << std::endl;
                      printVector(gradRunningCost);
                      for(int i=0;i<n+m;i++){
                        dJdx[j*(n+m)+i] = gradRunningCost[i];
                      }
                    }
                    return dJdx;
                  }
                  
                  //From the problem adaptor
                  virtual Matrix Df(Vector x, Vector u)=0;
                  //Overload for conveniece
                  Matrix Df(Vector y, unsigned int index){
                    Vector x(n);
                    Vector u(m);
                    
                    for(int i=0;i<n;i++){
                      x[i] = y[index];
                      index++;
                    }
                    for(int i=0;i<n;i++){
                      u[i] = y[index];
                      index++;
                    }
                    return Df(x,u);
                  }
                  
                  //Numerical approximation to first variation diff constraint with respect to (x,u)
                  Matrix Dphi(Vector y){
                    
                    int n = stateDim;
                    int m = controlDim;
                    int N = numSteps;
                    double dt = y[(N+1)*(n+m)];
                    double count = 0;
                    Matrix dgdy(n*N,(N+1)*(n+m));
                    for(int p=0;p<N;p++){
                      Vector xdot = f(y,p*(n+m));
                      for(int q=0;q<n;q++){
                        dgdy[p*n+q][(N+1)*(n+m)-1] = xdot[q];
                        for(int j=0;j<N+1;j++){
                          Matrix dfdx = Df(y,j*(m+n));
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
                  
                  //User provides gradient of state constraints
                  virtual Matrix Dh(Vector x, Vector u)=0;
                  //Overload for convenience
                  Matrix Dh(Vector y, unsigned int index){
                    Vector x(n);
                    Vector u(m);
                    
                    for(int i=0;i<n;i++){
                      x[i] = y[index];
                      index++;
                    }
                    for(int i=0;i<n;i++){
                      u[i] = y[index];
                      index++;
                    }
                    
                    return Dh(x,u);
                  }
                  
                  //Apply the user specified constraint to each timestep
                  Vector H(Vector y){
                    Vector H_;
                    Vector::const_iterator x_first,x_last,u_first,u_last;
                    for(int j=0;j<N+1;j++){

                      x_first = y.begin()+j*(n+m);
                      x_last = x_first+n-1;
                      u_first = x_last+1;
                      u_last = u_first+m-1;
                      
                      Vector x(x_first,x_last+1);
                      Vector u(u_first,u_last+1);
                      Vector h_t( h(x,u) );
                      H_.insert(H_.end(),h_t.begin(),h_t.end());
                    }
                    return H_;
                  }
                  
                  //Essentially the first variation of the state-control constraints
                  Matrix DH(Vector y){
                    std::size_t k = (h(y,0)).size();
                    std::cout << "k" << k << std::endl;
                    Matrix dHdx(k*(N+1),(n+m)*(N+1));
                    for(int p=0;p<N+1;p++){
                      for(int q=0;q<k;q++){
                        Matrix dhdx = Dh(y,p*(n+m));
                        for(int j=0;j<N+1;j++){
                          for(int i=0;i<n+m;i++){
                            if(p==j){dHdx[p*k+q][j*(n+m)+i]=dhdx[q][i];}
                          }
                        }
                      }
                    }
                    return dHdx;
                  }
                  
                  double timeStep(){return dt;}
                  unsigned int steps(){return numSteps;}
                  unsigned int state_dimension(){return stateDim;}
                  unsigned int control_dimension(){return controlDim;}
};

#endif