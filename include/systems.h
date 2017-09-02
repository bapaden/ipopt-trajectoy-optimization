#ifndef SYSTEMS_H
#define SYSTEMS_H

#include <ipopt_adaptor.h>

/*
 *Simple example:
 *Dynamics: dx/dt = -x + u 
 *Boundary Conditions: x(0)=1 x(tf)0
 *Cost: Integral of g(x(t),u(t))=1 + x^2 + u^2 
 *Constraints: h_1(x(t),u(t) = 0.5-u(t) h_2(t) = x(t)-1
 */

class FirstOrderSystem : public TPBVP {
public:
  FirstOrderSystem(unsigned int numSteps, double dt) : TPBVP(1,1,numSteps,dt){}

  //~~~Cost, Dynamics, Constraints, Boundary conditions~~~//
  
  //Define the running cost g(x(t),u(t)) whose integral is the cost
  double g(Vector x, Vector u, double dt) final {
    return 1;
  }
  
  //System dynamics dx/dt = f(x(t),u(t)) 
  Vector  f(Vector  x, Vector  u) final {
    Vector  xdot(n);
    xdot[0] = u[0];
    return xdot;
  }
  
  //State and control constraints h(x(t),u(t))<=0 
  Vector  h(Vector  x, Vector  u) final {
    assert(x.size()==n and u.size()==m);
    Vector  h;
    h.push_back(-2.0-x[0]);
    h.push_back(-u[0]-2.0);
    
    return h;
  }
  
  //Boundary values
  Matrix bv(){
    enum {initialCondition=0,terminalCondition=1}; 
    Matrix boundaryValues(2,n);
    boundaryValues[initialCondition] = Vector({1.0});
    boundaryValues[terminalCondition] = Vector({0.0});
    return boundaryValues;
  }
  
  //Linear interpolation
  Matrix initialTraj(){
    Matrix boundary = bv();
    Matrix init(N+1,m+n);
    
//     Vector control(m);
//     Vector ic = boundary[0];
//     ic.insert(ic.end(),control.begin(),control.end());
//     Vector tc = boundary[1];
//     tc.insert(tc.end(),control.begin(),control.end());
//     
//     for(int i=0;i<N+1;i++){
//       double lambda = double(i)/double(N);
//       init[i]=lambda*ic + (1.0-lambda)*tc; 
//     }
    
//     std::cout << " Initial guess is " << std::endl;
//     printMatrix(init);
    
    return init;
  }
  
  //~~~Derivatives of problem system attributes~~~///
  
  //Gradient of running cost
  Vector Dg(Vector x, Vector u) final {
    return Vector({0.0});
  }
  
  //Jacobian of the dynamics with respect to (x,u). Df = [df/dx, df/du]
  Matrix Df(Vector x, Vector u) final {
    Matrix Jacobian(n,n+m);
    Jacobian[0][0] = 0.0;
    Jacobian[0][1] = 1.0;
    return Jacobian;
  }
  
  //Jacobian of the state-control constraints
  Matrix Dh(Vector  x, Vector  u) final {
    std::size_t numConstraints = Vector (h(x,u)).size();
    Matrix dhdx(numConstraints,n+m);
    dhdx[0][0] = -1;//d/dx h_1(x,u)
    dhdx[0][1] = 0; //d/du h_1(x,u)
    dhdx[1][0] = 0; //d/dx h_2(x,u)
    dhdx[1][1] = -1;//d/du h_2(x,u)
    
    return dhdx;
  }
};

class Pendulum : public TPBVP {
  
public:
  Pendulum(unsigned int numSteps, double dt) : TPBVP(2, 1,numSteps,dt) {}
  
  //xdot = f(x,u) - dynamical model of a pendulum
  Vector  f(Vector  x, Vector  u) final {
    Vector  xdot(n);
    xdot[0] = x[1];
    xdot[1] = -sin(x[0])+u[0];
    
    return xdot;
  }

  //Jac = [df/dx df/du] - jacobian of f evaluated at (y[index],...,y[index+stateDim+controlDim-1])
Matrix Df(Vector x, Vector u) final {
    
    //create 2x3 matrix for jacobian
    Matrix Jacobian(n,n+m);
    
    Jacobian[0][0] = 0; 
    Jacobian[0][1] = 1; 
    Jacobian[0][2] = 0;
    Jacobian[1][0] = cos(x[0]);//first state 
    Jacobian[1][1] = 0;
    Jacobian[1][2] = 1;
    
    return Jacobian;
  }
  
  //no constraints
  Vector  h(Vector  x, Vector  u){
    return Vector (0);
  }
  Matrix Dh(Vector x, Vector u){
    return Matrix(0,0);
  }
  
  //Boundary values
  Matrix bv() final {
    enum {initialCondition=0,terminalCondition=1}; 
    Matrix boundaryValues(2,n);
    boundaryValues[initialCondition] = Vector({M_PI,0.0});
    boundaryValues[terminalCondition] = Vector({0.0,0.0});
    return boundaryValues;
  }
};

#endif