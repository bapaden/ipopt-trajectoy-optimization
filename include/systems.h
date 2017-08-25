#ifndef SYSTEMS_H
#define SYSTEMS_H

#include <trajectories.h>

class FirstOrderSystem : public DynamicalSystem {
public:
  FirstOrderSystem(unsigned int numSteps, double dt) : DynamicalSystem(1,1,numSteps,dt){}
  
  std::vector<double> vectorField(std::vector<double> x, std::vector<double> u) final {
    std::vector<double> xdot(stateDim);
    xdot[0] = -x[0]+u[0];
    return xdot;
  }
  
  Matrix systemJacobian(DecisionVar y, unsigned int index) final {
    Matrix Jacobian(stateDim,stateDim+controlDim);
    Jacobian[0][0] = -1.0;
    Jacobian[0][1] = 1.0;
    return Jacobian;
  }
  
  std::vector<double> signalConstraint(std::vector<double> x, std::vector<double> u) final {
    assert(x.size()==n and u.size()==m);
    std::vector<double> h;
    h.push_back(0.25-x[0]);
//     double h0=0.25-x[0];
//     double h1=-0.5-u[0];
//     double h2=u[0]-0.5;
    return h;
  }
};

class Pendulum : public DynamicalSystem {
  
public:
  Pendulum(unsigned int numSteps, double dt) : DynamicalSystem(2, 1,numSteps,dt) {}
  
  //xdot = f(x,u) - dynamical model of a pendulum
  std::vector<double> vectorField(std::vector<double> x, std::vector<double> u) final {
    std::vector<double> xdot(stateDim);
    xdot[0] = x[1];
    xdot[1] = -sin(x[0])+u[0];
    
    return xdot;
  }

  //Jac = [df/dx df/du] - jacobian of f evaluated at (y[index],...,y[index+stateDim+controlDim-1])
Matrix systemJacobian(DecisionVar y, unsigned int index) final {
    
    //create 2x3 matrix for jacobian
    Matrix Jacobian(stateDim,stateDim+controlDim);
    
    Jacobian[0][0] = 0; 
    Jacobian[0][1] = 1; 
    Jacobian[0][2] = 0;
    Jacobian[1][0] = cos(y[index]);//first state 
    Jacobian[1][1] = 0;
    Jacobian[1][2] = 1;
    
    return Jacobian;
  }
  
  //no constraints
  std::vector<double> signalConstraint(std::vector<double> x, std::vector<double> u){
    return std::vector<double>(0);
  }
};


#endif