#include <systems.h>

//This unit test simulates the dynamics of a pendulum with arbitrary control and tests that the optimization constraints are satisfied with the simulated result

// #include "trajectories.h"

//control signal for the simulation u(t)
Vector control(const double& t){
  Vector input({0.01});
  return input;
}

int main(){
  std::cout << "---Testing numerical approximation constraint jacobian---" << std::endl;
  
  //Numerical integration parameters
  double dt = 0.5;
  double t0 = 0.0;
  int num_steps = 3;
  FirstOrderSystem sys(num_steps,dt);
  
  //Initial state, control, and time 
  Vector x({1.0});
  double t = 0;
  Vector u = control(t0);
  
  //The y is the the traj and control stacked up into a vector with params at the end: y=(x_0 u_0 x_1 u_1 ... dt)
  Vector y({x});
  y.push_back(u[0]);
  
  //Do the simulation with Euler's method
  for(int i=1;i<num_steps+1;i++){
    
    //one time-step
    u = control(t);
    x = x + dt*sys.f(x,u);
    t+=dt;
    
    //Append the new state and control to the vector y
    for(int j=0;j<sys.state_dimension();j++){
      y.push_back(x[j]);
    }
    for(int j=0;j<sys.control_dimension();j++){
      y.push_back(u[0]);
    }
  }
  
  //Append the simulation parameters to y
  y.push_back(dt);
  
  //Evaluate the constraint for the variable y (should be zero)
  std::cout << "trajectory and control " << std::endl;
  printVector(y);
  
  std::vector<double> residual = sys.phi(y);
  std::cout << "dynamic feasibility residual " << std::endl;
  printVector(residual);
  
  Matrix dgdy = sys.Dphi(y);
  std::cout << "first variation of dynamic constraint" << std::endl;
  printMatrix(dgdy);
  
  std::vector<double> signal_constraints = sys.H(y);
  std::cout << "state and input constraints " << std::endl;
  printVector(signal_constraints);
  
  Matrix constraint_jacobian = sys.DH(y);
  std::cout<< "first variataion of state-control constraint" << std::endl;
  printMatrix(constraint_jacobian);
  
  double cost = sys.J(y);
  std::cout << " cost is " << cost << std::endl;
  
  Vector gradCost = sys.DJ(y);
  std::cout << " gradient of cost is " << std::endl;
  printVector(gradCost);
  
  
}