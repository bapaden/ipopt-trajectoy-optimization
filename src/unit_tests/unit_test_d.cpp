#include <systems.h>

//This unit test simulates the dynamics of a pendulum with arbitrary control and tests that the optimization constraints are satisfied with the simulated result

// #include "trajectories.h"

//control signal for the simulation u(t)
Vector control(const double& t){
  Vector input({0.01});
  return input;
}

int main(){
  
  //Numerical integration parameters
  double dt = 1;
  double t0 = 0.0;
  int num_steps = 1;
  FirstOrderSystem sys(num_steps,dt);
  
  //Initial state, control, and time 
  Vector x({1.0});
  double t = 0;
  Vector u = control(t0);
  

  Vector y({1.0, -2.0, 0.0, 0.0, 0.5});
  
  double cost(sys.J(y));
  std::cout << "Cost is " << std::endl;
  std::cout << "J=" << cost << std::endl;
  
  Vector grad_f(sys.DJ(y));
  std::cout << "Gradient of cost: " << std::endl;
  printVector(grad_f);
  
  Vector g(sys.H(y));
  std::cout << "ipopt constraint functions " << std::endl << "g=" << std::endl;
  printVector(g);
  
  Matrix grad_g(sys.DH(y));
  std::cout << "constraint gradient: " << std::endl;
  printMatrix(grad_g);
  
  
}