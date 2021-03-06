//This unit test simulates the dynamics of a pendulum with arbitrary control and tests that the optimization constraints are satisfied with the simulated result

#include <systems.h>

//control signal for the simulation u(t)
Control control(const double& t){
  Control input({0.0});
  return input;
}

int main(){
  std::cout << "---Testing dynamic constraint evaluation---" << std::endl;
  
  //Numerical integration parameters
  double dt = 0.03;
  double t0 = 0.0;
  int num_steps = 200;
  Pendulum myPend(num_steps,dt);
  
  //Initial state, control, and time 
  State x({M_PI/2.0,0.0});
  double t = 0;
  Control u = control(t0);

  //The y is the the traj and control stacked up into a vector with params at the end: y=(x_0 u_0 x_1 u_1 ... dt)
  DecisionVar y({x});
  y.push_back(u[0]);
  
  //Do the simulation with Euler's method
  for(int i=1;i<num_steps+1;i++){
    
    //one time-step
    u = control(t);
    x = x + dt*myPend.vectorField(x,u);
    t+=dt;
    
    //Append the new state and control to the vector y
    for(int j=0;j<myPend.state_dimension();j++){
      y.push_back(x[j]);
    }
    for(int j=0;j<myPend.control_dimension();j++){
      y.push_back(u[0]);
    }
  }
  
  //Append the simulation parameters to y
  y.push_back(dt);

  //Evaluate the constraint for the variable y (should be zero)
  DecisionVar residual = myPend.constraintResidual(y);
  bool passed = true;
  for(auto x : residual){
    if(std::fabs(x)>1e-12){
      passed = false;
      std::cout << "\n\nUNIT TEST FAILED!\n\n" << std::endl;
      break;
    }
  }
  if(passed==true){std::cout << "\n\nUNIT TEST PASSED!\n\n" << std::endl;}
}