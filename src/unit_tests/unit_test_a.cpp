//This unit test checks that the DynamicalSystem class is working properly

#include "trajectories.h"

int main(){
  std::cout << "---Testing pendulum dynamics and jacobian calculation---" << std::endl;
  
  Pendulum overratedExample(10,0.1);
  
  State x({M_PI, 1.0});
  Control u({1.0});
  
  State xdot = overratedExample.vectorField(x,u);
  Matrix jacobian = overratedExample.systemJacobian(x,u);
  
  printVector(xdot);//Should be (1,1)
  printMatrix(jacobian);//Should be [0,1,0;-1,0,1]  
}