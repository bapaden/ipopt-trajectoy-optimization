//This unit test checks that the DynamicalSystem class is working properly

#include <systems.h>

int main(){
  std::cout << "---Testing pendulum dynamics and jacobian calculation---" << std::endl;
  
  Pendulum overratedExample(10,0.1);
  
  State x({M_PI, 1.0});
  Control u({1.0});
  
  DecisionVar y;
  y.insert(y.end(),x.begin(),x.end());
  y.insert(y.end(),u.begin(),u.end());
  
  std::cout << "(x,u)=";
  printVector(y);
  
  State xdot = overratedExample.vectorField(x,u);
  std::cout << "xdot=" << std::endl;
  printVector(xdot);
  std::cout << "\ndf/dx=" <<std::endl;
  Matrix jacobian = overratedExample.systemJacobian(y,0);
  printMatrix(jacobian);//Should be [0,1,0;-1,0,1] 
  
}