// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-09

#ifndef OPT_INTERFACE_HPP
#define OPT_INTERFACE_HPP

#include <IpTNLP.hpp>
#include <ipopt_adaptor.h>
#include <limits>

// using namespace Ipopt;

/** C++ Example NLP for interfacing a problem with IPOPT.
 *  HS071_NLP implements a C++ example of problem 71 of the
 *  Hock-Schittkowski test suite. This example is designed to go
 *  along with the tutorial document and show how to interface
 *  with IPOPT through the TNLP interface. 
 *
 * Problem hs071 looks like this
 *
 *     min   x1*x4*(x1 + x2 + x3)  +  x3
 *     s.t.  x1*x2*x3*x4                   >=  25
 *           x1**2 + x2**2 + x3**2 + x4**2  =  40
 *           1 <=  x1,x2,x3,x4  <= 5
 *
 *     Starting point:
 *        x = (1, 5, 5, 1)
 *
 *     Optimal solution:
 *        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
 *
 *
 */
class TrajOpt : public Ipopt::TNLP
{
  TPBVP* myProb;
  tester t;
  int n;//number of states
  int m;//number of controls
  int N;//number of time steps
public:
  /** default constructor */
  TrajOpt(TPBVP* myProb_);
  
  /** default destructor */
  virtual ~TrajOpt();
  
  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Ipopt::Index& n, 
                            Ipopt::Index& m, 
                            Ipopt::Index& nnz_jac_g,
                            Ipopt::Index& nnz_h_lag, 
                            Ipopt::TNLP::IndexStyleEnum& index_style);
  
  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Ipopt::Index n,
                               Ipopt::Number* x_l, 
                               Ipopt::Number* x_u,
                               Ipopt::Index m, 
                               Ipopt::Number* g_l, 
                               Ipopt::Number* g_u);
  
  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Ipopt::Index n, 
                                  bool init_x, 
                                  Ipopt::Number* x,
                                  bool init_z, 
                                  Ipopt::Number* z_L, 
                                  Ipopt::Number* z_U,
                                  Ipopt::Index m, 
                                  bool init_lambda,
                                  Ipopt::Number* lambda);
  
  /** Method to return the objective value */
  virtual bool eval_f(Ipopt::Index n, 
                      const Ipopt::Number* x, 
                      bool new_x, 
                      Ipopt::Number& obj_value);
  
  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Ipopt::Index n, 
                           const Ipopt::Number* x, 
                           bool new_x, 
                           Ipopt::Number* grad_f);
  
  /** Method to return the constraint residuals */
  virtual bool eval_g(Ipopt::Index n, 
                      const Ipopt::Number* x, 
                      bool new_x, 
                      Ipopt::Index m, 
                      Ipopt::Number* g);
  
  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Ipopt::Index n, 
                          const Ipopt::Number* x, 
                          bool new_x,
                          Ipopt::Index m, 
                          Ipopt::Index nele_jac, 
                          Ipopt::Index* iRow, 
                          Ipopt::Index *jCol,
                          Ipopt::Number* values);
  
  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Ipopt::Index n, 
                      const Ipopt::Number* x, 
                      bool new_x,
                      Ipopt::Number obj_factor, 
                      Ipopt::Index m, 
                      const Ipopt::Number* lambda,
                      bool new_lambda, 
                      Ipopt::Index nele_hess, 
                      Ipopt::Index* iRow,
                      Ipopt::Index* jCol, 
                      Ipopt::Number* values);
  
  //@}
  
  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(Ipopt::SolverReturn status,
                                 Ipopt::Index n, 
                                 const Ipopt::Number* x, 
                                 const Ipopt::Number* z_L, 
                                 const Ipopt::Number* z_U,
                                 Ipopt::Index m, 
                                 const Ipopt::Number* g, 
                                 const Ipopt::Number* lambda,
                                 Ipopt::Number obj_value,
                                 const Ipopt::IpoptData* ip_data,
                                 Ipopt::IpoptCalculatedQuantities* ip_cq);
  //@}
  
private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  TrajOpt();
  TrajOpt(const TrajOpt&);
  TrajOpt& operator=(const TrajOpt&);
  //@}
};

// constructor
TrajOpt::TrajOpt(TPBVP* myProb_){
  myProb=myProb_;
  n=myProb->stateDim();
  m=myProb->controlDim();
  N=myProb->numSteps();
}

//destructor
TrajOpt::~TrajOpt()
{}

// returns the size of the problem
bool TrajOpt::get_nlp_info(Ipopt::Index& numVars, 
                           Ipopt::Index& numCons, 
                           Ipopt::Index& nnz_jac_g,
                           Ipopt::Index& nnz_h_lag, 
                           Ipopt::TNLP::IndexStyleEnum& index_style)
{
  // The total number of decision variables
  numVars = (N+1)*(n+m)+1;
  
  // Constraints for the dynamics and state-control constraint
  numCons = N*n+(myProb->H(Vector(numVars))).size();
  printf("numVars in get_nlp_info: %d \n",numVars);
  printf("numCons in get_nlp_info: %d \n",numCons);
  
  
  //Number of nonzero entries in dgdz
  nnz_jac_g = numVars*numCons;
  printf("nnz_jac_g in get_nlp_info: %d ",nnz_jac_g);
  
  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = numVars*numVars;
  
  // use the C style indexing (0-based)
  index_style = Ipopt::TNLP::C_STYLE;
  
  printf("numVars in get_nlp_info: %d \n",numVars);
  printf("numCons in get_nlp_info: %d \n",numCons);
  
  return true;
}

// returns the variable bounds
bool TrajOpt::get_bounds_info(Ipopt::Index numVars, 
                              Ipopt::Number* y_l, 
                              Ipopt::Number* y_u,
                              Ipopt::Index numCons, 
                              Ipopt::Number* g_l, 
                              Ipopt::Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(numVars == (N+1)*(n+m)+1);
  assert(numCons == N*n+(myProb->H(Vector(numVars))).size());
  
  printf("numVars in get_bounds_info: %d \n",numVars);
  printf("numCons in get_bounds_info: %d \n",numCons);
  
  
  for(int i=0;i<numVars;i++){
    y_l[i]=std::numeric_limits<double>::min();//+infinity
    y_u[i]=std::numeric_limits<double>::max();//-infinity
  }
  
  // Initial condition
  printf("\ninitial condition");
  for(int i=0; i<n; i++) {
    std::cout << std::endl << (myProb->bv())[0][i] << ",";
    y_l[i] = y_u[i] = (myProb->bv())[0][i];
  }
  
  // Terminal condition
  printf("\nterminal condition");
  for(int i=0;i<n; i++){
    std::cout <<std::endl << (myProb->bv())[1][i] << "," ;
    y_l[N*(n+m)+i] = y_u[N*(n+m)+i] = (myProb->bv())[1][i];
  }
  
  // Equality constraints for dynamic feasibility
  for(int i=0;i<n*N;i++){
    g_l[i]=0.0;
    g_u[i]=0.0;
  }
  
  //constrants for h(x,u)<=0 
  for(int i=N*n;i<numCons;i++){
    g_l[i] = -2.0e9;//-std::numeric_limits<double>::max();
    g_u[i] = 0.0;
  }
  
  for(int i=0;i<numCons;i++){
    std::cout << std::endl << "(g_l[" << i << "],g_u[" << i << "])=(" << g_l[i] << "," << g_u[i] << std::endl;  
  }
  return true;
}

// returns the initial point for the problem
bool TrajOpt::get_starting_point(Ipopt::Index numVars, 
                                 bool init_x, 
                                 Ipopt::Number* y,
                                 bool init_z, 
                                 Ipopt::Number* z_L, 
                                 Ipopt::Number* z_U,
                                 Ipopt::Index numCons, 
                                 bool init_lambda,
                                 Ipopt::Number* lambda)
{
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);
  
  // initialize to the given starting point
  Vector init = myProb->initialGuess();
  
  //Write into y
  for(int i=0;i<init.size();i++){
    y[i]=init[i];
  }
  
  return true;
}

// returns the value of the objective function
bool TrajOpt::eval_f(Ipopt::Index numVars, 
                     const Ipopt::Number* y_, 
                     bool new_x, Ipopt::Number& obj_value)
{
  //   dataVec.insert(dataVec.end(), &dataArray[0], &dataArray[dataArraySize]);
  Vector y(&y_[0],&y_[numVars]);
  //   y.insert(y.end(), &y_[0], &y_[numVars-1]); 
  obj_value = myProb->J(y);
  
  std::cout << " Cost: " << obj_value << std::endl;
  std::cout << " Decisio Var " << std::endl;
  printVector(y);
  
  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool TrajOpt::eval_grad_f(Ipopt::Index numVars, const Ipopt::Number* y_, bool new_y, Ipopt::Number* grad_f){
  
  Vector y(&y_[0], &y_[numVars]);
  //   y.insert(y.end(), &y_[0], &y_[numVars-1]); 
  //   std::cout << "what is y" << std::endl;
  //   printVector(y);
  //Never tried this... new cpp standard needed probably 
  Vector DJ = myProb->DJ(y);
  grad_f = &DJ[0];
  std::cout << " Cost gradient: " << std::endl;
  printVector(DJ);
  
  
  return true;
}

// return the value of the constraints: g(x)
bool TrajOpt::eval_g(Ipopt::Index numVars, 
                     const Ipopt::Number* y_, 
                     bool new_y, 
                     Ipopt::Index numCons, 
                     Ipopt::Number* g)
{
  //   std::cout << "numVars " << numVars << std::endl;
  Vector y(&y_[0], &y_[numVars]);
  //   y.insert(y.begin(), &y_[0], &y_[numVars]); 
  //   std::cout << "y_ " << y_[0] << "," << y_[1] << "," << y_[2] << "," << y_[3] << "," << y_[4] << std::endl;
  
  //Constraints related to dynamic feasibility
  //   std::cout << "A: ";
  Vector A = myProb->phi(y);
  //Point-wise constraints on state and control
  Vector B = myProb->H(y);
  A.insert(A.end(), B.begin(), B.end());
  //   std::cout << "\n\nDecision Variable: " << std::endl;
  //   printVector(y);
  std::cout << " Constraints: " << std::endl;
  printVector(A);
  //   printf("\n\n");
  
  g = &A[0];
  
  return true;
}

// return the structure or values of the jacobian
bool TrajOpt::eval_jac_g(Ipopt::Index numVars, 
                         const Ipopt::Number* y_,
                         bool new_y,
                         Ipopt::Index numCons, 
                         Ipopt::Index nele_jac, 
                         Ipopt::Index* iRow, 
                         Ipopt::Index *jCol,
                         Ipopt::Number* values)
{
  if (values == NULL) {
    //     std::cout << "Getting sparsity pattern" << std::endl;
    int k=0;
    // this particular jacobian is dense
    for(int i=0;i<5;i++){//TODO not 5
      for(int j=0;j<(N+1)*(n+m)+1;j++){
        iRow[k]=i; jCol[k]=j; 
        //         std::cout << "index: " << k << "  (i,j)=(" << i << "," << j << ")" << std::endl;
        k++;
      }
    }
  }
  else {
    // return the values of the jacobian of the constraints
    
    Vector y(&y_[0], &y_[numVars]);
    //     y.insert(y.end(), &y_[0], &y_[numVars-1]); 
    
    
    nele_jac=numCons*numVars;
    
    Matrix jac = myProb->Dphi(y);
    
    jac.append_bottom( myProb->DH(y) );
    assert(numVars == jac[0].size() );
    
    int k=0;
    std::cout << " Jacobian:" << jac.size() << std::endl; 
    printMatrix(jac);
    for(int i=0;i<jac[0].size();i++){
      for(int j=0;j<jac.size();j++){
        values[k]=jac[i][j];
        
        //         std::cout << "value: " << values[k] << " index: " << k << "Position: (i,j)=(" << i << "," << j << ")" << std::endl;
        k++;
      }
    }
  }
  
  return true;
}

//return the structure or values of the hessian
bool TrajOpt::eval_h(Ipopt::Index numVars, 
                     const Ipopt::Number* y_, 
                     bool new_y,
                     Ipopt::Number obj_factor, 
                     Ipopt::Index numCons, 
                     const Ipopt::Number* lambda,
                     bool new_lambda, 
                     Ipopt::Index nele_hess, 
                     Ipopt::Index* iRow,
                     Ipopt::Index* jCol, 
                     Ipopt::Number* values)
{
  //   std::cout << "eval_h" << std::endl;
  
  return false;//Too lazy to compute this. Just use bfgs
}

void TrajOpt::finalize_solution(Ipopt::SolverReturn status,
                                Ipopt::Index n, 
                                const Ipopt::Number* x, 
                                const Ipopt::Number* z_L, 
                                const Ipopt::Number* z_U,
                                Ipopt::Index m,
                                const Ipopt::Number* g, 
                                const Ipopt::Number* lambda,
                                Ipopt::Number obj_value,
                                const Ipopt::IpoptData* ip_data,
                                Ipopt::IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.
  
  // For this example, we write the solution to the console
  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Ipopt::Index i=0; i<n; i++) {
    std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }
  
  
  
  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  for (Ipopt::Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  }
  for (Ipopt::Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }
  
  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;
  
  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Ipopt::Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }
}

#endif
