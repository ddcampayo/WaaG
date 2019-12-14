#ifndef _LINEAR_H_
#define _LINEAR_H_

#include"pParticles.h"

#ifdef CHOLMOD
 #include <Eigen/CholmodSupport>
#else
// #include <Eigen/IterativeLinearSolvers>
 #include <Eigen/SparseCholesky>
//#include <Eigen/SparseQR>
#endif

#include <iostream>
#include <vector>

#include <unsupported/Eigen/SparseExtra>

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
//using Eigen::ConjugateGradient;

//const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];");

class linear {
 public:  // constructor
 linear(Triangulation& TT) : T(TT) {}

  void fill_Delta();
  void fill_Delta_DD( const FT dt = 0);

  void solve_for_weights(  const FT dt = 0 );
  void solve_for_weights_centroid(  const FT dt = 0 );
  void solve_for_moments( );
  void w_equation( );
  void w_equation2( );
  void w_equation3( );
  void p_equation(const FT dt );
  void s_equation(const FT dt );
  void u_star( void );
  void reset_p( void );
  void reset_s( void );
  void u_add_press_grad( const FT dt ) ;
  void u_add_s_grad( const FT dt ) ;
  void u_add_w_grad( const FT dt ) ;
  void u_add_grads( const FT dt ) ;
  void u_add_press_grad_MM_w( const FT dt );

  void DD_scalar_vfield(const vfield_list::take from , const sfield_list::take to );  
  VectorXd DD_scalar_vfield(const vfield_list::take from );
  void DD_times_sfield(const sfield_list::take from ,
		       VectorXd& Dx,VectorXd& Dy);
  void MM_times_sfield(const sfield_list::take from ,
		       VectorXd& Dx,VectorXd& Dy);
  void copy(const sfield_list::take from, sfield_list::take to  );

  void dd2_stats( void ) ;

  
private:

  
  Triangulation& T; // Its triangulation

  typedef   SparseMatrix<double>  SpMat;
  typedef Eigen::Triplet<double> triplet;

  SpMat Delta;
  SpMat GG;
  SpMat DDx, DDy;
  SpMat LL;
  SpMat MMx, MMy;
  SpMat NN;

  SpMat EE;

  
  VectorXd field_to_vctr(const sfield_list::take sf );
  void vctr_to_field(const VectorXd& vv, const sfield_list::take sf  );
  void vfield_to_vctrs(const vfield_list::take vf , VectorXd& vx, VectorXd& vy ) ;
  void vctrs_to_vfield(const VectorXd& vx, const VectorXd& vy , const vfield_list::take vf ) ;

#define DIRECT_SOLVER

#ifdef DIRECT_SOLVER
#ifdef CHOLMOD
  Eigen::CholmodSupernodalLLT<SpMat> Delta_solver;
#else
  Eigen::SimplicialLDLT<SpMat> Delta_solver;
  Eigen::SimplicialLDLT<SpMat> EE_solver;
  Eigen::SimplicialLDLT<SpMat> GG_solver;
  Eigen::SimplicialLDLT<SpMat> LL_solver;
  Eigen::SimplicialLDLT<SpMat> NN_solver;
#endif
#else
  Eigen::BiCGSTAB<SpMat> solver_stiffp1;
#endif


};


#endif
