#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure from the s field

void linear::p_equation_from_s( ) {

  cout << "Solving pressure equation " << endl;
  
  VectorXd s  = field_to_vctr( sfield_list::s );

  VectorXd gamma_s  = GG.transpose() * s;

  VectorXd p;

  p =  Delta_solver.solve( -gamma_s );

  vctr_to_field( p  ,  sfield_list::p ) ;

  return;
}

