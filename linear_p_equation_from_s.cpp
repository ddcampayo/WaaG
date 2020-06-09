#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure from the s field

void linear::p_equation_from_s( ) {

  cout << "Solving pressure equation " << endl;
  
  VectorXd s  = field_to_vctr( sfield_list::s );

  VectorXd gamma_s  = GG.transpose() * s;

  VectorXd p;

  p =  Delta_solver.solve( gamma_s );

  vctr_to_field( p  ,  sfield_list::p ) ;

  return;
}



// Solve for s from the pressure field

void linear::s_equation_from_p( ) {

  cout << "Solving pressure equation " << endl;
  
  VectorXd p  = field_to_vctr( sfield_list::p );

  VectorXd Delta_p  = Delta.transpose() * p;

  VectorXd s;

  s =  GG_solver.solve( Delta_p );

  vctr_to_field( s  ,  sfield_list::s ) ;

  return;
}

