#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure from the s field

void linear::p_equation_from_s( const FT dt ) {

  cout << "Solving pressure equation " << endl;
  
  VectorXd s  = field_to_vctr( sfield_list::s );

  VectorXd gamma_s  = GG.transpose() * s;

  VectorXd p;

  p =  Delta_solver.solve( gamma_s );

  //  vctr_to_field( p / (dt*dt)  ,  sfield_list::p ) ;
  vctr_to_field( p ,  sfield_list::p ) ;

  return;
}



// Solve for s from the pressure field

void linear::s_equation_from_p( const FT dt ) {

  cout << "Solving s equation " << endl;
  
  VectorXd p  = field_to_vctr( sfield_list::p );

  VectorXd Delta_p  = Delta.transpose() * p; // Delta is symmetric anyway ...

  VectorXd s;

  s =  GG_solver.solve( Delta_p );

  // offset:
  //  FT s_min = s.minCoeff(); 
  //  s = s.array() - s_min;

  vctr_to_field( s ,  sfield_list::s ) ;

  return;
}

