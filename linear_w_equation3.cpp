//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for weights.
// Similar to solve_for_weights, but not iterative
// Similar to w_equation

void linear::w_equation3( ) {

  cout << "Solving weight equation v3 " << endl;
  
  //  volumes( T );

  //  copy_weights( T ) ;

  VectorXd I  = field_to_vctr( sfield_list::I ) ;
  VectorXd I0  = field_to_vctr( sfield_list::I0 ) ;

  int N = I.size();

  VectorXd DI = I0 - I ;

  FT DI_sigma =  DI.array().square().sum() / N;
  FT I_mean  =  I.array().sum() / N ;

  cout << " w field  "
       << " rel DI std dev: " << sqrt( DI_sigma ) / I_mean
       << endl;

  VectorXd Dw = GG_solver.solve( DI );

  VectorXd w0  = field_to_vctr( sfield_list::w ) ;

  vctr_to_field( w0 + Dw ,  sfield_list::w ) ;

  //  vctr_to_field( Dw ,  sfield_list::w ) ;

  return;
}



