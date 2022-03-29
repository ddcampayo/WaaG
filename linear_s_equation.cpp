//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for s field

void linear::s_equation(const FT dt ) {

  cout << "Solving s equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly  

  //  volumes( T );

  VectorXd I  = field_to_vctr( sfield_list::I ) ;
  VectorXd I0  = field_to_vctr( sfield_list::I0 ) ;

  int N = I.size();

  VectorXd DI = I.array() - I0.array()  ;

  FT DI_sigma =  DI.array().square().sum() / N;
  FT I_mean  =  I.array().sum() / N ;

  cout << " s field  "
       << " rel DI std dev: " << sqrt( DI_sigma ) / I_mean
       << endl;

  VectorXd Ds  =  NN_solver.solve( DI );

  VectorXd s0  = field_to_vctr( sfield_list::s ) ;

  //  vctr_to_field( s0 + Ds / ( ddt * ddt) , sfield_list::s  ) ;
  vctr_to_field( Ds / ( ddt * ddt) , sfield_list::s  ) ;

  //  vctr_to_field( vol , sfield_list::vol0 );

  return;
}


void linear::s_equation_p(const FT dt ) {

  cout << "Solving s equation " << endl;
  
  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd MMUstar  =  MM_scalar_vfield( vfield_list::Ustar );
  
  VectorXd p  = field_to_vctr( sfield_list::p );

  VectorXd NLp = NL * p;

  VectorXd s =  NN_solver.solve( - MMUstar / dt + NLp );
  
  vctr_to_field( s ,  sfield_list::s ) ;

  return;
}




