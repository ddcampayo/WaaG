//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for s field

void linear::p_equation(const FT dt ) {

  cout << "Solving s equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly
  

  //  volumes( T );

  VectorXd I  = field_to_vctr( sfield_list::I ) ;

  //FT target_vol_val =  simu.meanV() ;

  //  FT target_vol_val = vol.array().sum() / FT( vol.size() );
  //VectorXd Dvol = vol.array() - target_vol_val  ;

  VectorXd I0  = field_to_vctr( sfield_list::vol0 ) ;
  VectorXd DI = I.array() - I0.array()  ;

  FT DI_sigma =  DI.array().square().sum() ; // / FT( vol.size() );
  FT DI_mean  =  I.array().square().sum() ; // / FT( vol.size() );

  cout << " s field  "
       << " rel DI std dev: " << sqrt( DI_sigma / DI_mean )
       << endl;

  VectorXd Ds  =  NN_solver.solve( DI );

  VectorXd s0  = field_to_vctr( sfield_list::s ) ;

  vctr_to_field( s0 + Ds / ( ddt * ddt) , sfield_list::s  ) ;

  //  vctr_to_field( vol , sfield_list::vol0 );

  return;
}



// Not a "gradient" at all, but named so in parallel with the p
// counterpart
void linear::u_add_s_grad( const FT dt ) {

  VectorXd gradsx,gradsy;

  MM_times_sfield( sfield_list::s  ,  gradsx, gradsy);

  VectorXd I  = field_to_vctr( sfield_list::I );
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  U_x = Ustar_x.array() - ddt * gradsx.array() / vol.array()  ;
  U_y = Ustar_y.array() - ddt * gradsy.array() / vol.array() ;

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}
