//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure
// The famous pressure Poisson equation

void linear::p_equation(const FT dt ) {

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly
  
  // A
  // Approximate Laplacian ~ Delta 
  //  VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );
  // VectorXd p =  Delta_solver.solve( divUstar );
  // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  // vctr_to_field( -0.5 * p / ddt ,  sfield_list::p ) ;

  // B
  //  Laplacian as div of grad :
  // VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );
  //  VectorXd p =  LL_solver.solve( divUstar );
  //  vctr_to_field( p / ddt ,  sfield_list::p ) ;

  // C
  // As B, but Dvol source. This is an _iterative_ procedure,
  // yielding a pressure change

  //  volumes( T );

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  //FT target_vol_val =  simu.meanV() ;

  //  FT target_vol_val = vol.array().sum() / FT( vol.size() );
  //VectorXd Dvol = vol.array() - target_vol_val  ;

  VectorXd vol0  = field_to_vctr( sfield_list::vol0 ) ;

  VectorXd Dvol = vol.array() - vol0.array()  ;

  int N = vol.size();

  FT Dvol_sigma =  Dvol.array().square().sum() / N ; // / FT( vol.size() );
  FT Dvol_mean  =  vol.array().sum() / N ; // / FT( vol.size() );

  cout << "Pressure  "
       << " rel Dvol std dev: " << sqrt( Dvol_sigma ) / Dvol_mean 
       << endl;

  VectorXd Dp  =  LL_solver.solve( Dvol );

  VectorXd p0  = field_to_vctr( sfield_list::p ) ;
  
  //vctr_to_field( p0 + Dp / ( ddt * ddt) , sfield_list::p  ) ;

  vctr_to_field(  Dp / ( ddt * ddt) , sfield_list::p  ) ;

  //  vctr_to_field( vol , sfield_list::vol0 );

  return;
}




void linear::u_add_press_grad( const FT dt ) {

  VectorXd gradPx,gradPy;

  DD_times_sfield( sfield_list::p  ,  gradPx, gradPy);

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  U_x = Ustar_x.array() - ddt * gradPx.array() / vol.array()  ;
  U_y = Ustar_y.array() - ddt * gradPy.array() / vol.array() ;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}



void linear::u_add_grads( const FT dt ) {

  VectorXd gradPx,gradPy;
  VectorXd gradsx,gradsy;

  DD_times_sfield( sfield_list::p  ,  gradPx, gradPy);

  MM_times_sfield( sfield_list::s  ,  gradsx, gradsy);

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine

  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly
  VectorXd grad_x = gradPx + gradsx ;
  VectorXd grad_y = gradPy + gradsy ;
  
  U_x = Ustar_x.array() - ddt * grad_x.array() / vol.array()  ;
  U_y = Ustar_y.array() - ddt * grad_y.array() / vol.array() ;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}
