//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for weights.
// Similar to solve_for_weights, but not iterative
// Similar to w_equation

void linear::w_equation2( ) {

  cout << "Solving weight equation v2 " << endl;
  
  //  volumes( T );

  //  copy_weights( T ) ;

  VectorXd vol0  = field_to_vctr( sfield_list::vol0 ) ;
  VectorXd vol   = field_to_vctr( sfield_list::vol ) ;

  FT totV= vol.sum();

  int N = vol.size();

  FT meanV = totV/ FT( N );

  FT vol_sigma =  ( vol.array() - meanV ).square().sum() / FT( N );

  FT target_vol_val = meanV;

  cout << "  vol mean : " << meanV;
  cout << "  vol variance : " << vol_sigma << endl;

    //  FT target_v = simu.meanV();
  //  fill_Delta_DD();
  
  // given by volume departure from initial
  VectorXd Dvol  =  vol0 - vol ;

  // given by departure from mean.-
  // VectorXd Dvol  =  target_vol_val - vol.array() ;

  VectorXd Dw = Delta_solver.solve( Dvol );

  //  VectorXd w0  = field_to_vctr( sfield_list::w0 ) ;

  // vctr_to_field( w0 + Dw ,  sfield_list::w ) ;
  vctr_to_field( Dw ,  sfield_list::w ) ;

  return;
}




void linear::u_add_w_grad( const FT dt ) {

  VectorXd gradPx,gradPy;

  DD_times_sfield( sfield_list::w  ,  gradPx, gradPy);

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

