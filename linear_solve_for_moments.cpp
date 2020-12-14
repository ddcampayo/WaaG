//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Iterative process to adjust moments so that
// the weighted Voronoi diagram has cells with
// fixed moments of inertia (aka second moments of mass)

// TODO: badly named

void linear::solve_for_moments( ) {

  cout << "Equalizing moments " << endl;
  
  volumes( T );
  copy_weights( T );
  VectorXd I0 = field_to_vctr( sfield_list::I0 ) ;

  //  VectorXd target_vol( vol ) ;
  //  target_vol.setConstant( target_vol_val );
  
  const int max_iter = 20;
  const FT threshold = 1e-6;
  const FT mixing = 0.5; // 1: only new iter; 0: only old

  int iter=0;

  for( ; iter< max_iter ; iter++) {

    VectorXd I  = field_to_vctr( sfield_list::I ) ;

    FT diff = (I - I0 ).norm() / I.norm(); // relative!

    cout << "I rel diff: " << diff<< endl;

    cout << "  solving for  ws, iter : " << iter << endl;
    
    // FT totI    = I.sum();
    // int N      = I.size();
    // FT meanI   = totI / FT( N );
    // FT I_sigma =  ( I.array() - meanI ).square().sum() / FT( N );

    // FT target_I = meanI;
    // VectorXd DI  =  mixing * ( target_I - I.array() );

    VectorXd DI = mixing * ( I0  - I ) ;

    fill_Delta_DD();

    VectorXd Dw = GG_solver.solve( DI );

    //    copy_weights( T ) ;

    VectorXd w0  = field_to_vctr( sfield_list::w ) ;

    vctr_to_field( w0 + Dw ,  sfield_list::w ) ;

    volumes( T ); // ??

    move_weights( T );

    volumes( T );

    //    VectorXd I0( I );

    //    I  = field_to_vctr( sfield_list::I ) ;


    if( diff < threshold ) break;
  }

  cout << "Equalized moments after " << iter << " iterations" << endl;

  return;
}


void linear::u_add_press_grad_MM_w( const FT dt ) {

  VectorXd gradPx,gradPy;

  DD_times_sfield( sfield_list::p  ,  gradPx, gradPy);

  VectorXd w = field_to_vctr( sfield_list::w );

  VectorXd MMw_x, MMw_y;

  MM_times_sfield( sfield_list::w  ,  MMw_x, MMw_y ) ;
  
    VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  U_x = Ustar_x.array() - ddt * ( MMw_x.array() +  gradPx.array() ) / vol.array()  ;
  U_y = Ustar_y.array() - ddt * ( MMw_y.array() +  gradPy.array() ) / vol.array()  ;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}

