//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Iterative process to adjust weights so that
// the weighted Voronoi diagram has constant volumes

void linear::solve_for_weights_centroid( const FT dt ) {

  cout << "Making centroidal " << endl;
  
  volumes( T );
  copy_weights( T );
  //  VectorXd vol0 = field_to_vctr( sfield_list::vol0 ) ;

  //  VectorXd target_vol( vol ) ;
  //  target_vol.setConstant( target_vol_val );
  
  const int max_iter = 100;
  const FT threshold = 1e-16;
  const FT mixing = 1;
  
  for( int iter=0 ; iter< max_iter ; iter++) {

    VectorXd dd2  = field_to_vctr( sfield_list::dd2 ) ;

    // x component only . - 
    // VectorXd dd2, ddy;
    // vfield_to_vctrs( vfield_list::dd , dd2, ddy ) ;

    FT tot_dd2 = dd2.sum();

    int N = dd2.size();

    FT mean_dd2 = tot_dd2/ FT( N );
  
    FT dd2_sigma =  ( dd2.array() - mean_dd2 ).square().sum() / FT( N );

    // 0 target
    //    FT target_vol_val = meanV; 

    cout << "It:  " << iter;
    cout << "  dd2 mean : " << mean_dd2;
    cout << "  dd2 variance : " << dd2_sigma << endl;

    //  FT target_v = simu.meanV();

    VectorXd D_dd2 = - dd2.array()  ;

    fill_Delta_DD( dt );

    VectorXd Dw = EE_solver.solve( D_dd2 );

    //    copy_weights( T ) ;

    VectorXd w0  = field_to_vctr( sfield_list::w ) ;

    vctr_to_field( w0 + mixing * Dw ,  sfield_list::w ) ;
    //vctr_to_field( Dw ,  sfield_list::w ) ;

    volumes( T ); // ??

    move_weights( T );

    volumes( T );

    VectorXd dd2_0( dd2 );

    dd2  = field_to_vctr( sfield_list::dd2 ) ;

    FT diff = (dd2 - dd2_0 ).norm();

    cout << "dd2 diff: " << diff<< endl;
    
    tot_dd2 = dd2.sum();
    mean_dd2 = tot_dd2/ FT( N );  
    dd2_sigma =  ( dd2.array() - mean_dd2 ).square().sum() / FT( N );

    cout << "  dd2 mean : " << mean_dd2;
    cout << "  dd2 variance : " << dd2_sigma << endl;

    if( diff < threshold ) break;
  }

  return;
}
