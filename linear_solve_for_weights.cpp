//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Iterative process to adjust weights so that
// the weighted Voronoi diagram has equal volumes

void linear::solve_for_weights( ) {

  cout << "Equalizing volumes " << endl;
  
  volumes( T );
  copy_weights( T );
  //  VectorXd vol0 = field_to_vctr( sfield_list::vol0 ) ;

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  FT totV= vol.sum();

  int N = vol.size();

  FT meanV = totV/ FT( N );
  
  FT vol_sigma =  ( vol.array() - meanV ).square().sum() / FT( N );

  FT target_vol_val = meanV;

  //  VectorXd target_vol( vol ) ;
  //  target_vol.setConstant( target_vol_val );
  
  const int max_iter = 100;
  const FT threshold = 1e-16;
  const FT mixing = 1;
  
  for( int iter=0 ; iter< max_iter ; iter++) {
    cout << "It:  " << iter;
    cout << "  vol mean : " << meanV;
    cout << "  vol variance : " << vol_sigma << endl;

    //  FT target_v = simu.meanV();

    VectorXd Dvol = mixing * ( target_vol_val  - vol.array()  ) ;

    fill_Delta_DD();

    VectorXd Dw = Delta_solver.solve( Dvol );

    //    copy_weights( T ) ;

    VectorXd w0  = field_to_vctr( sfield_list::w ) ;

    vctr_to_field( w0 + Dw ,  sfield_list::w ) ;

    volumes( T ); // ??

    move_weights( T );

    volumes( T );

    VectorXd vol0( vol );

    vol  = field_to_vctr( sfield_list::vol ) ;

    FT diff = (vol - vol0 ).norm();

    cout << "vol diff: " << diff<< endl;

    totV= vol.sum();
    meanV = totV/ FT( N ); 
    vol_sigma =  ( vol.array() - meanV ).square().sum() / FT( N );

    cout << "  vol mean : " << meanV;
    cout << "  vol variance : " << vol_sigma << endl;

    if( diff < threshold ) break;
  }

  return;
}
