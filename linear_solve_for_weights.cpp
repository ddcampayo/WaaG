//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Iterative process to adjust weights so that
// the weighted Voronoi diagram has constant volumes

void linear::solve_for_weights( const FT dt ) {

  //  cout << "Equalizing volumes " << endl;
  
  volumes( T );
  copy_weights( T );
  VectorXd vol0 = field_to_vctr( sfield_list::vol0 ) ;

  //  VectorXd target_vol( vol ) ;
  //  target_vol.setConstant( target_vol_val );
  
  const int max_iter = 100;
  const FT threshold = 1e-5;
  const FT mixing = 1; // 1: only new iter; 0: only old

  FT diff , vol_sigma  , meanV ;

  int iter=0;

  for( ; iter< max_iter ; iter++) {

    VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

    FT totV= vol.sum();

    int N = vol.size();

    meanV = totV/ FT( N );
  
    vol_sigma =  ( vol.array() - meanV ).square().sum() / FT( N );

    FT diff = (vol - vol0 ).norm() / vol.norm(); // relative!

    cout << "Volumes rel diff: " << diff<< endl;

    cout << "   solving for  weights, iter : " << iter << endl;

    //    FT target_vol_val = meanV;

    // cout << "It:  " << iter;
    // cout << "  vol mean : " << meanV;
    // cout << "  vol variance : " << vol_sigma << endl;

    //  FT target_v = simu.meanV();

    //    VectorXd Dvol = mixing * ( target_vol_val  - vol.array()  ) ;
    VectorXd Dvol = mixing * ( vol0  - vol ) ;

    fill_Delta_DD( dt );

    VectorXd Dw = Delta_solver.solve( Dvol );

    //    copy_weights( T ) ;

    VectorXd w0  = field_to_vctr( sfield_list::w ) ;

    vctr_to_field( w0 + Dw ,  sfield_list::w ) ;

    volumes( T ); // ??

    move_weights( T );

    volumes( T );

    //    VectorXd vol0( vol );

    //    vol  = field_to_vctr( sfield_list::vol ) ;

    //    diff = (vol - vol0 ).norm();

    //    cout << "vol diff: " << diff<< endl;

    totV= vol.sum();
    meanV = totV/ FT( N ); 
    vol_sigma =  ( vol.array() - meanV ).square().sum() / FT( N );

    //    cout << "  vol mean : " << meanV;
    //    cout << "  vol variance : " << vol_sigma << endl;

    if( diff < threshold ) break;
  }

  cout << "Volumes equalized after " << iter << " iterations " << endl;

  cout << "vol diff: " << diff<< "   ";
  cout << "  vol mean : " << meanV << "    ";
  cout << "  vol variance : " << vol_sigma << endl;
 
  return;
}
