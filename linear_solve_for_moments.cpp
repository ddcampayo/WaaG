//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Iterative process to adjust moments so that
// the weighted Voronoi diagram has constant moment

void linear::solve_for_moments( ) {

  cout << "Equalizing moments " << endl;
  
  volumes( T );
  copy_weights( T );
  VectorXd I0 = field_to_vctr( sfield_list::I0 ) ;

  //  VectorXd target_vol( vol ) ;
  //  target_vol.setConstant( target_vol_val );
  
  const int max_iter = 100;
  const FT threshold = 1e-16;
  const FT mixing = 1;
  
  for( int iter=0 ; iter< max_iter ; iter++) {

    VectorXd I  = field_to_vctr( sfield_list::I ) ;

    VectorXd DI = mixing * ( I0  - I ) ;

    fill_Delta_DD();

    VectorXd Dw = GG_solver.solve( DI );

    //    copy_weights( T ) ;

    VectorXd w0  = field_to_vctr( sfield_list::w ) ;

    vctr_to_field( w0 + Dw ,  sfield_list::w ) ;

    volumes( T ); // ??

    move_weights( T );

    volumes( T );

    VectorXd I0( I );

    I  = field_to_vctr( sfield_list::I ) ;

    FT diff = (I - I0 ).norm();

    cout << "I diff: " << diff<< endl;

    cout << "  solving for  I, iter : " << iter << endl;

    if( diff < threshold ) break;
  }

  cout << "Equalized moments " << endl;

  return;
}
