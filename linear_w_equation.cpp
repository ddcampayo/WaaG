//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for weights.
// Similar to solve_for_weights, but not iterative

void linear::w_equation( ) {

  cout << "Solving weight equation " << endl;
  
  volumes( T );

  //  copy_weights( T ) ;

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  FT totV= vol.sum();

  int N = vol.size();

  FT meanV = totV/ FT( N );
  
  FT vol_sigma =  ( vol.array() - meanV ).square().sum() / FT( N );

  FT target_vol_val = meanV;

  cout << "  vol mean : " << meanV;
  cout << "  vol variance : " << vol_sigma << endl;

    //  FT target_v = simu.meanV();
  fill_Delta_DD();
  
  // div of displacement:
  //  VectorXd divDr  =  DD_scalar_vfield( vfield_list::Dr );
  //  VectorXd Dw = Delta_solver.solve( divDr );

//  VectorXd vol00  = field_to_vctr( sfield_list::vol0 ) ;
//  VectorXd Dw = Delta_solver.solve( vol - vol00 );

  // given by volume departure
  VectorXd Dvol  =  target_vol_val - vol.array() ;
  VectorXd Dw = Delta_solver.solve( Dvol );

  //  cout << "Dw " << Dw << endl;
  
    //    copy_weights( T ) ;

  VectorXd w0  = field_to_vctr( sfield_list::w0 ) ;
  
  vctr_to_field( w0 + Dw ,  sfield_list::w ) ;

//    volumes( T ); // ??

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

  return;
}
