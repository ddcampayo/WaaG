//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure
// The famous pressure Poisson equation

void linear::p_equation(const FT dt ) {

  cout << "Solving pressure equation " << endl;
  
  fill_Delta_DD(); // This may be important -- or not

  VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  // A
  // Approximate Laplacian ~ Delta 
  // VectorXd p =  Delta_solver.solve( divUstar );
  // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  // vctr_to_field( -0.5 * p / ddt ,  sfield_list::p ) ;

  // B
  //  Laplacian as div of grad :
  // VectorXd p =  LL_solver.solve( divUstar );
  // vctr_to_field( p / ddt ,  sfield_list::p ) ;

  // C
  // As B, but Dvol source

  volumes( T );

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  FT totV= vol.sum();
  int N = vol.size();
  FT meanV = totV/ FT( N );
 
  FT target_vol_val = meanV;

  for(int iter=0 ; iter < 4 ; iter++) {
    VectorXd Dvol = vol.array() - target_vol_val  ;
    VectorXd p0 =  field_to_vctr( sfield_list::p ) ;
    VectorXd p  =  LL_solver.solve( Dvol );

    vctr_to_field( p0 + p / ( ddt * ddt) ,  sfield_list::p ) ;

    cout << "Pressure iter " << iter
	 << " rel Dvol" << Dvol.sum()/Dvol.size()
	 << endl;
  }
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

  U_x = Ustar_x.array() - ddt * gradPx.array() / vol.array();
  U_y = Ustar_y.array() - ddt * gradPy.array() / vol.array();
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}
