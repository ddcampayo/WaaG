//#define PRESSURE_PPE_DIV_SOURCE
// TODO: better define different PPE functions, independent

//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"

// Solve for pressure
// The famous pressure Poisson equation

void linear::p_equation(const FT dt , const bool ws ) {

  // choose!!

  //   p_equation_divgrad_div_source(dt,ws);
  //   p_equation_lapl_div_source(dt);
  p_equation_lapl_Dvol_source( dt );

  return;
}

//  Laplacian = - Delta / (2 V), div v as source

void linear::p_equation_lapl_div_source(const FT dt ){

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  
  VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );
  VectorXd p =  Delta_solver.solve( divUstar );
  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  vctr_to_field( -0.5 * p / ddt ,  sfield_list::p ) ;

  return;
}


// same as above, but div u calculated as in FEM

void linear::p_equation_lapl_div_source_fem(const FT dt ){

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  
  VectorXd divUstar  =  DD_scalar_vfield_fem( vfield_list::Ustar );
  VectorXd p =  Delta_solver.solve( divUstar );
  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  vctr_to_field( -0.5 * p / ddt ,  sfield_list::p ) ;

  return;
}




//  Laplacian = - Delta / (2 V), Delta vol as source

void linear::p_equation_lapl_Dvol_source(const FT dt ){

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  VectorXd vol0  = field_to_vctr( sfield_list::vol0 ) ;

  VectorXd Dvol = vol.array() - vol0.array()  ;

  VectorXd Dp =  Delta_solver.solve( Dvol );

  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  vctr_to_field( -0.5 * Dp / ( ddt * ddt) , sfield_list::p  ) ;

  
  return;
}


// //  Laplacian = - Delta / (2 V), Delta FEM vol as source

// void linear::p_equation_lapl_Dvol_source_fem(const FT dt ){

//   cout << "Solving pressure equation " << endl;
  
//   //  fill_Delta_DD(); // This may be important -- or not

//   FT ddt = dt;
//   if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

//   VectorXd vol  = field_to_vctr( sfield_list::Dvol ) ;

//   VectorXd vol0  = field_to_vctr( sfield_list::Dvol0 ) ;

//   VectorXd Dvol = vol.array() - vol0.array()  ;

//   VectorXd Dp =  Delta_solver.solve( Dvol );

//   // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
//   vctr_to_field( -0.5 * Dp / ( ddt * ddt) , sfield_list::p  ) ;

  
//   return;
// }



//  Laplacian = div of grad, div v as source

void linear::p_equation_divgrad_div_source(const FT dt , const bool ws ) {

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );
  VectorXd p;

  // Possible correction due to w  field .-
  if( ws ) {
    VectorXd w  = field_to_vctr( sfield_list::w );
    VectorXd w0 = field_to_vctr( sfield_list::w0 );

    VectorXd Delta_w = Delta * ( w - w0);
    p =  LL_solver.solve( divUstar + Delta_w / ddt );
  }
  else
    p =  LL_solver.solve( divUstar  );
  
  vctr_to_field( p / ddt ,  sfield_list::p ) ;

  return;
}


//  Laplacian = div of grad, D vol as source

void linear::p_equation_divgrad_Dvol_source(const FT dt , const bool ws ) {

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  
  // diagnostics on volumes.-

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  VectorXd vol0  = field_to_vctr( sfield_list::vol0 ) ;

  VectorXd Dvol = vol.array() - vol0.array()  ;

  int N = vol.size();

  FT Dvol_sigma =  Dvol.array().square().sum() / N ; // / FT( vol.size() );
  FT Dvol_mean  =  vol.array().sum() / N ; // / FT( vol.size() );

  cout << "Pressure  "
       << " rel Dvol std dev: " << sqrt( Dvol_sigma ) / Dvol_mean 
       << endl;

  // C1: LL Laplacian
  VectorXd Dp  =  LL_solver.solve( Dvol );

  vctr_to_field( Dp / ( ddt * ddt) , sfield_list::p  ) ;


  return;
}



void linear::p_equation_s(const FT dt ) {

  cout << "Solving pressure equation " << endl;
  
  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );
  VectorXd p;

  VectorXd s  = field_to_vctr( sfield_list::s );

  VectorXd LNs = LN * s;

  p =  LL_solver.solve( divUstar / dt + LNs );
  
  vctr_to_field( p ,  sfield_list::p ) ;

  return;
}


void linear::u_add_press_grad( const FT dt ) {

  VectorXd VgradPx, VgradPy;

  DD_times_sfield( sfield_list::p  ,  VgradPx, VgradPy);

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine

  VectorXd gradPx = gradPx.array() / vol.array()  ;
  VectorXd gradPy = gradPy.array() / vol.array()  ;

  vctrs_to_vfield( gradPx , gradPy , vfield_list::gradp );
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly


  // There's a (-1) x (-1) for historical reasons:
  // (-1) in the  definition of grad_ij as -(1/V) D_ij,
  // (-1) in -grad(p) in the Euler equation

  U_x = Ustar_x - ddt * gradPx;
  U_y = Ustar_y - ddt * gradPy;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}


void linear::u_add_press_grad_fem( const FT dt ) {

  VectorXd VgradPx,VgradPy;

  DD_times_sfield_fem( sfield_list::p  ,  VgradPx, VgradPy);

  VectorXd vol  = field_to_vctr( sfield_list::Dvol );
  // perhaps mean vol would be just fine

  VectorXd gradPx = VgradPx.array() / vol.array()  ;
  VectorXd gradPy = VgradPy.array() / vol.array()  ;

  vctrs_to_vfield( gradPx , gradPy , vfield_list::gradp );
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly


  // There's a (-1) x (-1) for historical reasons:
  // (-1) in the  definition of grad_ij as -(1/V) D_ij,
  // (-1) in -grad(p) in the Euler equation
  U_x = Ustar_x - ddt * gradPx;
  U_y = Ustar_y - ddt * gradPy;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}



#ifdef ALWAYS_FALSE

// Version A
void linear::u_add_press_grad_wdot(  const FT dt ) {

  FT ddt = dt;

  VectorXd p = field_to_vctr( sfield_list::p );
  VectorXd w = field_to_vctr( sfield_list::w );
  VectorXd w0= field_to_vctr( sfield_list::w0 );

  VectorXd p_bar = p - ( 0.5 / ( ddt * ddt) ) * ( w - w0 );
  //  VectorXd p_bar = p ; //  - ( 0.5 / ( 4 * ddt * ddt) ) * ( w - w0 );
  //VectorXd p_bar =  0.5*( w - w0 ) / (ddt*ddt);
  
  VectorXd gradx = -DDx.transpose() * p_bar;
  VectorXd grady = -DDy.transpose() * p_bar;

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  // There's a (-1) x (-1) for historical reasons:
  // (-1) in the  definition of grad_ij as -(1/V) D_ij,
  // (-1) in -grad(p) in the Euler equation

  
  U_x = Ustar_x.array()  -  ddt * gradx.array() / vol.array()  ;
  U_y = Ustar_y.array()  -  ddt * grady.array() / vol.array()  ;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}

#endif

// Version B

void linear::u_add_press_grad_wdot(  const FT dt ) {

  FT ddt = dt;

  VectorXd p = field_to_vctr( sfield_list::p );
  VectorXd w = field_to_vctr( sfield_list::w );
  VectorXd w0= field_to_vctr( sfield_list::w0 );

  VectorXd gradx = -DDx.transpose() * p;
  VectorXd grady = -DDy.transpose() * p;

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine

  //  VectorXd dwdt =  - ( 0.5 / ( ddt * ddt) ) * ( w - w0 );
  VectorXd dwdt =  ( 0.5 / ( ddt * ddt) ) *  w ;
  
  VectorXd gradwx = -DDx_fem.transpose() * dwdt;
  VectorXd gradwy = -DDy_fem.transpose() * dwdt;

  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  // There's a (-1) x (-1) for historical reasons:
  // (-1) in the  definition of grad_ij as -(1/V) D_ij,
  // (-1) in -grad(p) in the Euler equation

  VectorXd total_gradx = gradx.array() + gradwx.array() ;
  VectorXd total_grady = grady.array() + gradwy.array() ;
  
  U_x = Ustar_x.array()  -  ddt * total_gradx.array()  / vol.array()  ;
  U_y = Ustar_y.array()  -  ddt * total_grady.array()  / vol.array()  ;

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
  VectorXd grad_x = - gradPx + gradsx ;
  VectorXd grad_y = - gradPy + gradsy ;

  U_x = Ustar_x.array() + ddt * grad_x.array() / vol.array()  ;
  U_y = Ustar_y.array() + ddt * grad_y.array() / vol.array() ;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}




void linear::om_add_press_grad( const FT dt ) {

  //  cout << "Omega equation \n";
  
  VectorXd gradPx,gradPy;

  DD_times_sfield( sfield_list::p  ,  gradPx, gradPy);

  VectorXd r_x, r_y;

  vfield_to_vctrs(  vfield_list::dd , r_x , r_y );

  VectorXd I  = field_to_vctr( sfield_list::I );

  VectorXd tau =
      r_x.array() * gradPy.array()
    - r_y.array() * gradPx.array() ; // torque of V grad p

  VectorXd om0  = field_to_vctr( sfield_list::om0 );

  FT ddt = dt;
  
  VectorXd om = om0.array() - ddt * tau.array() / I.array()  ;

  vctr_to_field( om ,  sfield_list::om ) ;

}



void linear::u_add_angular( void ) {

  //  cout << "adding angular velocity \n";

  VectorXd U_x, U_y;

  vfield_to_vctrs(  vfield_list::U , U_x, U_y );

  VectorXd r_x, r_y;

  vfield_to_vctrs(  vfield_list::dd , r_x , r_y );

  VectorXd om  = field_to_vctr( sfield_list::om );

  U_x = U_x.array() - om.array() * r_y.array() ;
  U_y = U_y.array() + om.array() * r_x.array() ;

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}
