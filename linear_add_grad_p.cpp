//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


void linear::u_add_press_grad( const FT dt ) {

  VectorXd VgradPx, VgradPy;

  DD_times_sfield( sfield_list::p  ,  VgradPx, VgradPy);

  VectorXd vol  = field_to_vctr( sfield_list::vol );
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
  
  vctrs_to_vfield( U_x, U_y , vfield_list::Ustar );

}


// Experimental grads that include the w field:



// Version A

void linear::u_add_press_grad_wdot( const FT dt , const FT beta ) {

  FT ddt = dt;

  VectorXd p = field_to_vctr( sfield_list::p );
  VectorXd w = field_to_vctr( sfield_list::w );
  VectorXd w0= field_to_vctr( sfield_list::w0 );

  VectorXd p_bar = p - ( 0.5 / ( ddt * ddt) ) * ( w - w0 );
  //  VectorXd p_bar = p ; //  - ( 0.5 / ( 4 * ddt * ddt) ) * ( w - w0 );
  //VectorXd p_bar =  0.5*( w - w0 ) / (ddt*ddt);

  VectorXd vol  = field_to_vctr( sfield_list::vol );

  VectorXd gradx = -(DDx.transpose() * p_bar ).array() / vol.array()  ;
  VectorXd grady = -(DDy.transpose() * p_bar ).array() / vol.array()  ;

  // perhaps mean vol would be just fine

  vctrs_add_to_vfield( gradx , grady , vfield_list::gradp );
   
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  // There's a (-1) x (-1) for historical reasons:
  // (-1) in the  definition of grad_ij as -(1/V) D_ij,
  // (-1) in -grad(p) in the Euler equation

  
  U_x = Ustar_x  -  ddt * gradx;
  U_y = Ustar_y  -  ddt * grady;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::Ustar );

}


#ifdef ALWAYS_FALSE

// Version B

void linear::u_add_press_grad_wdot( const FT dt , const FT beta ) {
  FT ddt = dt;

  VectorXd p = field_to_vctr( sfield_list::p );
  VectorXd w = field_to_vctr( sfield_list::w );
  //  VectorXd w0= field_to_vctr( sfield_list::w0 );

  // which??
  VectorXd vol  = field_to_vctr( sfield_list::vol );
  //  VectorXd Dvol  = field_to_vctr( sfield_list::Dvol );
  //  VectorXd vol  = field_to_vctr( sfield_list::Dvol );
  
  VectorXd gradPx = (-DDx.transpose() * p).array() / vol.array()   ;
  VectorXd gradPy = (-DDy.transpose() * p).array() / vol.array()   ;

  // perhaps mean vol would be just fine

  //  VectorXd dwdt =  - ( 0.5 / ( ddt * ddt) ) * ( w - w0 );
  VectorXd dwdt =  ( 0.5 / ( ddt * ddt) ) *  w ;

  
  VectorXd gradwx = (-DDx_fem.transpose() * dwdt).array() / vol.array()   ;
  VectorXd gradwy = (-DDy_fem.transpose() * dwdt).array() / vol.array()   ;

  VectorXd total_gradx =  gradPx.array() + gradwx.array() ;
  VectorXd total_grady =  gradPy.array() + gradwy.array() ;

  vctrs_to_vfield( total_gradx , total_grady , vfield_list::gradp );
   
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  // There's a (-1) x (-1) for historical reasons:
  // (-1) in the  definition of grad_ij as -(1/V) D_ij,
  // (-1) in -grad(p) in the Euler equation

  
  U_x = Ustar_x.array()  -  ddt * total_gradx.array();
  U_y = Ustar_y.array()  -  ddt * total_grady.array();

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}




// Version C

void linear::u_add_press_grad_wdot( const FT dt , const FT beta ) {
  FT ddt = dt;

  VectorXd p = field_to_vctr( sfield_list::p );
  VectorXd w = field_to_vctr( sfield_list::w );
  //  VectorXd w0= field_to_vctr( sfield_list::w0 );

  // which??
  VectorXd vol  = field_to_vctr( sfield_list::vol );
  //  VectorXd Dvol  = field_to_vctr( sfield_list::Dvol );
  //  VectorXd vol  = field_to_vctr( sfield_list::Dvol );
  
  VectorXd gradPx = (-DDx_fem.transpose() * p).array() / vol.array()   ;
  VectorXd gradPy = (-DDy_fem.transpose() * p).array() / vol.array()   ;

  // perhaps mean vol would be just fine

  //  VectorXd dwdt =  - ( 0.5 / ( ddt * ddt) ) * ( w - w0 );
  //  VectorXd dwdt =  ( 0.5 / ( ddt * ddt) ) *  w ;

  
  VectorXd gradwx = (-DDx.transpose() * w).array() / vol.array()   ;
  VectorXd gradwy = (-DDy.transpose() * w).array() / vol.array()   ;

  VectorXd total_gradx =  gradPx.array() + gradwx.array() ;
  VectorXd total_grady =  gradPy.array() + gradwy.array() ;

  vctrs_to_vfield( total_gradx , total_grady , vfield_list::gradp );
   
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  // There's a (-1) x (-1) for historical reasons:
  // (-1) in the  definition of grad_ij as -(1/V) D_ij,
  // (-1) in -grad(p) in the Euler equation

  
  U_x = Ustar_x.array()  -  ddt * total_gradx.array();
  U_y = Ustar_y.array()  -  ddt * total_grady.array();

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}


// Version D

void linear::u_add_press_grad_wdot( const FT dt , const FT beta ) {
  FT ddt = dt;

  VectorXd p = field_to_vctr( sfield_list::p );
  VectorXd w = field_to_vctr( sfield_list::w );
  //  VectorXd w0= field_to_vctr( sfield_list::w0 );

  // which??
  VectorXd vol  = field_to_vctr( sfield_list::vol );
  VectorXd Dvol  = field_to_vctr( sfield_list::Dvol );
  //  VectorXd vol  = field_to_vctr( sfield_list::Dvol );
  
  VectorXd gradPx = (-DDx_fem.transpose() * p).array() / Dvol.array()   ;
  VectorXd gradPy = (-DDy_fem.transpose() * p).array() / Dvol.array()   ;

  // perhaps mean vol would be just fine

  //  VectorXd dwdt =  - ( 0.5 / ( ddt * ddt) ) * ( w - w0 );
  //  VectorXd dwdt =  ( 0.5 / ( ddt * ddt) ) *  w ;

  
  VectorXd gradwx = beta*(DDx.transpose() * w).array() / vol.array()   ;
  VectorXd gradwy = beta*(DDy.transpose() * w).array() / vol.array()   ;

  VectorXd total_gradx =  gradPx + gradwx ;
  VectorXd total_grady =  gradPy + gradwy ;

  vctrs_add_to_vfield( total_gradx , total_grady , vfield_list::gradp );
   
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  // There's a (-1) x (-1) for historical reasons:
  // (-1) in the  definition of grad_ij as -(1/V) D_ij,
  // (-1) in -grad(p) in the Euler equation

  
  U_x = Ustar_x.array()  -  ddt * total_gradx.array();
  U_y = Ustar_y.array()  -  ddt * total_grady.array();

  vctrs_to_vfield( U_x, U_y , vfield_list::Ustar );

}




#endif



void linear::u_add_grads( const FT dt ) {

  VectorXd gradPx,gradPy;
  VectorXd gradsx,gradsy;

  DD_times_sfield( sfield_list::p  ,  gradPx, gradPy);

  MM_times_sfield( sfield_list::s  ,  gradsx, gradsy);

  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  FT ddt = dt;
//  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly
// standard sign for s:
//  VectorXd grad_x = -gradPx + gradsx ;
//  VectorXd grad_y = -gradPy + gradsy ;
// alternative sign for s:
  VectorXd grad_x = -gradPx - gradsx ;
  VectorXd grad_y = -gradPy - gradsy ;

  VectorXd vol  = field_to_vctr( sfield_list::vol );
  // perhaps mean vol would be just fine

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
