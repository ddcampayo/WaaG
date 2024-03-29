//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


void linear::test_operators( void ) {

  // this makes u = grad p, for debugging purposes

  fill_Delta_DD();

  VectorXd gradPx,gradPy;

  DD_times_sfield( sfield_list::p  ,  gradPx, gradPy);

  VectorXd vol  = field_to_vctr( sfield_list::vol );

  VectorXd U_x, U_y;

  U_x = gradPx.array() / vol.array()  ;
  U_y = gradPy.array() / vol.array() ;

  vctrs_to_vfield( U_x, U_y , vfield_list::U );
  
}



void linear::u_star( void ) {
  VectorXd usx, usy;
  
  //  vfield_to_vctrs( vfield_list::U0 , usx, usy );

  //  vctrs_to_vfield( usx, usy, vfield_list::Ustar );


   copy( vfield_list::U0 ,  vfield_list::Ustar );
  
  return;
}


void linear::reset_p( void ) {

  VectorXd p = field_to_vctr( sfield_list::p );

  vctr_to_field( 0*p , sfield_list::p );

  return;
}


void linear::reset_s( void ) {

  VectorXd s = field_to_vctr( sfield_list::s );

  vctr_to_field( 0*s , sfield_list::s );

  return;
}


void linear::copy(const sfield_list::take from, sfield_list::take to  ) {

  VectorXd p = field_to_vctr( from );

  vctr_to_field( p , to );

  return;
}


void linear::copy(const FT a , const sfield_list::take from, sfield_list::take to  ) {

  VectorXd p = a*field_to_vctr( from );

  vctr_to_field( p , to );

  return;
}

void linear::copy(const vfield_list::take from, vfield_list::take to  ) {

  VectorXd px, py;
  vfield_to_vctrs( from , px , py );

  vctrs_to_vfield( px , py , to );

  return;
}


