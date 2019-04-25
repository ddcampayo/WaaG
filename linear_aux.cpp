//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


void linear::u_star( void ) {
  VectorXd usx, usy;
  
  vfield_to_vctrs( vfield_list::U0 , usx, usy );

  vctrs_to_vfield( usx, usy, vfield_list::Ustar );

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


