//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


void linear::u_star( void ) {
  VectorXd usx, usy;
  
  vfield_to_vctrs( vfield_list::U , usx, usy );

  vctrs_to_vfield( usx, usy, vfield_list::Ustar );

  return;
}


