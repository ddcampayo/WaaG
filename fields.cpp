#include"pParticles.h"
#include"simu.h"

const double pi = M_PI;

Vector_2 field_rotation(const FT x,const FT y ) {

  const FT cut   = 0.4;
  const FT width = 0.1;
  
  FT r2= x*x + y*y;

  FT r = std::sqrt(r2);

  if(  r > cut  )
    return CGAL::NULL_VECTOR ;
  else {
    FT cross = 1; //(1 - std::tanh( ( r - cut) / width  )) / 2;

    return cross * 2 * M_PI * Vector_2( -y , x ) ;

  }

}


void set_vels_rotating(Triangulation& T) {

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    vit->U.set( field_rotation(x,y) );

  }

  return;
}


void set_vels_Lamb_Oseen(Triangulation& T) {

  const FT tiny = 1e-10;
  const FT Gamma_over_2pi = 1.0 / (2 * pi );
  const FT r_c = 0.1;
  const FT r2_c = r_c * r_c;
  
  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    FT r2 =  x*x + y*y ;
    FT r  =  std::sqrt( r2 ) ;
    Vector_2  u_theta( Vector_2( -y , x ) / ( r + tiny) ) ;
    FT amp = Gamma_over_2pi * ( 1 -  std::exp( - r2 / r2_c ) )  / ( r + tiny ) ;

    vit->U.set( amp * u_theta );
//    cout << r << "  " << amp << endl;
  }

  return;
}


Vector_2 Gresho_v( const FT x, const FT y) {

  const FT tiny = 1e-10;

  const FT rc1 = 0.2;
  const FT rc2 = 0.4;

  FT r2 =  x*x + y*y ;
  FT r  =  std::sqrt( r2 ) ;
  Vector_2  u_theta( Vector_2( -y , x ) / ( r + tiny) ) ;

  FT amp = 0;

  if( r <= rc1 ) amp = 5*r ;
  else if( r <= rc2 )  amp = 2 - 5*r;

  return amp * u_theta ;

}


void set_vels_Gresho(Triangulation& T) {
  
  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    vit->U.set( Gresho_v( x , y) ) ;

  }

  return;
}


FT L2_vel_Gresho( Triangulation& T) {

  FT L2=0;
  int nn=0;
  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    Vector_2 U0 = Gresho_v( x , y) ;
    Vector_2 U  = vit->U.val();
    L2 += ( U - U0 ).squared_length();
    ++nn;
    
  }

  return L2 / nn;
}
