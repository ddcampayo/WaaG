// Interpolate values of fields at point p0 from values at
// vertices of triangle where it is at

#include"pParticles.h"
//#include"simu.h"


// z component of vector product v1 x v2
FT vect_prod(const Vector_2& v1 , const Vector_2& v2 ) {
  return v1.x() * v2.y() - v1.y() * v2.x() ;
}


void FEM_hs(const Face_handle& fc, const Point& pp,
	    std::vector<Vertex_handle>& v,
	    std::vector<FT>& hh    ) {  

  for(int i=0; i < 3 ; ++i)    v[i] = fc->vertex(i);

  std::vector<Point> p(3);

  for(int i=0; i < 3 ; ++i)    p[i] = v[i]->point().point();

  Vector_2 v01 = Vector_2( p[0] , p[1] );
  Vector_2 v02 = Vector_2( p[0] , p[2] );
  Vector_2 rr  = Vector_2( p[0] , pp   );

  FT twice_A = vect_prod( v01 , v02 ) ;  // 2x area triangle

  hh[1] = vect_prod( rr , v02 ) / twice_A;

  hh[2] = vect_prod( v01 , rr ) / twice_A;

  hh[0] = 1 - hh[1] - hh[2];

  return;
}



Vector_2 values_at_v(const Triangulation& T , const Point& p0, const vfield_list::take v_field) {
  
  std::vector<Vertex_handle> v(3);
  std::vector<FT> hh(3);

  //  cout << "Locating p0 = " << p0 ;

  // caution: locates _weighted_ points. 

  wPoint wp0( p0 , 0 );

  Face_handle fc = T.locate( wp0 );
  //  cout << "   located" << endl;

  FEM_hs( fc , p0, v, hh);

  Vector_2 vv( 0 , 0 );  // accumulator

  for( int j = 0; j < 3 ; ++j )
    vv = vv + hh[ j ] * v[ j ]->vfield( v_field ).val() ;

  return vv;
}
