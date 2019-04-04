#include"pParticles.h"
#include"simu.h"


void create(Triangulation& T, const FT& LL) {

  const FT w0 = 0;

  typedef CGAL::Creator_uniform_2<FT,Point> Creator;

  int Nb = 20; //simu_N_side() ;
  typedef std::vector<Point> vctP;
  vctP points;

  int N=Nb*Nb;

  simu.set_no_of_particles(N);

  points.reserve(N);
  cout << N << " particles placed on square lattice" << endl;

  FT spacing=LL/FT(Nb+0);
  FT side=LL-2*spacing;

  points_on_square_grid_2(side/2.0, N, std::back_inserter(points),Creator());;


  //  Vertex_handle v0 = T.insert( wPoint( Point(0,1) ,  w0 ) );

  if(simu.perturb()) {
    CGAL::perturb_points_2(
			   points.begin(), points.end(),
			   simu.pert_rel()* spacing );
    cout << "each particle perturbed about " << simu.pert_rel()* spacing  << endl;
  }

  cout << "Inserting interior points" << endl;
    
  for( vctP::iterator vv= points.begin() ;
       vv != points.end() ;
       vv++) {
    Point p = *vv;
    wPoint wp = wPoint( p , w0 ) ;
    // 
    //    T.insert( wp );

    Vertex_handle vtx = T.push_back( wp );

  }

  cout << "Inserting frame " << endl;

  FT edge = LL /2 ;
  Point p;
  wPoint wp;
  Vertex_handle vtx;
  
  for(int i=0 ; i < Nb ; i++) {
    FT var = - side/2 + side* i/FT(Nb-1) ;

    p = Point( var , edge );
    vtx = T.push_back( wPoint( p , w0 ) );
    vtx->idx.set( -1 );

    p = Point( var , -edge );
    vtx = T.push_back( wPoint( p , w0 ) );
    vtx->idx.set( -1 );
 
    //    if(i>0) {
    p = Point(  edge , var );
    vtx = T.push_back( wPoint( p , w0 ) );
    vtx->idx.set( -1 );

    p = Point( -edge , var );
    vtx = T.push_back( wPoint( p , w0 ) );
    vtx->idx.set( -1 );

      //    }
  }

  // corners

  int signx = 1;
  int signy = 1;
  
  for(int i = 0 ; i < 2 ; i++) {
    signx *= -1 ;
    for(int j = 0 ; j < 2 ; j++) {
        signy *= -1 ;
	
  	p = Point( signx*edge , signy*edge );
  	wp = wPoint( p , w0 ) ;
  	vtx = T.push_back( wp );
  	vtx->idx.set( -1 );
    }
  }

  
  // for(F_v_it vit=T.finite_vertices_begin();
  //     vit != T.finite_vertices_end();
  //     vit++) {
  //   Point p=vit->point().point();
  // }
  
  backup( T );

  return;

}


