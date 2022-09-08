#include"pParticles.h"
#include"simu.h"
#include"data_kept.h"


void create(Triangulation& T, const FT& LL) {

  const FT w0 = 0;

  typedef CGAL::Creator_uniform_2<FT,Point> Creator;

  int Nb = 39; //simu_N_side() ;
//  int Nb = 49; //simu_N_side() ;
  typedef std::vector<Point> vctP;
  vctP points;

  int N=Nb*Nb;

  simu.set_no_of_particles(N);

  points.reserve(N);
  cout << N << " particles to be placed on square lattice" << endl;

  FT spacing = LL/FT(Nb);
  FT side    = LL - spacing;

  points_on_square_grid_2(side/2.0, N, std::back_inserter(points),Creator());;

  //  Vertex_handle v0 = T.insert( wPoint( Point(0,1) ,  w0 ) );


  if(simu.perturb()) {
  
    // perturbation with some velocity field

    for( vctP::iterator vv= points.begin() ;
	 vv != points.end() ;
	 vv++) {
      Point p = *vv;

      FT x=p.x();    FT y=p.y();

      p+= simu.pert_rel()* TG_v( x, y);

      *vv = p;
  
    }
    cout << "each particle perturbed about " << simu.pert_rel()   << endl;
  // // random perturbation
  
  //   CGAL::perturb_points_2(
  // 			   points.begin(), points.end(),
  // 			   simu.pert_rel()* spacing );
    //    cout << "each particle perturbed about " << simu.pert_rel()* spacing  << endl;
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

  backup( T );

  return;

}

void expand(Triangulation& T, const FT& LL) {

  std::vector<Vector_2> dirs(8);
  dirs[0] = Vector_2( 1 ,  0 );
  dirs[1] = Vector_2( 1 ,  1 );
  dirs[2] = Vector_2( 0 ,  1 );
  dirs[3] = Vector_2(-1 ,  1 );
  dirs[4] = Vector_2(-1 ,  0 );
  dirs[5] = Vector_2(-1 , -1);
  dirs[6] = Vector_2( 0 , -1 );
  dirs[7] = Vector_2( 1 , -1 );
  
  cout << "Copying inner square on surrounding " << endl;

  vector<data_kept> prev;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) {
    data_kept data(fv);

    data.idx = -1;
    Point r0 = fv->point().point(); 

    for( int i=0 ; i<8 ; i++ ) {
      Vector_2 disp = LL * dirs[i];
      Point rr = r0 +  disp; 

      data.pos = rr;

      prev.push_back (data);

    }

  }

  for(vector<data_kept>::iterator data=prev.begin();
      data!=prev.end();
      data++) {

    //    cout << "Inserting back at " << data->pos << endl ;
    
    Vertex_handle fv=T.insert( wPoint( data->pos , data->w )  );

    data->restore(fv);
    
  }

}
