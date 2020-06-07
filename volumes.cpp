#include"pParticles.h"
#include"simu.h"


// #define WARNING

typedef vector<Point> vvP;


// Compute Voronoi volumes (i.e. areas, in 2D)

void volumes(Triangulation& T) {

  const FT threshold=1e-10;
  const FT threshold2= threshold*threshold;
  FT totalV = 0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {
    fv->vol.reset();
    fv->I.reset();
    fv->centroid.reset( );
  }

#ifdef FEM

  // Volumes, FEM shape functions (Delaunay areas)
  
  for(F_f_it ff=T.finite_faces_begin();
      ff!=T.finite_faces_end();
      ff++)    {

    Vertex_handle v0 = ff->vertex(0);
    Vertex_handle v1 = ff->vertex(1);
    Vertex_handle v2 = ff->vertex(2);

    Point p0 = v0->point().point();
    Point p1 = v1->point().point();
    Point p2 = v2->point().point();

    Triangle tr( p0 ,  p1 , p2 );

    FT area = std::fabs( tr.area() );

    v0->vol += area / 3.0;
    v1->vol += area / 3.0;
    v2->vol += area / 3.0;
  }

#else

  // Volumes, Voronoi cells

  for(F_e_it fe=T.finite_edges_begin();
      fe!=T.finite_edges_end();
      fe++)    {

    Face_handle f =  fe -> first ;
    int i0 = fe -> second;

    Vertex_handle vi = f->vertex( (i0+1) % 3);
    Point pi = vi->point().point();

    Vertex_handle vj = f->vertex( (i0+2) % 3);
    Point pj = vj->point().point();
    
    CGAL::Object o = T.dual(fe);

    const Segment * Vor_segment = CGAL::object_cast<Segment>( &o );

    if (! Vor_segment ) continue;

    Point p1 = Vor_segment->source() ;
    Point p2 = Vor_segment->target() ;
      
    Triangle tri( pi ,  p1 , p2 );
    Triangle trj( pj ,  p1 , p2 );

    FT ar_i = std::fabs( tri.area() );
    FT ar_j = std::fabs( trj.area() );

    vi->vol +=  ar_i ;
    vj->vol +=  ar_j ;

    totalV += ar_i;
    totalV += ar_j;
    {
      Vector_2 vi1 = p1 - pi ;
      Vector_2 vi2 = p2 - pi ;

      vi->I += ar_i / 6 * 
	( vi1.squared_length() +
	  vi2.squared_length() +
	  vi1 * vi2 );
    }

    {
      Vector_2 vj1 = p1 - pj ;
      Vector_2 vj2 = p2 - pj ;
    
      vj->I += ar_j / 6 * 
	( vj1.squared_length() +
	  vj2.squared_length() +
	  vj1 * vj2 );
    }

    // CGAL::ORIGIN needed because points cannot be added, or multiplied
    Vector_2 tri_ctr_v = CGAL::centroid( tri ) - CGAL::ORIGIN;
    Vector_2 trj_ctr_v = CGAL::centroid( trj ) - CGAL::ORIGIN;

    vi->centroid =  vi->centroid.val() + ar_i * tri_ctr_v;
    vj->centroid =  vj->centroid.val() + ar_j * trj_ctr_v;

  }

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    FT a = fv->vol.val();

    Vector_2 ctr_v = fv->centroid.val() - CGAL::ORIGIN;

    fv->centroid = CGAL::ORIGIN + ctr_v/a;

    Point cc = fv->centroid.val();
    Point p  = fv->point().point();

    Vector_2 dA =  a *( cc - p );

    fv->dd.set( dA );

    fv->dd2.set(  dA.squared_length()  );

  }

#endif

  
  //  cout << "Volumes: total = " << totalV << " ; ";

  FT  inner_V = 0;
  int inner   = 0;
  
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)
    if (fv->idx() > -1 ) {
      inner_V +=  fv->vol.val();
      ++inner;
    }

  //  cout << "inner = " << inner_V << " ; ";

  simu.set_totalV( totalV );
  simu.set_innerV( inner_V );

  //  cout << "mean inner V = " << simu.meanV() << endl;

  
}
