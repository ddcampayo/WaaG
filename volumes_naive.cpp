#include"pParticles.h"
#include"simu.h"

// Compute Voronoi volumes (i.e. areas, in 2D)
//
//
void volumes(Triangulation& T) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {
    fv->vol.reset();
    //    fv->vol.set(0);
  }

  F_e_it eit = T.finite_edges_begin();

  FT totalV=0;

  for ( ; eit !=T.finite_edges_end(); ++eit) {

    CGAL::Object o = T.dual(eit);

    const Segment * Vor_segment = CGAL::object_cast<Segment>( &o );

    if (! Vor_segment ) continue;

    //    cout << " l0 = " << std::sqrt(Vor_segment->squared_length()) << endl;
    
    Point p1 = Vor_segment->source() ;
    Point p2 = Vor_segment->target() ;
    
    Face_handle f =  eit -> first ;
    int i = eit -> second;

    Vertex_handle v12 = f->vertex( (i+1) % 3);
    Vertex_handle v21 = f->vertex( (i+2) % 3);

    Point p12 = v12->point().point() ;
    Point p21 = v21->point().point() ;

    Triangle tr_1_2_12( p1 , p2 , p12);

    FT vol12 = std::fabs(tr_1_2_12.area());

    v12->vol += vol12;

    totalV += vol12;

    //    cout << " vol0 = " << vol12 << endl;

//    continue;

    Triangle tr_1_2_21( p1 , p2 , p21); // Notice ordering, area may be negative otherwise

    FT vol21 = std::fabs(tr_1_2_21.area());

    v21->vol += vol21;

    totalV += vol21;

  }

  cout << "Volumes computed; total V = " << totalV << " ; ";

  totalV = 0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)
    if (fv->idx() > 0 )
      totalV +=  fv->vol.val();

  cout << "inner V = " << totalV << " ; ";

  simu.set_totalV( totalV );

  cout << "mean V = " << simu.meanV() << endl;
  
  
}
