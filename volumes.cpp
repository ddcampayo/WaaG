#include"pParticles.h"
#include"simu.h"


// #define WARNING

typedef vector<Point> vvP;


// second moment of area wrt p -- Ix + Iy
// https://en.wikipedia.org/wiki/Second_moment_of_area#Any_polygon
FT moi( const Point& p , const vvP& vs ) {

  FT mm = 0;
  int N = vs.size();
  for(int i=0 ; i < N ; i++ ) {
    Vector_2 pi = vs[  i         ] - p ;
    Vector_2 pj = vs[ (i+1 ) % N ] - p ;

    FT xi= pi.x();
    FT yi= pi.y();

    FT xj= pj.x();
    FT yj= pj.y();

    // ugly.-
    mm += std::fabs(
		    (  xi * xi + yi * yi +
		       xj * xj + yj * yj +
		       xi * xj + yi * yj  ) *
		    ( xi * yj - xj * yi ) / 12  );

    // elegant.-
    //    FT area =  ( xi * yj - xj * yi ) / 2 ;    
    // mm += std::fabs(
    // 		    ( pi.squared_length() +
    // 		      pj.squared_length() +
    // 		      pi * pj ) * area / 6  );
  }

  return mm;
}

// Compute Voronoi volumes (i.e. areas, in 2D)
// This is a combination of "naive" (for areas)
// and "abstruse" (for moments of area)

void volumes(Triangulation& T) {

  const FT threshold=1e-10;
  const FT threshold2= threshold*threshold;
  FT totalV = 0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  
    fv->vol.reset();


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

    FT ar_i = tri.area();
    FT ar_j = trj.area();

    vi->vol += std::fabs( ar_i );
    vj->vol += std::fabs( ar_j );

    totalV += ar_i;
    totalV += ar_j;

  }
     

  // Other: polygons, barycenters ...
  
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {
    //    fv->vol.reset();

    Edge_circulator edge  = T.incident_edges( fv );

    Edge_circulator first = edge;

    first--; // avoid last one entirely
 
    vvP poly_vertices;

    int nn=1;

    // Collecting of Voronoi vertices --- TODO: clearly improvable

    Point p0 = fv->point().point();

    do {

      if(T.is_infinite( edge ) ) {
	++edge;
	continue;
      }
	
      CGAL::Object o = T.dual(edge);

      const Segment * Vor_segment = CGAL::object_cast<Segment>( &o );

      if (! Vor_segment ) {
	++edge;
	continue;
      }

      Point p1 = Vor_segment->source() ;
      Point p2 = Vor_segment->target() ;
      
      //Triangle tr( p0 ,  p1 , p2 );
      //FT tr_area = std::fabs( tr.area() );
      //area += tr_area;

      //      totalV += tr_area;

      if( nn == 1 ) {
	poly_vertices.push_back( p1 ) ;
	poly_vertices.push_back( p2 ) ;
      }
      else if( nn == 2) {
	Point p_new;
	if ( p1 == poly_vertices[1] )  p_new = p2;
	else if ( p2 == poly_vertices[1] ) p_new = p1;
	else if ( p1 == poly_vertices[0] ) {
	  poly_vertices[0] = poly_vertices[1];
	  poly_vertices[1] = p1;
	  p_new = p2;
	}
	else if ( p2 == poly_vertices[0] ) {
	  poly_vertices[0] = poly_vertices[1];
	  poly_vertices[1] = p2;
	  p_new = p1;
	}
	else{
	  // cout << "WARNING: conflict in polygon building, at n="
	  // 	   << nn <<endl;
	  // cout << fv->point().point() << endl ;
	  // cout << p1 << endl ;
	  // cout << p2 << endl ;

	  // cout  << endl ;
	  
	  // for( auto pp : poly_vertices )
	  //   cout << pp << endl; 

	  // cout  << endl ;
	  // //	  cout << poly_vertices << endl ;
	  
	  ++edge;
	  continue;
	} 
	poly_vertices.push_back( p_new ) ;
      }
      else if ( p2 == poly_vertices[nn-1] ) 
	poly_vertices.push_back( p1 ) ;
      else if ( p1 == poly_vertices[nn-1] )
	poly_vertices.push_back( p2 ) ;
      else{
	// cout << "WARNING: conflict in polygon building, at n="
	// 	   << nn <<endl;
	//   cout << fv->point().point() << endl ;
	//   cout << p1 << endl ;
	//   cout << p2 << endl ;
	//   //	  cout << poly_vertices << endl ;

	  ++edge;
	  continue;
      }

      // else {
      // 	if( p1 == poly_vertices[nn-1] ) 
      // 	  poly_vertices.push_back( p2 ) ;
      // 	else   if( p2 == poly_vertices[nn-1] ) 
      // 	  poly_vertices.push_back( p1 ) ;
      // 	else{ cout << "WARNING: conflict in polygon building, at n="
      // 		   << nn <<endl; }

      
      ++nn;
      ++edge;
    } while ( edge != first);


    int NN = poly_vertices.size();
    if( NN == 0 ) continue;

    vvP poly_vertices2;

    //    poly_vertices2.push_back( poly_vertices[0] );

    for( int i= 0 ; i < NN ; i++ ) {

      Vector_2 dd( poly_vertices[ (i + 1) % NN ] - poly_vertices[ i ] );

      if( dd.squared_length() > threshold2 )
	poly_vertices2.push_back( poly_vertices[i] );

    }
    
    Polygon poly( poly_vertices2.begin() , poly_vertices2.end() );

    // if (!poly.is_convex() ) {
    //    cout << "WARNING: non-convex polygon" <<endl;
    //    cout << fv->point().point() << endl ;
       
    //    for( auto pp : poly_vertices2 ) // C++11
    // 	 cout << pp << endl; 
    // }

    //    int idx = fv->idx() ; // debugging
    
    fv->set_poly( poly );

   //    FT
    //area = poly.area();
    //fv->vol.set( area );

    FT area = fv->vol.val();
    
    if( area < threshold2 ) continue;
    
    Point c2 = CGAL::centroid(poly_vertices2.begin() , poly_vertices2.end(),
				CGAL::Dimension_tag<0>());
    
    fv->centroid.set( c2 );

    totalV += area;

    Point p_i = fv->point().point();
    
    fv->I.set( moi( p_i  ,  poly_vertices2  ) );

    Vector_2 d1 =  c2 - p_i ;

    Vector_2 d1A = area * d1;
 
    fv->dd.set( d1A );

    fv->dd2.set(  d1A.squared_length()  );

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
