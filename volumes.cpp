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
//
//
void volumes(Triangulation& T) {

  const FT threshold=1e-10;
  const FT threshold2= threshold*threshold;
  FT totalV = 0;
  
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {
    //    fv->vol.reset();

    Edge_circulator edge  = T.incident_edges( fv );

    Edge_circulator first = edge;

    first--; // avoid last one entirely
    

    vvP poly_vertices;

    int nn=1;

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
    poly_vertices2.push_back( poly_vertices[0] );

    for( int i= 1 ; i < NN ; i++ ) {
      Vector_2 dd( poly_vertices[ i ] - poly_vertices[ i-1 ] );
      if(
	 dd.squared_length() > threshold2 )
	poly_vertices2.push_back( poly_vertices[i] );
    }
    
    Polygon poly( poly_vertices2.begin() , poly_vertices2.end() );

    // if (!poly.is_convex() ) {
    //    cout << "WARNING: non-convex polygon" <<endl;
    //    cout << fv->point().point() << endl ;
       
    //    for( auto pp : poly_vertices2 ) // C++11
    // 	 cout << pp << endl; 
    // }
       
    fv->set_poly( poly );

    FT area = poly.area();

    fv->vol.set( area );

    if( area < threshold2 ) continue;
    
    Point c2 = CGAL::centroid(poly_vertices2.begin() , poly_vertices2.end(),
				CGAL::Dimension_tag<0>());
    
    fv->centroid.set( c2 );

    totalV += area;

    fv->I.set( moi( fv->point().point() ,  poly_vertices2  ) );
    
  }

  cout << "Volumes: total = " << totalV << " ; ";

  FT  inner_V = 0;
  int inner   = 0;
  
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)
    if (fv->idx() > 0 ) {
      inner_V +=  fv->vol.val();
      ++inner;
    }

  cout << "inner = " << inner_V << " ; ";

  simu.set_totalV( totalV );
  simu.set_innerV( inner_V );

  cout << "mean inner V = " << simu.meanV() << endl;

  
}
