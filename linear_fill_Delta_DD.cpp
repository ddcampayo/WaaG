#include"linear.h"

void linear::fill_Delta_DD(void){ 

  std::cout << " Filling Delta _and_ capital D matrices" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet> aa, ax, ay; // , bb ;            // list of non-zeros coefficients

  typedef std::map<int,FT> diag_map;
  diag_map   dd; // diagonal
  diag_map   dd_x, dd_y;

  int N=1;

  F_e_it eit = T.finite_edges_begin();

  for ( ; eit !=T.finite_edges_end(); ++eit) {

    Face_handle f =  eit -> first ;
    int i0 = eit -> second;

    Vertex_handle vi = f->vertex( (i0+1) % 3);
    Vertex_handle vj = f->vertex( (i0+2) % 3);

    int i = vi->idx();
    int j = vj->idx();

    //    if( (i < 0 ) || ( j < 0) ) continue;
    
    CGAL::Object o = T.dual(eit);

    const Segment * Vor_segment = CGAL::object_cast<Segment>( &o );

    if (! Vor_segment ) continue;

    FT Aij = std::sqrt( Vor_segment->squared_length() );

    //    cout << " A0 = " << Aij << endl;
    
    Point pi = vi->point().point();
    Point pj = vj->point().point();

    Vector_2 eij = pj - pi;
    
    FT lij = std::sqrt( eij.squared_length() );
    
    FT ddelta = - 0.5 * Aij / lij;

    Point bij = CGAL::midpoint( Vor_segment->source() , Vor_segment->target() );

    // nope!
    //    Point qj = vj->centroid();

    Vector_2 DDij = Aij / lij * ( pj - bij);
    Vector_2 DDji = Aij / lij * ( pi - bij);
    
    if( (i >= 0 ) && ( j >= 0) ) {
      aa.push_back( triplet( i, j,  ddelta ));
      aa.push_back( triplet( j, i,  ddelta ));

      ax.push_back( triplet( i, j,  DDij.x() ));
      ay.push_back( triplet( i, j,  DDij.y() ));

      ax.push_back( triplet( j, i,  DDji.x() ));
      ay.push_back( triplet( j, i,  DDji.y() ));
    }

    if (i >= 0 ) {
      dd[ i ]  -= ddelta;
      dd_x[ i ] -= DDij.x();
      dd_y[ i ] -= DDij.y();
    }
    if (j >= 0 ) {
      dd[ j ]  -= ddelta;
      dd_x[ j ] -= DDji.x();
      dd_y[ j ] -= DDji.y();
    }

    // if( (i!=0) && (j!=0) ) {
    //   aa.push_back( triplet(i - 1 , j -1 ,  ddelta ));
    //   aa.push_back( triplet(j - 1 , i -1 ,  ddelta ));

    //   ++N;
    // }

    if( i+1 > N ) { N = i+1 ; } // keep maximum

    //    cout << i << "  " << j << "  " << ddelta << endl;
  }



  for( diag_map::const_iterator it = dd.begin(); it != dd.end(); ++it ) {
    int i = it->first;
    FT diag = it->second;
    aa.push_back( triplet( i, i,  diag ));

    //    cout << i << "  " << i << "  " << diag << endl;

  }

  for( diag_map::const_iterator it = dd_x.begin(); it != dd_x.end(); ++it ) {
    int i    = it->first;
    FT diagx = it->second;
    ax.push_back( triplet( i, i,  diagx ));
  }

  for( diag_map::const_iterator it = dd_y.begin(); it != dd_y.end(); ++it ) {
    int i    = it->first;
    FT diagy = it->second;
    ay.push_back( triplet( i, i,  diagy ));
  }

  Delta.resize( N , N );

  Delta.setFromTriplets(aa.begin(), aa.end());
  // std::cout << " Filled Delta  matrix" << std::endl;
  // cout << "matrix size " << Delta.rows() << " times " << Delta.cols() << endl;

  DDx.resize( N , N );

  DDx.setFromTriplets(ax.begin(), ax.end());
  // std::cout << " Filled DDx  matrix" << std::endl;
  // cout << "matrix size " << DDx.rows() << " times " << DDx.cols() << endl;

  DDy.resize( N , N );

  DDy.setFromTriplets(ay.begin(), ay.end());
  // std::cout << " Filled DDy  matrix" << std::endl;
  // cout << "matrix size " << DDy.rows() << " times " << DDy.cols() << endl;

  Delta_solver.compute( Delta );

  if(Delta_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing Delta " << //minus 1
      " matrix\n";
  }

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;
  VectorXd inv_vol  = 1.0 / vol.array() ;
  
  LL =
    - DDx * inv_vol.asDiagonal() * DDx.transpose()
    - DDy * inv_vol.asDiagonal() * DDy.transpose();

  LL_solver.compute( LL );

  if( LL_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing LL matrix " << endl;
  }

  
  return;

}

