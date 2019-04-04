#include"linear.h"

void linear::fill_Delta(void){ 

  std::cout << " Filling Delta matrix" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet> aa; // , bb ;            // list of non-zeros coefficients

  typedef std::map<int,FT> diag_map;
  diag_map   dd; // diagonal

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

    if( (i >= 0 ) && ( j >= 0) ) {
      aa.push_back( triplet( i, j,  ddelta ));
      aa.push_back( triplet( j, i,  ddelta ));
    }

    if (i >= 0 )    dd[ i ] -= ddelta;
    if (j >= 0 )    dd[ j ] -= ddelta;

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
  
  Delta.resize( N , N );

  Delta.setFromTriplets(aa.begin(), aa.end());
  std::cout << " Filled Delta  matrix" << std::endl;
  cout << "matrix size " << Delta.rows() << " times " << Delta.cols() << endl;

  //  bool s = saveMarket(Delta, "Delta");
  //bool s = saveMarket( MatrixXd( Delta ), "Delta");
  //  cout << MatrixXd( Delta );

//   // Only non-direct iterative solvers
// #ifndef DIRECT_SOLVER
//   solver_stiffp1.setTolerance( TOL );
//   solver_stiffp1.setMaxIterations( 40 * N );
// #endif

  Delta_solver.compute( Delta );

  if(Delta_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing Delta " << //minus 1
      " matrix\n";
  }

  return;

}

