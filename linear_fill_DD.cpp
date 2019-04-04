#include"linear.h"

void linear::fill_Delta_DD(void){ 

  std::cout << " Filling Delta _and_ capital D matrices" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet> aa, ax, ay; // , bb ;            // list of non-zeros coefficients

  typedef std::map<int,FT> diag_map;
  diag_map   dd; // diagonal
  diag_map   ddx;
  diag_map   ddy;

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

    Point bij = CGAL::midpoint( Vor_segment );

    Point qj = vj->centroid();

    FT DDx =  Aij / lij * ( qj.x() - bij.x() );
    FT DDy =  Aij / lij * ( qj.y() - bij.y() );
    
    if( (i >= 0 ) && ( j >= 0) ) {
      aa.push_back( triplet( i, j,  ddelta ));
      aa.push_back( triplet( j, i,  ddelta ));

      aax.push_back( triplet( i, j,  DDx ));
      aax.push_back( triplet( j, i,  DDx ));

      aay.push_back( triplet( i, j,  DDy ));
      aay.push_back( triplet( j, i,  DDy ));
      
    }

    if (i >= 0 ) {
      dd[ i ]  -= ddelta;
      ddx[ i ] -= DDx;
      ddy[ i ] -= DDy;
    }
    if (j >= 0 ) {
      dd[ j ]  -= ddelta;
      ddx[ j ] -= DDx;
      ddy[ j ] -= DDy;
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

  for( diag_map::const_iterator it = ddx.begin(); it != ddx.end(); ++it ) {
    int i    = it->first;
    FT diagx = it->second;
    aax.push_back( triplet( i, i,  diagx ));
  }

  for( diag_map::const_iterator it = ddy.begin(); it != ddy.end(); ++it ) {
    int i    = it->first;
    FT diagy = it->second;
    aay.push_back( triplet( i, i,  diagy ));
  }

  Delta.resize( N , N );

  Delta.setFromTriplets(aa.begin(), aa.end());
  std::cout << " Filled Delta  matrix" << std::endl;
  cout << "matrix size " << Delta.rows() << " times " << Delta.cols() << endl;

  DDx.resize( N , N );

  DDx.setFromTriplets(aax.begin(), aax.end());
  std::cout << " Filled DDx  matrix" << std::endl;
  cout << "matrix size " << DDx.rows() << " times " << DDx.cols() << endl;

  DDy.resize( N , N );

  DDy.setFromTriplets(aay.begin(), aay.end());
  std::cout << " Filled DDy  matrix" << std::endl;
  cout << "matrix size " << DDy.rows() << " times " << DDy.cols() << endl;
  
  
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

  DDx_solver.compute( DDx );

  if(DDx_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing DDx " << //minus 1
      " matrix\n";
  }

  DDy_solver.compute( DDy );

  if(DDy_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing DDy " << //minus 1
      " matrix\n";
  }

  

  
  return;

}

