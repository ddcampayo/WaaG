#include"linear.h"

void linear::fill_Delta_DD(void){ 

  std::cout << " Filling Delta _and_ capital D matrices" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet> aa, ax, ay, mx, my ; // , bb ;            // list of non-zeros coefficients

  typedef std::map<int,FT> diag_map;
  diag_map   dd; // diagonal
  diag_map   dd_x, dd_y;
  diag_map   dm_x, dm_y;

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

    Vector_2 rr_ij_j = pj - bij;
    Vector_2 rr_ij_i = pi - bij;

    Vector_2 DDij = Aij / lij * rr_ij_j; // ( pj - bij);
    Vector_2 DDji = Aij / lij * rr_ij_i; // ( pi - bij);

    FT r2_ij_j = rr_ij_j.squared_length();  // (these two are the same on Voronoi)
    FT r2_ij_i = rr_ij_i.squared_length();
    
    // todo: maybe define I = Aij*Aij/12, to ease notation
    Vector_2 MMij = Aij / lij * (
				 ( r2_ij_j +  Aij*Aij / 4 ) * rr_ij_j
				 - ( Aij*Aij / 12 ) * eij
				 // i.e.  + ( Aij*Aij / 12 ) * eji
				 );

    Vector_2 MMji = Aij / lij * (
				 ( r2_ij_i +  Aij*Aij / 4 ) * rr_ij_i
				 + ( Aij*Aij / 12 ) * eij
				 );


    if( (i >= 0 ) && ( j >= 0) ) {
      aa.push_back( triplet( i, j,  ddelta ));
      aa.push_back( triplet( j, i,  ddelta ));

      ax.push_back( triplet( i, j,  DDij.x() ));
      ay.push_back( triplet( i, j,  DDij.y() ));

      ax.push_back( triplet( j, i,  DDji.x() ));
      ay.push_back( triplet( j, i,  DDji.y() ));

      mx.push_back( triplet( i, j,  MMij.x() ));
      my.push_back( triplet( i, j,  MMij.y() ));

      mx.push_back( triplet( j, i,  MMji.x() ));
      my.push_back( triplet( j, i,  MMji.y() ));

    }

    // diagonal terms

    if (i >= 0 ) {
      dd[ i ]  -= ddelta;
//      dd_x[ i ] -= DDij.x();
//      dd_y[ i ] -= DDij.y();
      dd_x[ i ] -= DDji.x();
      dd_y[ i ] -= DDji.y();

//      dm_x[ i ] -= MMij.x();
//      dm_y[ i ] -= MMij.y();
      dm_x[ i ] -= MMji.x();
      dm_y[ i ] -= MMji.y();
 
    }
    if (j >= 0 ) {
      dd[ j ]  -= ddelta;
//      dd_x[ j ] -= DDji.x();
//      dd_y[ j ] -= DDji.y();
      dd_x[ j ] -= DDij.x();
      dd_y[ j ] -= DDij.y();

//      dm_x[ j ] -= MMji.x();
//      dm_y[ j ] -= MMji.y();
      dm_x[ j ] -= MMij.x();
      dm_y[ j ] -= MMij.y();

    }

    // if( (i!=0) && (j!=0) ) {
    //   aa.push_back( triplet(i - 1 , j -1 ,  ddelta ));
    //   aa.push_back( triplet(j - 1 , i -1 ,  ddelta ));

    //   ++N;
    // }

    if( i+1 > N ) { N = i+1 ; } // keep maximum

    //    cout << i << "  " << j << "  " << ddelta << endl;
  }



  // include "spring" term in M
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();

    if( idx < 0 ) continue;

    FT vol = fv->vol();
    Point ri = fv->point().point();
    Point bi = fv->centroid.val();

    Vector_2 dd = 2 * vol * ( bi - ri ) ;
    dm_x[ idx ] -= dd.x();
    dm_y[ idx ] -= dd.y();

  }

  
  // Add diagonal terms .-
  
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


  for( diag_map::const_iterator it = dm_x.begin(); it != dm_x.end(); ++it ) {
    int i    = it->first;
    FT diagx = it->second;
    mx.push_back( triplet( i, i,  diagx ));
  }

  
  for( diag_map::const_iterator it = dm_y.begin(); it != dm_y.end(); ++it ) {
    int i    = it->first;
    FT diagy = it->second;
    my.push_back( triplet( i, i,  diagy ));
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

  MMx.resize( N , N );
  MMx.setFromTriplets(mx.begin(), mx.end());

  MMy.resize( N , N );
  MMy.setFromTriplets(my.begin(), my.end());
  

  // set up solvers .-
  
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


  NN =
    - MMx * inv_vol.asDiagonal() * MMx.transpose()
    - MMy * inv_vol.asDiagonal() * MMy.transpose();

  NN_solver.compute( NN );

  if( NN_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing NN matrix " << endl;
  }
  
  
  return;

}

