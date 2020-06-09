// Function to fill the various matrices involved. All of
// them. May be split in the future


// Glossary:

// Delta_ij = d V_i / d w_j   (inf.  change of volume of cell i due to change in weight of cell j)

// DD_ij =  d V_i / d r_j  (inf.  change of volume of cell i due to change in position of cell j)  , a matrix of vectors, thus stored as DD_ij_x and DD_ij_y

// LL =  DD (1/V) DD^t, involved in Ralphson-Newton methods

// Gamma_ij = d I_i / d w_j , aka GG  (inf.  change of moment of inertia of cell i due to change in weight of cell j)

// MM_ij = d I_i / d r_j   (inf.  change of moment of inertia of cell i due to change in position of cell j)

// NN =  MM (1/V) MM^t, involved in Ralphson-Newton methods


#include"linear.h"


void linear::fill_Delta_DD( const FT dt ) {

  //  std::cout << " Filling Delta _and_ capital D matrices" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet>
    aa, ax, ay, mx, my,
    ee ,
    gg ; // , bb ;            // list of non-zeros coefficients

  typedef std::map<int,FT> diag_map;
  diag_map   dd; // diagonal
  diag_map   dd_g;
  diag_map   dd_e;
  diag_map   dd_x, dd_y;
  diag_map   dm_x, dm_y;

  int N=1;

  F_e_it eit = T.finite_edges_begin();

  for ( ; eit !=T.finite_edges_end(); ++eit) {

    Face_handle f =  eit -> first ;
    int i0 = eit -> second;

#ifdef FEM
    int i3 = i0 ;
    Vertex_handle v3  = f->vertex( i3 );
    Point p3 = v3->point().point() ;
    
    Vertex_handle v1 = f->vertex( (i3+1) % 3);
    Point p1 = v1->point().point() ;

    Vertex_handle v3p = T.mirror_vertex( f , i3 );
    Point p3p = v3p->point().point() ;

    Vector_2 v_3p_3 = p3 - p3p ;

    // negative right angle turn
    Vector_2 v_3p_3_perp = Vector_2( v_3p_3.y() , -v_3p_3.x() );

    //    Triangle tr( v1->point().point() , v3p->point().point() , v3->point().point() );

    //    CGAL::Orientation ori = tr.orientation();
    //    if( ori == CGAL::NEGATIVE ) v33_perp = -v33_perp;

    Vector_2 v_1_3p = p3p - p1 ;

    
    CGAL::Orientation ori = CGAL::orientation(  v_1_3p , v_3p_3 );
    if( ori == CGAL::RIGHT_TURN ) v_3p_3_perp = -v_3p_3_perp;

    Vector_2 DDij = v_3p_3_perp / 6.0 ;
    Vector_2 DDji = -DDij;

#endif
    
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

    // experimental
    //    Point mij = pi + 0.5 * eij;  // mid-point between i and j
    
    FT lij2 =  eij.squared_length() ;
    FT lij = std::sqrt( lij2 );
  
    FT ddelta = - 0.5 * Aij / lij;

    Point bij = CGAL::midpoint( Vor_segment->source() , Vor_segment->target() );

    // nope!
    //    Point qj = vj->centroid();

    Vector_2 rr_ij_j = pj - bij;
    Vector_2 rr_ij_i = pi - bij;

#ifndef FEM
    // regular
    Vector_2 DDij = Aij / lij * rr_ij_j; // ( pj - bij);
    Vector_2 DDji = Aij / lij * rr_ij_i; // ( pi - bij);
    //    //experimental
    //    Vector_2 DDij = Aij / lij * (pj - mij) ; // ( pj - bij);
    //    Vector_2 DDji = Aij / lij * (pi - mij) ; // ( pi - bij);
#endif
    
    // dd**2
    FT Eij = Aij / lij * ( vi->dd.val() * rr_ij_i );
    FT Eji = Aij / lij * ( vj->dd.val() * rr_ij_j );

    // dd.x only:
    //    FT Eij = Aij / lij * rr_ij_i.x() / 2.0 ;
    //    FT Eji = Aij / lij * rr_ij_j.x() / 2.0 ;
    
    FT r2_ij_j = rr_ij_j.squared_length();  // (these two are the same on Voronoi)
    FT r2_ij_i = rr_ij_i.squared_length();

    Vector_2 rr_ij_j_perp =   ( ( rr_ij_j * eij  ) / lij2 ) * eij ;
    Vector_2 rr_ij_j_para =     rr_ij_j  -  rr_ij_j_perp;

    Vector_2 rr_ij_i_perp =   ( ( rr_ij_i * eij  ) / lij2 ) * eij ;
    Vector_2 rr_ij_i_para =     rr_ij_i  -  rr_ij_i_perp; // ==  rr_ij_j_para, actually !
 
    FT I = Aij*Aij/12 ;
    
    // Voronoi-only:
    // Vector_2 MMij = Aij / lij * (
    // 				 ( r2_ij_j +  Aij*Aij / 4 ) * rr_ij_j
    // 				 - ( Aij*Aij / 12 ) * eij
    // 				 );
    // // debug 
    // // cout
    // //   << " MMij :  "
    // //   <<  Aij*Aij / 4  * rr_ij_j - ( Aij*Aij / 12 ) * eij
    // //   << "  "
    // //   <<  Aij*Aij / 12  * rr_ij_j + ( Aij*Aij / 6 ) * rr_ij_i_para
    // //   << endl;


    // Vector_2 MMji = Aij / lij * (
    // 				 ( r2_ij_i +  Aij*Aij / 4 ) * rr_ij_i
    // 				 + ( Aij*Aij / 12 ) * eij
    // 				 );

    
    // General:
    Vector_2 MMij = Aij / lij * (
    				 ( r2_ij_i +  I ) * rr_ij_j
    				 + 2 *I * rr_ij_i_para
    				 );

    Vector_2 MMji = Aij / lij * (
    				 ( r2_ij_j +  I ) * rr_ij_i
    				 + 2* I * rr_ij_j_para
    				 );

    Vector_2 MMii =-Aij / lij * (
    				 ( r2_ij_i +  I ) * rr_ij_i
    				 + 2 * I * rr_ij_i_para
    				 );

    Vector_2 MMjj =-Aij / lij * (
    				 ( r2_ij_j +  I ) * rr_ij_j
    				 + 2 * I * rr_ij_j_para
    				 );

    FT gamma_ij = ddelta * ( I + r2_ij_i  );
    FT gamma_ji = ddelta * ( I + r2_ij_j  );

    if( (i >= 0 ) && ( j >= 0) ) {
      aa.push_back( triplet( i, j,  ddelta ));
      aa.push_back( triplet( j, i,  ddelta ));

      gg.push_back( triplet( i, j,  gamma_ij ));
      gg.push_back( triplet( j, i,  gamma_ji ));

      ax.push_back( triplet( i, j,  DDij.x() ));
      ay.push_back( triplet( i, j,  DDij.y() ));

      ax.push_back( triplet( j, i,  DDji.x() ));
      ay.push_back( triplet( j, i,  DDji.y() ));

      mx.push_back( triplet( i, j,  MMij.x() ));
      my.push_back( triplet( i, j,  MMij.y() ));

      mx.push_back( triplet( j, i,  MMji.x() ));
      my.push_back( triplet( j, i,  MMji.y() ));

      ee.push_back( triplet( i, j,  Eij ) );
      ee.push_back( triplet( j, i,  Eji ) );
      
    }

    // diagonal terms

    if (i >= 0 ) {
      dd[ i ]  -= ddelta;

      dd_g[ i ]  -= gamma_ji ;

      dd_e[ i ]  -= Eij ;
      
//      dd_x[ i ] -= DDij.x();
//      dd_y[ i ] -= DDij.y();

      dd_x[ i ] -= DDji.x();
      dd_y[ i ] -= DDji.y();

      // Voronoi-only:

      // dm_x[ i ] -= MMji.x();
      // dm_y[ i ] -= MMji.y();
      
      // General:
      
      dm_x[ i ] += MMii.x();
      dm_y[ i ] += MMii.y();
 
    }

    if (j >= 0 ) {
      dd[ j ]  -= ddelta;

      dd_g[ j ]  -= gamma_ij ;

      dd_e[ j ]  -= Eji ;

//      dd_x[ j ] -= DDji.x();
//      dd_y[ j ] -= DDji.y();
      dd_x[ j ] -= DDij.x();
      dd_y[ j ] -= DDij.y();

      // Voronoi-only:

      // dm_x[ j ] -= MMij.x();
      // dm_y[ j ] -= MMij.y();

      // General:

      dm_x[ j ] += MMjj.x();
      dm_y[ j ] += MMjj.y();

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

    FT vol = fv->vol.val();
    Point ri = fv->point().point();
    Point bi = fv->centroid.val();

    Vector_2 dd = 2 * vol * ( bi - ri ) ;

    // signs ???
    dm_x[ idx ] += dd.x();
    dm_y[ idx ] += dd.y();

  }

  
  // Add diagonal terms .-
  
  for( diag_map::const_iterator it = dd.begin(); it != dd.end(); ++it ) {
    int i = it->first;
    FT diag = it->second;
    aa.push_back( triplet( i, i,  diag ));

    //    cout << i << "  " << i << "  " << diag << endl;

  }


  for( auto it : dd_g ) {
    int i = it.first;
    FT diag = it.second;
    gg.push_back( triplet( i, i,  diag ));
  }


  for( auto it : dd_e ) {
    int i = it.first;
    FT diag = it.second;
    ee.push_back( triplet( i, i,  diag ));
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

  GG.resize( N , N );
  GG.setFromTriplets(gg.begin(), gg.end());

  EE.resize( N , N );
  EE.setFromTriplets(ee.begin(), ee.end());

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
  
  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;
  VectorXd inv_vol  = 1.0 / vol.array() ;

  LL =
    - DDx * inv_vol.asDiagonal() * DDx.transpose()
    - DDy * inv_vol.asDiagonal() * DDy.transpose();

  LL_solver.compute( LL );

  if( LL_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing LL matrix " << endl;
  }

  // special.-  experimental

  if( dt > 1e-8) Delta -= dt*dt*LL;

  Delta_solver.compute( Delta );

  if(Delta_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing Delta " << //minus 1
      " matrix\n";
  }

  GG = GG.transpose();

  GG_solver.compute( GG );

  if(GG_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing Gamma " << //minus 1
      " matrix\n";
  }

  
  NN =
    - MMx * inv_vol.asDiagonal() * MMx.transpose()
    - MMy * inv_vol.asDiagonal() * MMy.transpose();

  NN_solver.compute( NN );

  if( NN_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing NN matrix " << endl;
  }
  

  EE_solver.compute( EE );

  if(EE_solver.info()!=Eigen::Success)
    std::cout << "Failure decomposing Epsilon matrix\n";
  
  return;

}

