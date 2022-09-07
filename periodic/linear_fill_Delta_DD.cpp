// Function to fill the various matrices involved. All of
// them. May be split in the future

// Glossary:

// Delta_ij = d V_i / d w_j   (inf.  change of volume of cell i due to change in weight of cell j)

// DD_ij =  d V_i / d r_j, aka "L"  (inf.  change of volume of cell i due to change in position of cell j), a matrix of vectors, thus stored as DD_ij_x and DD_ij_y

// LL = - DD (1/V) DD^t, involved in Ralphson-Newton methods ( notice sign )

// Gamma_ij = d I_i / d w_j , aka "GG"  (inf.  change of moment of inertia of cell i due to change in weight of cell j)

// MM_ij = d I_i / d r_j , aka "N"  (inf.  change of moment of inertia of cell i due to change in position of cell j), a matrix of vectors, thus stored as MM_ij_x and MM_ij_y

// More info at the end of this file !!


#include"linear.h"


void linear::fill_Delta_DD( const FT dt ) {

  //  std::cout << " Filling Delta _and_ capital D matrices" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet>
    aa, aa0 ,ax, ay, mx, my,
    ax_fem, ay_fem,
    ee ,
    gg ; // , bb ;            // list of non-zeros coefficients

  typedef std::map<int,FT> diag_map;
  diag_map   dd, dd0; // diagonal
  diag_map   dd_g;
  diag_map   dd_e;
  diag_map   dd_x, dd_y;
  diag_map   dd_x_fem, dd_y_fem;
  diag_map   dm_x, dm_y;

  int N=1;

  F_e_it eit = T.finite_edges_begin();

  for ( ; eit !=T.finite_edges_end(); ++eit) {

    Face_handle f =  eit -> first ;
    int i0 = eit -> second;

    //   #ifdef FEM

    Vector_2 DDij_fem, DDji_fem ;

    FT ddelta0;  // Voronoi ddelta (regardless of weights!!!)

    {
      int i3 = i0 ;
      Vertex_handle v3  = f->vertex( i3 );
      Point p3 = v3->point().point() ;
    
      Vertex_handle v1 = f->vertex( (i3+1) % 3);
      Point p1 = v1->point().point() ;

      Vertex_handle v2 = f->vertex( (i3+2) % 3);
      Point p2 = v2->point().point() ;

      Vertex_handle v3p = T.mirror_vertex( f , i3 );
      Point p3p = v3p->point().point() ;

      Vector_2 v_3p_3 = p3 - p3p ;

      // negative right angle turn
      Vector_2 v_3p_3_perp = Vector_2( v_3p_3.y() , -v_3p_3.x() );

      //    Triangle tr( v1->point().point() , v3p->point().point() , v3->point().point() );

      //    CGAL::Orientation ori = tr.orientation();
      //    if( ori == CGAL::NEGATIVE ) v33_perp = -v33_perp;

      // Vector_2 v_1_3p = p3p - p1 ;

      // CGAL::Orientation ori = CGAL::orientation(  v_1_3p , v_3p_3 );
      // if( ori == CGAL::RIGHT_TURN ) v_3p_3_perp = -v_3p_3_perp;

      DDij_fem = v_3p_3_perp / 6.0 ;

      ///////
      //experimental: project onto line connecting i and j:
      // does not seem to work too great, perhaps because ...

      // Vector_2 e12 = p2 - p1; // == eij later

      // FT l12_2 =  e12.squared_length() ;
      
      // DDij_fem = ( ( DDij_fem * e12 ) / l12_2 ) * e12;

      ///////

      // ... interaction respects 3rd law
      
      DDji_fem = -DDij_fem ;

      Point p0 = p3;

      FT ll0= CGAL::squared_distance( p1 , p2 );
      FT ll1= CGAL::squared_distance( p2 , p0 );
      FT ll2= CGAL::squared_distance( p0 , p1 );

      FT area = CGAL::area(p0, p1, p2) ; // std::fabs( tr.area() );

      Point p0p = p3p;
      FT ll1p= CGAL::squared_distance( p0p , p2 );
      FT ll2p= CGAL::squared_distance( p1  , p0p );

      FT areap = CGAL::area(p0p, p2, p1) ;
      
      ddelta0 =
	( ll1  + ll2  - ll0 ) / ( 8 * area  )+
	( ll1p + ll2p - ll0 ) / ( 8 * areap );

    }
    
    // #endif
    
    Vertex_handle vi = f->vertex( (i0+1) % 3);
    Vertex_handle vj = f->vertex( (i0+2) % 3);

    // equiv. v1 and v2 before !
    
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

    //#ifndef FEM
    // regular
    Vector_2 DDij = Aij / lij * rr_ij_j; // ( pj - bij);
    Vector_2 DDji = Aij / lij * rr_ij_i; // ( pi - bij);
    //    //experimental
    //    Vector_2 DDij = Aij / lij * (pj - mij) ; // ( pj - bij);
    //    Vector_2 DDji = Aij / lij * (pi - mij) ; // ( pi - bij);
    //experimental, equivalent: project onto line connecting i and j:
    DDij = ( (DDij*eij) / lij2 ) * eij;
    DDji = ( (DDji*eij) / lij2 ) * eij;
    //#endif
    
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


    //// WORK IN PROGRESS!!
    if(i >= 0 ) {

      int j0 =  vj->idx0();  // takes care of periodic BCs

      aa.push_back( triplet( i, j,  ddelta ));
      aa.push_back( triplet( j, i,  ddelta ));

      aa0.push_back( triplet( i, j,  ddelta0 ));
      aa0.push_back( triplet( j, i,  ddelta0 ));
      
      gg.push_back( triplet( i, j,  gamma_ij ));
      gg.push_back( triplet( j, i,  gamma_ji ));

      ax.push_back( triplet( i, j,  DDij.x() ));
      ay.push_back( triplet( i, j,  DDij.y() ));

      ax.push_back( triplet( j, i,  DDji.x() ));
      ay.push_back( triplet( j, i,  DDji.y() ));

      ax_fem.push_back( triplet( i, j,  DDij_fem.x() ));
      ay_fem.push_back( triplet( i, j,  DDij_fem.y() ));

      ax_fem.push_back( triplet( j, i,  DDji_fem.x() ));
      ay_fem.push_back( triplet( j, i,  DDji_fem.y() ));
      
      mx.push_back( triplet( i, j,  MMij.x() ));
      my.push_back( triplet( i, j,  MMij.y() ));

      mx.push_back( triplet( j, i,  MMji.x() ));
      my.push_back( triplet( j, i,  MMji.y() ));

      ee.push_back( triplet( i, j,  Eij ) );
      ee.push_back( triplet( j, i,  Eji ) );
      
    }

    
    if( (i >= 0 ) && ( j >= 0) ) {
      aa.push_back( triplet( i, j,  ddelta ));
      aa.push_back( triplet( j, i,  ddelta ));

      aa0.push_back( triplet( i, j,  ddelta0 ));
      aa0.push_back( triplet( j, i,  ddelta0 ));
      
      gg.push_back( triplet( i, j,  gamma_ij ));
      gg.push_back( triplet( j, i,  gamma_ji ));

      ax.push_back( triplet( i, j,  DDij.x() ));
      ay.push_back( triplet( i, j,  DDij.y() ));

      ax.push_back( triplet( j, i,  DDji.x() ));
      ay.push_back( triplet( j, i,  DDji.y() ));

      ax_fem.push_back( triplet( i, j,  DDij_fem.x() ));
      ay_fem.push_back( triplet( i, j,  DDij_fem.y() ));

      ax_fem.push_back( triplet( j, i,  DDji_fem.x() ));
      ay_fem.push_back( triplet( j, i,  DDji_fem.y() ));
      
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
      dd0[ i ]  -= ddelta0;

      dd_g[ i ]  -= gamma_ij ; // NOT gamma_ji (I think)

      dd_e[ i ]  -= Eij ;
      
//      dd_x[ i ] -= DDij.x();
//      dd_y[ i ] -= DDij.y();

      dd_x[ i ] -= DDji.x();
      dd_y[ i ] -= DDji.y();

      dd_x_fem[ i ] -= DDji_fem.x();
      dd_y_fem[ i ] -= DDji_fem.y();
      

      // Voronoi-only:

      // dm_x[ i ] -= MMji.x();
      // dm_y[ i ] -= MMji.y();
      
      // General:
      
      dm_x[ i ] += MMii.x();
      dm_y[ i ] += MMii.y();
 
    }

    if (j >= 0 ) {
      dd[ j ]  -= ddelta;
      dd0[ j ]  -= ddelta0;

      dd_g[ j ]  -= gamma_ji ;

      dd_e[ j ]  -= Eji ;

//      dd_x[ j ] -= DDji.x();
//      dd_y[ j ] -= DDji.y();

      dd_x[ j ] -= DDij.x();
      dd_y[ j ] -= DDij.y();

      dd_x_fem[ j ] -= DDij_fem.x();
      dd_y_fem[ j ] -= DDij_fem.y();

      
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



  // include "spring" term in M, dI_i / dr_i = ...  - 2 V_i (b_i - r_i)
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();

    if( idx < 0 ) continue;

    FT vol = fv->vol.val();
    Point ri = fv->point().point();
    Point bi = fv->centroid.val();

    Vector_2 dd = 2 * vol * ( ri - bi ) ;

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

  for( auto it : dd0 ) {
    int i = it.first;
    FT diag = it.second;
    aa0.push_back( triplet( i, i,  diag ));
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


  for( diag_map::const_iterator it = dd_x_fem.begin(); it != dd_x_fem.end(); ++it ) {
    int i    = it->first;
    FT diagx = it->second;
    ax_fem.push_back( triplet( i, i,  diagx ));
  }

  
  for( diag_map::const_iterator it = dd_y_fem.begin(); it != dd_y_fem.end(); ++it ) {
    int i    = it->first;
    FT diagy = it->second;
    ay_fem.push_back( triplet( i, i,  diagy ));
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

  Delta0.resize( N , N );
  Delta0.setFromTriplets(aa0.begin(), aa0.end());

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

  DDx_fem.resize( N , N );
  DDx_fem.setFromTriplets(ax_fem.begin(), ax_fem.end());

  DDy_fem.resize( N , N );
  DDy_fem.setFromTriplets(ay_fem.begin(), ay_fem.end());
  
  MMx.resize( N , N );
  MMx.setFromTriplets(mx.begin(), mx.end());

  MMy.resize( N , N );
  MMy.setFromTriplets(my.begin(), my.end());
  

  // set up solvers .-
  
  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;
  VectorXd inv_vol  = 1.0 / vol.array() ;


  // regular.-
   LL =
     - DDx * inv_vol.asDiagonal() * DDx.transpose()
     - DDy * inv_vol.asDiagonal() * DDy.transpose();

  // rotation.-

//  VectorXd I  = field_to_vctr( sfield_list::I ) ;
//  VectorXd inv_I  = 1.0 / I.array() ;

//  VectorXd x,y;
//  vfield_to_vctrs( vfield_list::dd , x , y ) ;

//  VectorXd x2 = x.array()*x.array();
//  VectorXd y2 = y.array()*y.array();
//  VectorXd xy = x.array()*y.array();

//  LL =
//    - DDx * inv_vol.asDiagonal() * DDx.transpose()
//    - DDy * inv_vol.asDiagonal() * DDy.transpose()
//    - DDx * inv_I.asDiagonal() * (
//				  y2.asDiagonal() * DDx.transpose()
//				  -xy.asDiagonal() * DDy.transpose() )
//    - DDy * inv_I.asDiagonal() * (
//				  x2.asDiagonal() * DDy.transpose()
//				  -xy.asDiagonal() * DDx.transpose() );
  
  LL_solver.compute( LL );

  if( LL_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing LL matrix " << endl;
  }

  // special.-  experimental

  if( dt > 1e-4) Delta -= dt*dt*LL;

  Delta_solver.compute( Delta );

  if(Delta_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing Delta " << //minus 1
      " matrix\n";
  }

  Delta0_solver.compute( Delta0 );

  if(Delta0_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing Delta0 " << //minus 1
      " matrix\n";
  }
  
  GG_solver.compute(  GG.transpose() );

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

  LN =
    - DDx * inv_vol.asDiagonal() * MMx.transpose()
    - DDy * inv_vol.asDiagonal() * MMy.transpose();

  NL =
    - MMx * inv_vol.asDiagonal() * DDx.transpose()
    - MMy * inv_vol.asDiagonal() * DDy.transpose();

  return;

}

// Related variables :

// aa     ->   Delta  (collector)
// ddelta ->   aa     (off-diagonal terms,  = dV_i/dw_j =dV_j/dw_i )
// dd     ->   aa     (diagonal terms,  = dV_i/dw_i)

// gg       ->   GG  (collector)
// gamma_ij ->   gg  (off-diagonal terms,  = dI_i/dw_j)
// gamma_ji ->   gg  (off-diagonal terms,  = dI_j/dw_i)
// dd_g     ->   gg  (diagonal terms,  = dI_i/dw_i)

// ax       ->   DDx (collector)
// DDij.x() ->   ax  (off-diagonal terms,  = dV_i/dr_j)
// DDji.x() ->   ax  (off-diagonal terms,  = dV_j/dr_i)
// dd_x     ->   ax  (diagonal terms,  = dV_i/dr_i)

// ay       ->   DDy (collector) , ... etc ...

// mx       ->   MMx (collector)
// MMij.x() ->   mx  (off-diagonal terms,  = dI_i/dr_j)
// MMji.x() ->   mx  (off-diagonal terms,  = dI_j/dr_i)
// dm_x     ->   mx   (diagonal terms,  = dI_i/dr_i)
// MMii.x() ->   dm_x (diagonal terms,  = dI_i/dr_i)
// MMjj.x() ->   dm_x (diagonal terms,  = dI_j/dr_j)
// spring   ->   dm_x (diagonal terms, spring-like force = int_i  d/dr_i |r-r_i|^2 = 2 V_i (b_i - r_i) )




// Solvers :
// notice GG solver solves for  u = GG^t  v !
// notice LL = - DD (1/V) DD^t ( sign! )


// Experimental (others are not documented):

// NN =  - MM (1/V) MM^t, involved in Ralphson-Newton methods

// LN = -  DD (1/V) MM^t, involved in Ralphson-Newton methods

// NL = -  MM (1/V) DD^t, involved in Ralphson-Newton methods

