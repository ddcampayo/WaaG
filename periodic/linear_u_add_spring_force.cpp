#include"linear.h"
#include"fields_enum.h"

void linear::u_add_spring_force( const FT k , const FT dt ) {

  // TODO: this can be done through .dd, a vrtx member function
  // calculated in volumes.cpp
  
  // first, build a displacement vector
  
  typedef std::pair<int,FT> idx_val;
  std::vector<idx_val> idx_vals_x;
  std::vector<idx_val> idx_vals_y;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();

    if( idx < 0 ) continue;

    Point b_i = fv->centroid.val();

    Point r_i = fv->point().point();

    Vector_2 d_i = r_i - b_i ; //displacement from center-of-mass

    idx_vals_x.push_back( idx_val( idx , d_i.x()  ) );
    idx_vals_y.push_back( idx_val( idx , d_i.y()  ) );

  }

  int N = idx_vals_x.size();
  VectorXd disp_x( N );
  VectorXd disp_y( N );

  for(int nn=0; nn < N ; nn++)  {
    disp_x( idx_vals_x[nn].first ) = idx_vals_x[nn].second;
    disp_y( idx_vals_y[nn].first ) = idx_vals_y[nn].second;
  }


  VectorXd f_x = - k * disp_x;
  VectorXd f_y = - k * disp_y;

  // minus, to be interpreted as an "extra pressure gradient"
  vctrs_add_to_vfield(-f_x , -f_y , vfield_list::gradp );

  // velocity before forces are applied
  VectorXd U_x, U_y;

  // add spring force to u^star
  //  vfield_to_vctrs(  vfield_list::Ustar , U_x, U_y );

  // add spring force to u0  
  //  vfield_to_vctrs( vfield_list::U0  , U_x, U_y );

  // add force to u
  //vfield_to_vctrs(  vfield_list::U , U_x, U_y );
  //add force to u_star

  vfield_to_vctrs(  vfield_list::Ustar , U_x, U_y );

  VectorXd Unew_x, Unew_y;

  Unew_x = U_x + dt * f_x;
  Unew_y = U_y + dt * f_y;

  vctrs_to_vfield( Unew_x, Unew_y , vfield_list::Ustar );
  //  vctrs_to_vfield( Unew_x, Unew_y , vfield_list::U );

}


void linear::u_add_fem_force( const FT k , const FT dt ) {

  // TODO: this can be done through .dd, a vrtx member function
  // calculated in volumes.cpp
  
  // first, build a displacement vector
  
  typedef std::pair<int,FT> idx_val;
  std::vector<idx_val> idx_vals_x;
  std::vector<idx_val> idx_vals_y;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();

    if( idx < 0 ) continue;

    Vector_2 d_i = fv->FEM_disp.val();

    idx_vals_x.push_back( idx_val( idx , d_i.x()  ) );
    idx_vals_y.push_back( idx_val( idx , d_i.y()  ) );

  }

  int N = idx_vals_x.size();
  VectorXd disp_x( N );
  VectorXd disp_y( N );

  for(int nn=0; nn < N ; nn++)  {
    disp_x( idx_vals_x[nn].first ) = idx_vals_x[nn].second;
    disp_y( idx_vals_y[nn].first ) = idx_vals_y[nn].second;
  }

  VectorXd vol  = field_to_vctr( sfield_list::Dvol ) ;

  VectorXd f_x = k * disp_x.array() / vol.array()  ;
  VectorXd f_y = k * disp_y.array() / vol.array()  ;

  // minus, to be interpreted as an "extra pressure gradient"
  vctrs_add_to_vfield(-f_x , -f_y , vfield_list::gradp );

  // velocity before forces are applied
  VectorXd U_x, U_y;

  // add spring force to u^star
  //  vfield_to_vctrs(  vfield_list::Ustar , U_x, U_y );

  // add spring force to u0  
  //  vfield_to_vctrs( vfield_list::U0  , U_x, U_y );

  // add force to u
  //vfield_to_vctrs(  vfield_list::U , U_x, U_y );
  //add force to u_star

  vfield_to_vctrs(  vfield_list::Ustar , U_x, U_y );

  VectorXd Unew_x, Unew_y;

  Unew_x = U_x + dt * f_x;
  Unew_y = U_y + dt * f_y;

  vctrs_to_vfield( Unew_x, Unew_y , vfield_list::Ustar );
  //  vctrs_to_vfield( Unew_x, Unew_y , vfield_list::U );

}
