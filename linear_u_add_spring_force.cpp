#include"linear.h"
#include"fields_enum.h"

void linear::u_add_spring_force( const FT kdt ) {

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

  VectorXd Ustar_x, Ustar_y;

  // add spring force to u^star
  //  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  // add force to u
  vfield_to_vctrs(  vfield_list::U , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  U_x = Ustar_x.array() - kdt * disp_x.array();
  U_y = Ustar_y.array() - kdt * disp_y.array();

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}
