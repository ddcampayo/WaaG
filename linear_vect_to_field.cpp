#include"linear.h"
#include"fields_enum.h"

VectorXd linear::field_to_vctr(const sfield_list::take sf ) {

  typedef std::pair<int,FT> idx_val;
  std::vector<idx_val> idx_vals;
  //  std::vector<int> indices;

  //  std::vector<FT> ff;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();

    if( idx < 0 ) continue;

    idx_vals.push_back( idx_val( idx , fv->sfield(sf).val() ) );
    //    ff.push_back( fv->sfield(sf).val() );

  }

  int N = idx_vals.size();
  VectorXd vctr( N );

  //  vctr.resize(N);

  for(int nn=0; nn < N ; nn++)  {

    vctr( idx_vals[nn].first ) = idx_vals[nn].second;

    // cout << indices[nn] << " : "
    //  	 << ff[nn] << endl;
  }

  return vctr;

}


void linear::vfield_to_vctrs(const vfield_list::take vf , VectorXd& vx, VectorXd& vy ) {

  typedef std::pair<int,FT> idx_val;

  std::vector<idx_val> idx_vals_x;
  std::vector<idx_val> idx_vals_y;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();

    if( idx < 0 ) continue;

    Vector_2 vv = fv->vfield(vf).val();
    idx_vals_x.push_back( idx_val( idx , vv.x() ) );
    idx_vals_y.push_back( idx_val( idx , vv.y() ) );

  }

  int N = idx_vals_x.size();
  vx.resize( N );
  vy.resize( N );

  //  vctr.resize(N);

  for(int nn=0; nn < N ; nn++)  {

    vx( idx_vals_x[nn].first ) = idx_vals_x[nn].second;
    vy( idx_vals_y[nn].first ) = idx_vals_y[nn].second;

  }

  return;

}


void linear::vctr_to_field(const VectorXd& vv, const sfield_list::take sf  ) {

  //  int nn=0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();
    if( idx < 0 ) continue;
      
    // check (only makes sense if vertices are re-labeled at every step)
    // if(index!=nn) {
    //   std::cout << "Error transfering from vector onto field " << scalarf << std::endl;
    // }

    // cout << idx << " ";
    // cout << vv[ idx ] << endl;
    
    fv->sfield(sf).set( vv[ idx ] );

    //    ++nn;

    // cout << index << " : "
    // 	 << vv[index] << endl;

  }
  return;
}


void linear::vctrs_to_vfield(const VectorXd& vx,
			   const VectorXd& vy,
			   const vfield_list::take vf  ) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int idx = fv->idx.val();
    if( idx < 0 ) continue;
      
    fv->vfield(vf).set( Vector_2 ( vx[ idx ] , vy[idx] ) );

  }
  return;
  
}

