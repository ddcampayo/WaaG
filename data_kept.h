struct data_kept {
  int idx;
  Point pos;  // some position
  Point r0;
  weight w;
  weight w0;
  FT vol0;
  Vector_2 Dr;
  Vector_2 U;
  Vector_2 U0;
  Vector_2 Ustar;
  FT p;
  
  data_kept(const F_v_it fv) {
    idx = fv->idx();
    r0 = fv->r0.val();
    vol0 = fv->vol0.val();
    w = fv->w.val();
    w0 = fv->w0.val();
    Dr = fv->Dr.val();
    U = fv->U.val();
    U0 = fv->U0.val();
    Ustar = fv->Ustar.val();
    p = fv->p.val();
  }

  void restore(Vertex_handle fv) {
    fv->idx.set( idx );
    fv->r0.set( r0 );
    fv->vol0.set( vol0 );
    fv->w0.set( w0 );
    fv->w.set( w );
    fv->Dr.set( Dr );
    fv->U.set( U );
    fv->U0.set( U0 );
    fv->Ustar.set( Ustar );
    fv->p.set( p );
  }

};
