struct data_kept {
  int idx;
  Point pos;  // some position
  Point r0;
  weight w, w0;
  FT vol0;
  Vector_2 Dr;
  Vector_2 U, U0;
  Vector_2 Ustar;
  FT p, p0;

  FT s;
  FT I0;

  Vector_2 dd;
  FT dd2;

  FT om, om0;
  
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
    p0 = fv->p0.val();

    s  = fv->s.val();
    I0 = fv->I0.val();

    dd  =  fv->dd.val();
    dd2 = fv->dd2.val(); ;

    om = fv->om.val();
    om0 = fv->om0.val();

  }

  void restore(Vertex_handle fv) {
    fv->idx.set( idx );

    fv->r0.set( r0 );

    fv->vol0.set( vol0 );

    fv->w.set( w );
    fv->w0.set( w0 );

    fv->Dr.set( Dr );

    fv->U.set( U );
    fv->U0.set( U0 );
    fv->Ustar.set( Ustar );

    fv->p.set( p );
    fv->p0.set( p0 );

    fv->s.set( s );
    fv->I0.set( I0 );

    fv->dd.set( dd ) ;
    fv->dd2.set( dd2 ); ;

    fv->om.set( om );
    fv->om0.set( om0 );

  }

};
