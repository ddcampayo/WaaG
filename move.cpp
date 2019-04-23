#include"pParticles.h"
#include"simu.h"
#include"data_kept.h"


void copy_weights( Triangulation& T ) {
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) {
    //    int idx = fv->idx();

    //    if(idx < 0 ) continue;

    fv->w.set( fv->point().weight() ); 

  }
}




void backup( Triangulation& T ) {
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) {
    fv->vol0.set( fv->vol.val()  );
    fv->w0.set(   fv->w()  );
    fv->p0.set(   fv->p()  );
    fv->I0.set(   fv->I()  );
    fv->U0.set(   fv->U()  );
    fv->r0.set(   fv->point().point() );
  }
}


// // now in linear::
// void u_star(Triangulation& T) {
//   for(F_v_it fv=T.finite_vertices_begin();
//       fv!=T.finite_vertices_end();
//       fv++) {
//     fv->Ustar.set( fv->U.val()  );
//   }

// }


void update_full_vel( Triangulation& T ) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) {
    Vector_2 U0    = fv->U0.val();
    Vector_2 Uhalf = fv->U.val();
    fv->U.set(  2 * Uhalf - U0  );
  }
}


FT move(Triangulation& T, const FT dt , FT& dd0 ) {

  //  cout << "Moving nodes ... " << endl;
  
  copy_weights( T );
  
  vector<data_kept> prev;

  FT dd2=0;

  //  bool first=true;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) {

    int idx = fv->idx();
    
    data_kept data(fv);

    if(idx < 0 ) {
      
      data.pos = fv->point().point(); 

      prev.push_back (data);

      continue;

    }
    
    Vector_2  vel = fv->U();

    Vector_2 disp = dt * vel;

    Point rnow=fv->point().point(); // current point

    Point r0=fv->r0(); // starting point

    Point rnew= r0 + disp;

    Vector_2 disp2 = rnew - rnow;

    FT rel_disp = sqrt(disp2.squared_length() ) / simu.h();

    FT rel_disp0= sqrt( disp.squared_length() ) / simu.h();

    // if(first) {
    //  cout
    //    << "r0 " << r0 << "  "
    //    << "rnow " << rnow << "  "
    //    << "rnew " << rnew << "  "
    //    << "disp2 " << disp2 << "  "
    //    << " idx " << fv->idx() << " "
    //    << "rel_disp " << rel_disp
    //    << endl ;
    //  first=false;
    // }

    dd2 += rel_disp;

    dd0 += rel_disp0;

    data.pos = rnew;

    data.Dr = disp;

    prev.push_back (data);

  }

  // cout << " moved. Relative displacement: " <<
  //   sqrt(dd2)/simu.no_of_particles()/simu.h()   ;

  dd2 /= simu.no_of_particles()*simu.h();
  dd0 /= simu.no_of_particles()*simu.h();

  //  cout << " . Mean displacement: " << dd2 << endl ;

  T.clear(); // clears the triangulation !!

  for(vector<data_kept>::iterator data=prev.begin();
      data!=prev.end();
      data++) {

    //    cout << "Inserting back at " << data->pos << endl ;
    
    Vertex_handle fv=T.insert(   wPoint( data->pos , data->w )  );

    data->restore(fv);

  }

  return dd2;
}



FT move_from_centroid(Triangulation& T, const FT dt ) {

  cout << "Moving nodes from centroids ... " << endl;
  
  copy_weights( T );
  
  vector<data_kept> prev;

  FT dd=0;

  //  bool first=true;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) {

    int idx = fv->idx();
    
    data_kept data(fv);

    if(idx < 0 ) {
      
      data.pos = fv->point().point(); 

      prev.push_back (data);

      continue;

    }
    
    Vector_2  vel = fv->U();

    Vector_2 disp = dt * vel;

    Point r0 = fv->centroid.val();

    Point rnew= r0 + disp;

    FT rel_disp = sqrt(disp.squared_length() ) / simu.h();

    dd += rel_disp;

    data.pos = rnew;

    data.Dr = disp;

    prev.push_back (data);

  }
  
  cout << " moved. Relative displacement: " <<
    sqrt(dd)/simu.no_of_particles()/simu.h()   ;

  dd /= simu.no_of_particles();

  cout << " . Mean displacement: " << dd << endl ;

  T.clear(); // clears the triangulation !!

  for(vector<data_kept>::iterator data=prev.begin();
      data!=prev.end();
      data++) {

    //    cout << "Inserting back at " << data->pos << endl ;
    
    Vertex_handle fv=T.insert(   wPoint( data->pos , data->w )  );

    data->restore(fv);

  }

  return dd;
}





void move_weights( Triangulation& T )
{

  vector<data_kept> prev;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) {

    //    int idx = fv->idx();

    //    if(idx < 0 ) continue;
    
    data_kept data(fv);

    Point p = fv->point().point();

    data.pos = p;

    // FT new_w = fv->w();
    
    // data.w = new_w; // overwrite previous value!
    
    prev.push_back (data);

  }


  T.clear(); // clears the triangulation !!

  for(vector<data_kept>::iterator data=prev.begin();
      data!=prev.end();
      data++) {

    //    cout << "Inserting back at " << data->pos << endl ;

    Vertex_handle fv=T.insert(   wPoint( data->pos , data->w )  );

    data->restore(fv);

  }

}
