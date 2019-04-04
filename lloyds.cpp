#include"pParticles.h"
#include"simu.h"
#include"data_kept.h"


// Lloyds algorithm
//
//
FT lloyds(Triangulation& T) {

  cout << "Lloyds algorithm ... " << endl ;
  
  //  copy_weights( T );

  // volumes( T );

  //  copy_weights( T );

  vector<data_kept> prev;

  FT dd2=0;

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

    Point rnow=fv->point().point(); // current point

    Point rnew= fv->centroid.val();

    Vector_2 disp2 = rnew - rnow;

    FT rel_disp = sqrt(disp2.squared_length() ) / simu.h();

    dd2 += rel_disp;

    data.pos = rnew;

    data.Dr = disp2;

    prev.push_back (data);

  }

  dd2 /= simu.no_of_particles();

  cout << "Moved relative displacement: " <<
    sqrt(dd2)/simu.h()   ;


  cout << " . Mean displacement: " << dd2 << endl ;

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

