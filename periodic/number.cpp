#include"pParticles.h"
#include"simu.h"

void number(Triangulation& T) {

  int idx=0;

  int N=simu.no_of_particles();

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    if( vit->idx() < 0) continue;

    vit->idx.set( idx );

    ++idx;

  }

  return;
}
