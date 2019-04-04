// pParticles
// weigths are obtained dynamically, similarly to the pressure

#include"pParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {

  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

  //ysimu.do_perturb(0.01);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  linear algebra( T );

  // Init loop!
  
  const int max_iter = 100;
  const FT tol2 = 1e-3;
  int iter=0;

  for( ; iter < max_iter ; ++iter) {
  
    volumes( T ); 

    copy_weights( T ) ;

    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << iter << " dd = " << dd << endl;
    if( dd < tol2) break;

  }

  cout << "Init loop converged in " << iter << " steps " << endl;
  
  set_vels_Gresho( T );

  volumes( T ); 

  
  FT d0;
  FT dt=0.001;

  cin >> dt ;

  FT dt2 = dt / 2.0 ;

  simu.set_dt( dt );

  FT total_time = 1/( 2 * 3.14 * 0.2) ;

  algebra.solve_for_weights();

  std::ofstream log_file;
  log_file.open("main.log");

  do {
    simu.next_step();
    simu.advance_time( );

    backup( T );

    volumes( T ); 

    //  volumes( T );
    //  algebra.fill_Delta();
  
    algebra.u_star( );

    int iter = 0;

    copy_weights( T ) ;

    // half-step corrector loop
    for ( ; iter < 40 ; iter++) {
      cout << "Iter  " << iter << "  ";
  
      FT displ = move( T , dt2 , d0 );

      cout << " : disp " << displ << endl ;

      volumes( T ); 

      algebra.p_equation( dt2 );

      algebra.u_add_press_grad( dt2 );

      //algebra.w_equation();
      algebra.solve_for_weights();
      //volumes( T ); 

      if( displ < 1e-8) break;

    }

    FT displ = move( T , dt , d0 );
    //    algebra.w_equation();
    //algebra.solve_for_weights();
    volumes( T ); 

    update_full_vel( T );
    
    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  " << " iters  "
      << iter
      << endl ;

    
  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
