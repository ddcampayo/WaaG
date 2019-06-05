// pParticles
// weigths are obtained dynamically, similarly to the pressure

#include"pParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {

  // TODO: read parameter file
  
  const int init_iter = 100;
  const FT  init_tol2 = 1e-4;

  const int inner_iter= 20;
  const FT  inner_tol = 1e-5;

  const FT total_time =  1/( 2 * 3.14 * 0.2) ;

  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

  //simu.do_perturb(0.01);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  linear algebra( T );

  // Init loop!
  
  int iter=0;

  for( ; iter < init_iter ; ++iter) {
  
    volumes( T ); 

    copy_weights( T ) ;

    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << iter << " dd = " << dd << endl;
    if( dd < init_tol2) break;

  }

  copy_weights( T ) ;

  cout << "Init loop converged in " << iter << " steps " << endl;
  
  set_vels_Gresho( T );

  volumes( T ); 
  algebra.copy( sfield_list::vol,  sfield_list::vol0);
 
  FT d0;
  FT dt=0.001;

  cin >> dt ;

  simu.set_dt( dt );
  // half-step (for e.g. leapfrog)
  FT dt2 = dt / 2.0 ;

  //  algebra.solve_for_weights();

  draw( T , particle_file     );

  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  do {
    simu.next_step();
    simu.advance_time( );

    backup( T );

    algebra.reset_p();
 
    copy_weights( T ) ;
    
    //  volumes( T );
    //  algebra.fill_Delta();
  
    algebra.u_star( );

    int iter = 1;

    FT displ = 0; // move( T , dt2 , d0 );

    // half-step corrector loop
    for ( ; iter <= inner_iter ; iter++) {

      displ = move( T , dt , d0 );

      move_weights( T );

      //    volumes( T );

      cout
	<< "********" << endl
	<< "Iter  " << iter
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
	<< endl ;

      volumes( T ); 
 
      algebra.fill_Delta_DD( dt );
 
      // algebra.solve_for_weights();
      //      copy_weights( T ) ;

      algebra.w_equation();

      algebra.copy( sfield_list::w,  sfield_list::p);
    
      algebra.u_add_press_grad( dt2 );

      if( displ < inner_tol ) break;

    }
    //    algebra.u_add_press_grad( dt2 );

    cout
      << "Whole step  "
      << " : disp " << displ << endl ;

    //algebra.w_equation();
    //algebra.solve_for_weights();
    algebra.u_add_press_grad( dt );

    algebra.solve_for_weights();

    //    volumes( T ); 
    //    update_full_vel( T );

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << " iters = " << iter
      << " L2_vel =  " << L2_vel_Gresho(T)
      << endl ;

    
  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
