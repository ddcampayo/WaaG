// pParticles
// no weigths are used, basically

#include"pParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {


  // TODO: read parameter file
  
  const int init_iters = 0;
  const FT  init_tol2 = 1e-3;

  const int inner_iters= 100;
  const FT  inner_tol  = 1e-6;

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

  for( ; iter < init_iters ; ++iter) {
  
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
  
  FT d0;
  FT dt=0.001;

  cin >> dt ;

  FT dt2 = dt / 2.0 ;

  simu.set_dt( dt );

  //  algebra.solve_for_weights();

  draw( T , particle_file     );

  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  volumes( T ); 

  do {
    simu.next_step();
    simu.advance_time( );

    volumes( T ); 

    backup( T );
    
    copy_weights( T ) ;
    
    //  volumes( T );
    //  algebra.fill_Delta();
  
    int iter = 1;

    // half-step corrector loop
    for ( ; iter <= inner_iters ; iter++) {

      FT displ = move( T , dt2 , d0 );

      cout
	<< "********" << endl
	<< "Iter  " << iter
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
	<< endl ;

      if( displ < inner_tol ) break;

      volumes( T ); 

      algebra.fill_Delta_DD();

      algebra.u_star( );

      algebra.p_equation( dt2 );

      algebra.u_add_press_grad( dt2 );
      
      //algebra.w_equation();
      //algebra.solve_for_weights();
      //      volumes( T ); 

      
    }
    //    algebra.u_add_press_grad( dt2 );

    copy_weights( T ) ;

    FT displ = move( T , dt , d0 );

    cout
      << "Whole step  "
      << " : disp " << displ << endl ;

    //algebra.w_equation();
    //algebra.solve_for_weights();
    volumes( T ); 

    update_full_vel( T );

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
