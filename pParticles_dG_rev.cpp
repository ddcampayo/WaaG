// pParticles
// Attempt go beyond de Goes  et al.
// Power Particles: An incompressible fluid solver based on power diagrams

#include"pParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {

  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;
  
  cout << "Creating point cloud" << endl;

//  simu.do_perturb(0.1);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  linear algebra( T );


  // Init loop!
  
  const int init_max_iter = 100;
  const FT init_tol2 = 1e-6;

  const int inner_iters= 100;
  const FT  inner_tol  = 1e-6;

  int iter=0;


  for( ; iter < init_max_iter ; ++iter) {
  
    volumes( T ); 

    copy_weights( T ) ;

    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << iter << " dd = " << dd << endl;

    if( dd < init_tol2) break;

  }

  //algebra.solve_for_moments();

  cout << "Init loop converged in " << iter << " steps " << endl;
  
  set_vels_Gresho( T );

  volumes( T ); 

  algebra.copy( sfield_list::vol,  sfield_list::vol0);
  algebra.copy( sfield_list::I  ,  sfield_list::I0);

  FT d0;
  FT dt=0.001;

  cin >> dt ;

  simu.set_dt( dt );

  FT dt2 = dt / 2.0;
  
//  move_from_centroid( T , dt);

  draw( T , particle_file);
  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  FT total_time = 2 * M_PI * 0.2 ; // one whole turn

  do {
    simu.next_step();
    simu.advance_time( );

    backup( T );
      
    int iter = 1;

    algebra.u_star( );

    FT displ = 0; // move( T , dt2 , d0 );
    
    // full-step corrector loop

    for ( ; iter <= inner_iters ; iter++) {

      displ = move( T , dt , d0 );

      //move_from_centroid( T , dt);

      cout
	<< "********" << endl
	<< "Iter  " << iter
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
	<< endl ;

      //      volumes( T ); 

      //      copy_weights( T ) ;

      //algebra.solve_for_weights( dt );

      algebra.solve_for_moments();

      //algebra.fill_Delta_DD();
      // //      algebra.w_equation2();
      // algebra.w_equation3(); 
      // move_weights( T );
      //volumes( T ); 
      //algebra.fill_Delta_DD();
      //      copy_weights( T ) ;

      algebra.s_equation( dt );

      //algebra.p_equation( dt );

      algebra.u_add_s_grad( dt2 );
      //      algebra.u_add_press_grad( dt2 );
      //      algebra.u_add_press_grad_MM_w( dt2 );
      
      //algebra.u_add_grads( dt2 );

      //algebra.u_add_w_grad( dt2 );

      if( displ < inner_tol ) break;

    }

    algebra.solve_for_moments();
    
    cout
      << " ======= " << endl
      << "Whole step  "
      << " : disp = " << displ
      << " ; d0 = " << d0
      << endl 
      << " ======= " << endl;


    //    volumes( T ); 
      
    //    algebra.fill_Delta_DD();
    
    //    algebra.u_add_grads( dt );

    //    algebra.u_add_press_grad_MM_w( dt );
    //    algebra.u_add_press_grad( dt );
    algebra.u_add_s_grad( dt );

    //    algebra.u_add_w_grad( dt );

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << " iters = " << iter
      << " T =  " << kinetic_E(T)
      << " L2_vel =  " << L2_vel_Gresho(T)
      << endl ;

    
  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
