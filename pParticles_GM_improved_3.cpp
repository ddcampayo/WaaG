// pParticles
// Attempt to replicate

// Gallouët, T.O., Mérigot, Q.
// A Lagrangian Scheme à la Brenier for the Incompressible Euler Equations.
// Found Comput Math 18, 835–865 (2018).
// https://doi.org/10.1007/s10208-017-9355-y
// Simple explicit method

#include"pParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {
  
  const int init_iters = 20;
  const FT  init_tol2 = 1e-6;

  const int inner_iters= 10;
  const FT  inner_tol  = 1e-5;

  const  FT turn_time = 2 * M_PI * 0.2 ; // one whole turn
  
  //  const  FT total_time = turn_time; // once

  const  FT total_time = 2 * turn_time; // twice

  
  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

  simu.do_perturb(1e-3);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  volumes( T ); 
  linear algebra( T );
  algebra.copy( sfield_list::vol,  sfield_list::vol0);
  algebra.copy( sfield_list::I,  sfield_list::I0);

  // Init loop!
  
  int iter=1;

  for( ; iter <= init_iters ; ++iter) {
  
    //    volumes( T ); 
    //    copy_weights( T ) ;

    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << iter << " dd = " << dd << endl;
    if( dd < init_tol2) break;

  }

  copy_weights( T ) ;

  set_vels_Gresho( T );

  cout << "Init loop converged in " << iter << " steps " << endl;

  volumes( T ); 
  algebra.copy( sfield_list::vol,  sfield_list::vol0);

  FT d0;
  FT dt=0.001;

  cout << "Time step, dt = ";
  cin >> dt ;
  cout << endl << dt << endl;

  simu.set_dt( dt );
  FT dt2 = dt / 2.0 ;
  
  // Setting a spring period that includes several Dt, in
  // order spring forces be properly sampled
  
  //  FT spring_period = 10 * dt;

  FT spring_to_dt;
  cout << "Spring period / dt  = ";
  cin >> spring_to_dt;
  cout << endl << spring_to_dt << endl;

  // 31 dt is the value for G&M first simulation,
  // "Beltrami flow in the square"
  FT spring_period = spring_to_dt * dt;
//  FT spring_period = 80 * dt;
  FT omega = 2 * M_PI /  spring_period ;

  cout << " omega  = " << omega << endl ;

  FT spring = omega*omega; // factor that appears in the spring force

  draw( T , particle_file);
  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");
  log_file << " #  step   time   iters   kin_energy   L2_velocity " << endl;


  do {
    simu.next_step();
    simu.advance_time( );

    backup( T );
    int iter = 1;

    volumes( T );

    algebra.u_star( );

    FT displ = 0;

    for ( ; iter <= inner_iters ; iter++) {

      displ = move( T , dt2 , d0 );

      algebra.u_star( );

      cout
	<< "********" << endl
	<< "Iter  " << iter
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
      << endl ;

      volumes( T ); 
      algebra.fill_Delta_DD();

      algebra.clear_vfield( vfield_list::gradp );

      algebra.solve_for_weights();
      
      algebra.u_add_spring_force( spring , dt2 );

      algebra.p_equation_lapl_div_source_fem( dt2 );

      algebra.u_add_press_grad_fem( dt2 );

      //      algebra.fill_Delta_DD();

    // pressure proportional to weights: p = 1/2 spring w
      //      algebra.copy( 0.5 * spring ,  sfield_list::w ,  sfield_list::p);

    //  copy_weights( T ) ;
    //  volumes( T );


      //      algebra.copy( vfield_list::Ustar ,  vfield_list::U );

      algebra.copy( vfield_list::Ustar ,  vfield_list::U );
      
      if( displ < inner_tol ) break;
    }

    displ = move( T , dt , d0 );

    update_half_velocity( T );
    
    volumes( T ); 

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << iter-1 << " "
      << kinetic_E(T) << " "
      << L2_vel_Gresho(T) << " "
      << endl ;


  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
