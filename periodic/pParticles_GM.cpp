// pParticles
// no weigths are used, basically
// Pressure is used in order to enforce incompressibility

// half-step froggy predictor-corrector

#include"pParticles.h"
#include"linear.h"
#include"simu.h"


sim_data simu;

int main() {


  // TODO: read parameter file
  
  const int init_iters = -1;// 100;// 1000;
  const FT  init_tol2 = 1e-5;

  const int inner_iters= 10;
  const FT  inner_tol  = 1e-5;
  const  FT turn_time = 1 ;
  
  //  const  FT total_time = turn_time; // once

  const  FT total_time = 3 * turn_time; // three whole turns
  
  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

  simu.do_perturb(1e-2);
  //simu.do_perturb(0);
  create( T , 1.0 );
  number( T );
  expand( T , 1.0 );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );
  volumes( T ); 
  linear algebra( T );
  algebra.copy( sfield_list::vol,  sfield_list::vol0);

  set_vels_TG( T );

  // Init loop!
  
  int iter=1;

  for( ; iter < init_iters ; ++iter) {
  
    volumes( T ); 

    //    copy_weights( T ) ;

    //    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << iter << " dd = " << dd << endl;
    if( dd < init_tol2) break;

  }

  // volumes( T ); 
  // simu.set_dt( 0 );  
  // draw( T , particle_file     );
  // draw_diagram( T , diagram_file );  
  // return 0;

  copy_weights( T ) ;

  cout << "Init loop converged in " << iter << " steps " << endl;
  
  set_vels_TG( T );

  volumes( T ); 
  algebra.copy( sfield_list::vol,  sfield_list::vol0);
  
  FT d0;
  FT dt=0.001;

  cout << "Time step, dt = ";
  cin >> dt ;
  cout << endl << dt << endl;

  simu.set_dt( dt );
  
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
  
  // half-step (for e.g. leapfrog)
  FT dt2 = dt / 2.0 ;

  // whole step
  // FT dt2 = dt  ;

  //  algebra.solve_for_weights();

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");
  log_file << " #  step   time   iters   kin_energy   L2_velocity " << endl;

  do {
    simu.next_step();
    simu.advance_time( );

    backup( T );
    
    //    copy_weights( T ) ;
    
    //  volumes( T );
    //  algebra.fill_Delta();

    //    algebra.reset_s();
  
    int iter = 1;

//    algebra.fill_Delta_DD();

    algebra.u_star( );

    FT displ = 0; // move( T , dt2 , d0 );
    
    displ = move( T , dt , d0 );

    // frog
    //      displ = move( T , dt2 , d0 );

    algebra.u_star( );

    cout
      << "********" << endl
      << "Iter  " << iter
      << " . Moved from previous (rel.): " << displ <<
      " ; from original (rel.): " << d0
      << endl ;

    algebra.solve_for_weights();
    algebra.copy( 0.5 * spring ,  sfield_list::w ,  sfield_list::p);

    algebra.clear_vfield( vfield_list::gradp );

    algebra.u_add_spring_force( spring ,  dt );

    algebra.copy( vfield_list::Ustar ,  vfield_list::U );

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << iter-1 << " "
      << kinetic_E(T) << " "
      << L2_vel_TG(T) << " "
      << endl ;

  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
