// pParticles
// weigths are used to enforce incompressibility
// Pressure is determined some other way...
// An additional field, s, is used to enforce const moments of inertia

// s is obtained from p here

#undef PRESSURE_PPE
#include"pParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {

  // TODO: read better from parameter file
  
  int init_max_iters; cin >> init_max_iters; // = 40;
  FT  init_tol2 ; cin >> init_tol2 ; // = 1e-3; 

  int inner_max_iters; cin >> inner_max_iters; // = 10; 
  FT  disp_tol; cin >> disp_tol; //  = 1e-6;

  const  FT turn_time = 2 * M_PI * 0.2 ; // one whole turn
  
  //  const  FT total_time = turn_time; // once

  const  FT total_time = 2 * turn_time; // twice
  
  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

//  simu.do_perturb(0.1);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );
  volumes( T ); 
  linear algebra( T );
  algebra.copy( sfield_list::vol,  sfield_list::vol0);
  algebra.copy( sfield_list::I,  sfield_list::I0);

  // Init loop!

  int init_iter=0;
  
  for( ; init_iter < init_max_iters ; ++init_iter) {
  
    volumes( T ); 

    //    copy_weights( T ) ;
    //    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << init_iter << " dd = " << dd << endl;
    if( dd < init_tol2) break;

  }

  volumes( T ); 

  algebra.copy( sfield_list::vol,  sfield_list::vol0);
  algebra.copy( sfield_list::I,  sfield_list::I0);

  copy_weights( T ) ;
  
  //  copy_weights( T ) ;

  cout << "Init loop converged in " << init_iter << " steps " << endl;
  
  set_vels_Gresho( T );

  volumes( T );

  algebra.copy( sfield_list::vol,  sfield_list::vol0);
  algebra.copy( sfield_list::I,  sfield_list::I0);
  
  FT d0;
  FT dt=0.001;

  cin >> dt ;

  simu.set_dt( dt );

  // half-step leapfrog
  //  FT dt2 = dt / 2.0 ;

  FT dt2 = dt / 2.0 ;

  //  algebra.solve_for_weights();

  draw( T , particle_file     );

  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");
  log_file << " #  step   time   iters   kin_energy   L2_velocity " << endl;

  do {
    simu.next_step();
    simu.advance_time( );


    cout << "Time  " << simu.time() << endl;
    //    volumes( T ); 
    
    backup( T );

    //   displ = move( T , dt2 , d0 );

    //    algebra.reset_s();
    //   algebra.reset_p();

    algebra.u_star( );
 
    FT displ = 0 ;

    int in_iter = 1;
        
    for ( ; in_iter <= inner_max_iters ; in_iter++) {

      displ = move( T , dt2 , d0 );

      algebra.u_star( );

      cout
	<< "********" << endl
	<< "Iter  " << in_iter
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
	<< endl ;

      volumes( T ); 

      algebra.fill_Delta_DD();

      //      algebra.solve_for_weights();

      //algebra.solve_for_moments();
      //copy_weights( T ) ;
      //volumes( T ); 
      //algebra.fill_Delta_DD();

      algebra.clear_vfield( vfield_list::gradp );

      //algebra.p_equation_s( dt );
      
      //algebra.p_equation( dt );
      //algebra.p_equation_divgrad_div_source( dt2 );

      algebra.s_equation_from_p( 1 );

      algebra.u_add_s_grad( dt2 );

      algebra.p_equation_lapl_div_source( dt2 );

      //algebra.u_add_grads( dt2 );
      
      algebra.u_add_press_grad( dt2 );
      algebra.copy( vfield_list::Ustar ,  vfield_list::U );

      if( displ < disp_tol ) break;
 
    }
    
    //    algebra.u_star( );
//    algebra.p_equation_from_s( );

//    algebra.p_equation( dt );

//    algebra.u_add_press_grad( dt );

    displ = move( T , dt , d0 );

    update_half_velocity( T );

    volumes( T ); 
    
    cout
      << "Whole step  "
      << " : disp " << displ << endl ;



    // half-step:
//    update_full_vel( T );

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    
    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << in_iter-1 << " "
      << kinetic_E(T) << " "
      << L2_vel_Gresho(T) << " "
      << endl ;

  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
