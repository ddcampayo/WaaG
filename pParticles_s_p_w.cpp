// pParticles
// weigths are used to enforce incompressibility
// Pressure is determined some other way...
// An additional field, s, is used to enforce const moments of inertia

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

  int s_iters; cin >> s_iters; //= 10;
  int p_iters; cin >> p_iters; //= 10;

  FT total_time =  1/( 2 * 3.14 * 0.2) ;

  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

  simu.do_perturb(0.1);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  linear algebra( T );

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

  // whole step
  FT dt2 = dt / 2.0 ;

  //  algebra.solve_for_weights();

  draw( T , particle_file     );

  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  do {
    simu.next_step();
    simu.advance_time( );


    int in_iter = 1 , s_it = 1 , p_it= 1;
    
    //    volumes( T ); 
    
    backup( T );

    //   displ = move( T , dt2 , d0 );

    algebra.reset_s();
    algebra.reset_p();

    algebra.u_star( );
    
    FT displ = 0 ;

    
    for ( ; in_iter <= inner_max_iters ; in_iter++) {


      algebra.solve_for_weights();

      copy_weights( T ) ;

      volumes( T ); 

      algebra.fill_Delta_DD();

      algebra.p_equation( dt );

      algebra.u_add_press_grad( dt2 );
	
      displ = move( T , dt , d0 );

      cout
	<< "********" << endl
	<< "Inner Iter  " << p_it
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
	<< endl ;
	
	//	algebra.u_add_grads( dt2 );
	
      algebra.s_equation( dt );	

      algebra.u_add_s_grad( dt2 );

      if( displ < disp_tol ) break;
 
    }

    
    cout
      << "Whole step  "
      << " : disp " << displ << endl ;

    volumes( T ); 

    // half-step:
//    update_full_vel( T );

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << " iters = " << in_iter
      << " T =  " << kinetic_E(T)
      << " L2_vel =  " << L2_vel_Gresho(T)
      << endl ;



  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
