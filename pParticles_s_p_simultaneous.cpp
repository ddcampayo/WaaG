// pParticles
// no weigths are used, basically
// Pressure is used in order to enforce incompressibility
// An additional field, s, is used to enforce const moments of inertia

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

  //  simu.do_perturb(0.3);
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
    
    FT displ_p = 0 , displ_s = 0 ;
    
    for ( ; in_iter <= inner_max_iters ; in_iter++) {

      displ_p = move( T , dt , d0 );

      cout
	<< "********" << endl
	<< "P S Iter  " << in_iter
	<< " . Moved from previous (rel.): " << displ_p <<
	" ; from original (rel.): " << d0
	<< endl ;
	
      volumes( T ); 

      algebra.fill_Delta_DD();

      algebra.p_equation( dt );
      algebra.s_equation( dt );

      algebra.u_add_grads( dt2 );

      if(  displ_p < disp_tol ) break;

    }
    ///////////// end   p iter
    
    cout
      << "Whole step  "
      << " : disp " <<  displ_p << endl ;

    algebra.u_add_grads( dt );

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
