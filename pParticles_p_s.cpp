// pParticles
// no weigths are used, basically
// Pressure is used in order to enforce incompressibility
// An additional field, s, is used to enforce const moments of inertia

#include"pParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {

  // TODO: read parameter file
  
  const int init_max_iters = 0;
  const FT  init_tol2 = 1e-3;

  const int inner_max_iters = 10;
  const FT  disp_tol  = 1e-6;

  const int s_iters= 10;  
  const int p_iters= 4;
  const FT total_time = 2 * M_PI * 0.2 ; // one whole turn

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

  int init_iter=1;
  
  for( ; init_iter <= init_max_iters ; ++init_iter) {
  
    volumes( T ); 

    copy_weights( T ) ;

    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << init_iter << " dd = " << dd << endl;
    if( dd < init_tol2) break;

  }

  copy_weights( T ) ;

  cout << "Init loop converged in " << init_iter << " steps " << endl;
  
  set_vels_Gresho( T );

  volumes( T ); 
  
  FT d0;
  FT dt=0.001;

  cin >> dt ;

  simu.set_dt( dt );

  // half-step leapfrog
  //  FT dt2 = dt / 2.0 ;

  // whole step
  FT dt2 = dt  ;

  //  algebra.solve_for_weights();

  draw( T , particle_file     );

  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  do {
    simu.next_step();
    simu.advance_time( );

    FT displ;

    int in_iter = 1 , s_it = 1 , p_it= 1;
    
    for ( ; in_iter <= inner_max_iters ; in_iter++) {

      //   moment of inertia (s)   iteration

      volumes( T ); 

      backup( T );

      algebra.u_star( );
    
      s_it = 1;
     
      displ = move( T , dt2 , d0 );

      algebra.reset_s();

      for ( ; s_it <= s_iters ; s_it++) {
	volumes( T ); 

	algebra.fill_Delta_DD();

	algebra.s_equation( dt2 );

	// whole step, special 1st time
	if( simu.current_step() == 1 ){
	  algebra.u_add_s_grad( dt / 2 );
	}  else 
	  {
	    algebra.u_add_s_grad( dt2 );
	  }

	displ = move( T , dt2 , d0 );

	cout
	  << "********" << endl
	  << "S Iter  " << s_it
	  << " . Moved from previous (rel.): " << displ <<
	  " ; from original (rel.): " << d0
	  << endl ;

	//      if( displ < disp_tol ) break;

      }

      volumes( T ); 

      copy_weights( T ) ;

      //      displ = move( T , dt2 , d0 );

      ///////////// end   s iter


      //    //   volume (p)   iteration
    
      //    volumes( T ); 

      //    backup( T );
    
      //    copy_weights( T ) ;
    
    //  volumes( T );
    //  algebra.fill_Delta();


      volumes( T ); 

      backup( T );

      algebra.u_star( );
    
      p_it = 1;
     
      displ = move( T , dt2 , d0 );

      algebra.reset_p();

      // predictor loop
      for ( ; p_it <= p_iters ; p_it++) {
	volumes( T ); 

	algebra.fill_Delta_DD();

	algebra.p_equation( dt2 );

	// whole step, special 1st time
	if( simu.current_step() == 1 ){
	  algebra.u_add_press_grad( dt / 2 );
	}  else 
	  {
	    algebra.u_add_press_grad( dt2 );
	  }

	displ = move( T , dt2 , d0 );

	cout
	  << "********" << endl
	  << "P Iter  " << p_it
	  << " . Moved from previous (rel.): " << displ <<
	  " ; from original (rel.): " << d0
	  << endl ;

      }
 
      if( displ < disp_tol ) break;

      //algebra.w_equation();
      //algebra.solve_for_weights();
      //      volumes( T ); 

    }
    ///////////// end   p iter

    
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
      << " L2_vel =  " << L2_vel_Gresho(T)
      << endl ;

  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
