// pParticles
// no weigths are used, basically
// Pressure is used in order to enforce incompressibility

#include"pParticles.h"
#include"linear.h"
#include"simu.h"

sim_data simu;

int main() {


  // TODO: read parameter file
  
  const int init_iters = 0;
  const FT  init_tol2 = 1e-4;

  const int inner_iters= 100;
  const FT  inner_tol  = 1e-6;

  const  FT total_time = 2 * M_PI * 0.2 ; // one whole turn

  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

  
  simu.do_perturb(1e-2);
  
  create( T , 1.0 );
  number( T );

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );

  simu.next_step();
  simu.advance_time( );

  
  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  linear algebra( T );
  
  // Init loop!
  
  int iter=1;

  volumes( T ); 
    
  for( ; iter < init_iters ; ++iter) {

    //    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << iter << " dd = " << dd << endl;

    //    volumes( T ); 


    volumes( T ); 

    copy_weights( T ) ;
    
    algebra.dd2_stats( ) ;

    if( dd < init_tol2) break;

  }

  copy_weights( T ) ;

  cout << "Init loop converged in " << iter << " steps " << endl;

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );

  simu.next_step();
  simu.advance_time( );

  algebra.copy( sfield_list::vol,  sfield_list::vol0);

  
  algebra.solve_for_weights_centroid();

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );


  simu.next_step();
  simu.advance_time( );

  algebra.solve_for_weights();

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );
  
  return 0;



  
  set_vels_Gresho( T );

  volumes( T ); 
  algebra.copy( sfield_list::vol,  sfield_list::vol0);
  
  FT d0;
  FT dt=0.001;

  cin >> dt ;

  simu.set_dt( dt );

  // half-step (for e.g. leapfrog)
  FT dt2 = dt / 2.0 ;

  // whole step
  // FT dt2 = dt  ;

  //  algebra.solve_for_weights();

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  // special first iter.-
  // cout << " First iter, free ";
  // algebra.reset_p();
  // FT displ1 = move( T , dt , d0 );  
  // volumes( T ); 
  // simu.next_step();  simu.advance_time();
  // draw( T , particle_file     );
  // draw_diagram( T , diagram_file );
 
  do {
    simu.next_step();
    simu.advance_time( );

    backup( T );
    
    //    copy_weights( T ) ;
    
    //  volumes( T );
    //  algebra.fill_Delta();

    algebra.reset_s();
  
    int iter = 1;

//    algebra.fill_Delta_DD();

    algebra.u_star( );

    FT displ = 0; // move( T , dt2 , d0 );
    
    // full-step corrector loop

    for ( ; iter <= inner_iters ; iter++) {

      displ = move( T , dt , d0 );

      cout
	<< "********" << endl
	<< "Iter  " << iter
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
	<< endl ;

      volumes( T ); 
      
      algebra.fill_Delta_DD();

      // // whole step, special 1st time
      // if( simu.current_step() == 1 ){
      // 	algebra.p_equation( dt2 ); 
      // 	algebra.u_add_press_grad( dt2/2  );

      // }  else 
      // {
      // 	algebra.p_equation( dt ); 
      // 	algebra.u_add_press_grad( dt2 );
      // }

      algebra.p_equation( dt ); 

      algebra.u_add_press_grad( dt2 );

      if( displ < inner_tol ) break;

      ////// testing ...
      //      backup( T );
      //algebra.reset_p();
      //algebra.u_star( );
      ///////////
      

      //algebra.w_equation();
      //algebra.solve_for_weights();
      //      volumes( T ); 

      
    }
    //    algebra.u_add_press_grad( dt2 );
//    draw( T , particle_file     );
//    draw_diagram( T , diagram_file );
//    return 0;

//    copy_weights( T ) ;
    
//    displ = move( T , dt , d0 );
//    volumes( T ); 
      
    //    algebra.fill_Delta_DD();
    //    algebra.p_equation( dt2 );

    cout
      << "Whole step  "
      << " : disp " << displ << endl ;

    //algebra.w_equation();
    //algebra.solve_for_weights();

    //    volumes( T ); 

    // half-step:
    //  update_full_vel( T );
    algebra.u_add_press_grad( dt );

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
