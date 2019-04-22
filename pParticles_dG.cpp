// pParticles
// Attempt to replicate de Goes  et al.
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

  //ysimu.do_perturb(0.01);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  linear algebra( T );


  // Init loop!
  
  const int max_iter = 100;
  const FT tol2 = 1e-3;
  int iter=0;

  for( ; iter < max_iter ; ++iter) {
  
    volumes( T ); 

    copy_weights( T ) ;

    algebra.solve_for_weights();

    FT dd = lloyds( T ) ;

    cout << " init loop , iter " << iter << " dd = " << dd << endl;
    if( dd < tol2) break;

  }

  cout << "Init loop converged in " << iter << " steps " << endl;
  
  set_vels_Gresho( T );

  volumes( T ); 

  FT d0;
  FT dt=0.001;

  cin >> dt ;

  simu.set_dt( dt );
  
//  move_from_centroid( T , dt);

  draw( T , particle_file);
  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");

  FT total_time = 1/( 2 * 3.14 * 0.2) ; // one whole turn

  do {
    simu.next_step();
    simu.advance_time( );

    backup( T );

    move_from_centroid( T , dt);
     
    algebra.solve_for_weights();

    copy_weights( T ) ;

    volumes( T ); 
  
    algebra.u_star( );

    // int iter = 1;

    // for ( ; iter < 40 ; iter++) {
    //   cout << "Iter  " << iter << "  ";
  
    //    FT displ = move( T , dt , d0 );
    //    cout << " : disp " << displ << endl ;
      
    algebra.p_equation( dt );

    algebra.u_add_press_grad( dt );

    //volumes( T ); 

    //   if( displ < 1e-8) break;


    // }

    volumes( T ); 
    
    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << endl ;

    
  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
