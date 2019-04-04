// pParticles

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
  volumes( T );
  copy_weights( T ) ;

  //  set_vels_rotating( T );
//  set_vels_Lamb_Oseen( T );
  set_vels_Gresho( T );

  draw( T , particle_file);
  draw_diagram( T , diagram_file );

  FT d0;
  FT dt2=0.001;

  cin >> dt2 ;

  simu.set_dt( dt2 );

  simu.next_step();
  simu.advance_time( );

  backup( T );

  volumes( T ); 

  //  volumes( T );

  linear algebra( T );

  //  algebra.fill_Delta();
  algebra.solve_for_weights();

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );

  simu.next_step();
  simu.advance_time();

  algebra.u_star( );

  FT  displ = move( T , dt2 , d0 );

  algebra.w_equation();
  //algebra.solve_for_weights();
  //volumes( T ); 

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );

  simu.next_step();
  simu.advance_time();

  algebra.p_equation( dt2 );

  algebra.u_add_press_grad( dt2 );

  volumes( T ); 

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );

  //  algebra.u_star( );

  for (int iter = 0 ; iter < 40 ; iter++) {
    cout << "Iter  " << iter << "  ";
    simu.next_step();
    simu.advance_time();
  
    displ = move( T , dt2 , d0 );
    cout << " : disp " << displ << endl ;

    if( displ < 1e-8) break;
    
    algebra.w_equation();
    // algebra.solve_for_weights();
    //    volumes( T ); 

    algebra.p_equation( dt2 );

    algebra.u_add_press_grad( dt2 );

    volumes( T ); 

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );
  }
  
  std::ofstream log_file;

  log_file.open("main.log");

  log_file
    << simu.current_step() << "  "
    << simu.time() << "  " << endl ;

  log_file.close();


  return 0;

}
