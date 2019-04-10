#include"pParticles.h"
#include"simu.h"

void draw(Triangulation& T,  const std::string file_name  ) {


  //  cout << "draw on step  "<< simu.current_step() << endl;

  std::stringstream  mkdir;
  std::stringstream  dirname;

  dirname << simu.current_step();
  mkdir << "mkdir -p "  << dirname.str();
  //  cout << "running : " << mkdir.str() << endl;
  system(mkdir.str().c_str());

    // std::stringstream  cp_cfg;
    // cp_cfg << "cp -f simu.cfg " << dirname.str();
    // system(cp_cfg.str().c_str());


  std::stringstream  namefile;
  namefile << simu.current_step() << '/' << file_name;

  //  cout << "writing on file : " << namefile.str() << endl;
  std::ofstream main_data;
  main_data.open(namefile.str().c_str() );

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    Point p = vit->point().point();

    #include"printout.h"

  }

  main_data << endl;

  main_data.close();
}



void draw_diagram(Triangulation& T,  const std::string file_name  ) {

  //  cout << "draw on step  "<< simu.current_step() << endl;

  std::stringstream  mkdir;
  std::stringstream  dirname;
    
  dirname << simu.current_step();
  mkdir << "mkdir -p "  << dirname.str();
  //cout << "running : " << mkdir.str() << endl;
  system(mkdir.str().c_str());

    // std::stringstream  cp_cfg;
    // cp_cfg << "cp -f simu.cfg " << dirname.str();
    // system(cp_cfg.str().c_str());


  std::stringstream  namefile;
  namefile << simu.current_step() << '/' << file_name;

  //cout << "writing on file : " << namefile.str() << endl;
  std::ofstream main_data;
  main_data.open(namefile.str().c_str() );

  for (  F_e_it eit = T.finite_edges_begin() ;
	 eit !=T.finite_edges_end(); ++eit) {

    CGAL::Object o = T.dual(eit);

    const Segment * Vor_segment = CGAL::object_cast<Segment>( &o );

    if (! Vor_segment ) continue;

    //    cout << " l0 = " << std::sqrt(Vor_segment->squared_length()) << endl;
    
    Point p1 = Vor_segment->source() ;
    Point p2 = Vor_segment->target() ;

    main_data << p1 << endl;
    main_data << p2 << endl;
    main_data << endl;

    
  }

  main_data << endl;

  main_data.close();
}



