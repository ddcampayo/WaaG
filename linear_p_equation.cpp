//#define PRESSURE_PPE_DIV_SOURCE
// TODO: better define different PPE functions, independent

//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"

// Solve for pressure
// The famous pressure Poisson equation

void linear::p_equation(const FT dt , const bool ws ) {

  // choose!!

  //   p_equation_divgrad_div_source(dt,ws);
  p_equation_lapl_div_source(dt);
  // p_equation_lapl_Dvol_source( dt );

  return;
}

//  Laplacian = - Delta / (2 V), div v as source

void linear::p_equation_lapl_div_source(const FT dt ){

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  
  VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );
  VectorXd p =  Delta_solver.solve( divUstar );
  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  vctr_to_field( -0.5 * p / ddt ,  sfield_list::p ) ;

  return;
}


// same as above, but div u calculated as in FEM

void linear::p_equation_lapl_div_source_fem(const FT dt ){

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  
  VectorXd divUstar  =  DD_scalar_vfield_fem( vfield_list::Ustar );
  VectorXd p =  Delta_solver.solve( divUstar );
  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  vctr_to_field( -0.5 * p / ddt ,  sfield_list::p ) ;

  return;
}




//  Laplacian = - Delta / (2 V), Delta vol as source

void linear::p_equation_lapl_Dvol_source(const FT dt ){

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  VectorXd vol0  = field_to_vctr( sfield_list::vol0 ) ;

  VectorXd Dvol = vol.array() - vol0.array()  ;

  VectorXd Dp =  Delta_solver.solve( Dvol );

  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  vctr_to_field( -0.5 * Dp / ( ddt * ddt) , sfield_list::p  ) ;

  
  return;
}

//  Laplacian = - Delta / (2 V), Delta Delanay vol as source

void linear::p_equation_lapl_Dvol_source_fem(const FT dt ){

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd vol  = field_to_vctr( sfield_list::Dvol ) ;

  VectorXd vol0  = field_to_vctr( sfield_list::Dvol0 ) ;

  VectorXd Dvol = vol.array() - vol0.array()  ;

  VectorXd Dp =  Delta_solver.solve( Dvol );

  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  vctr_to_field( -0.5 * Dp / ( ddt * ddt) , sfield_list::p  ) ;

  
  return;
}



// //  Laplacian = - Delta / (2 V), Delta FEM vol as source

// void linear::p_equation_lapl_Dvol_source_fem(const FT dt ){

//   cout << "Solving pressure equation " << endl;
  
//   //  fill_Delta_DD(); // This may be important -- or not

//   FT ddt = dt;
//   if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

//   VectorXd vol  = field_to_vctr( sfield_list::Dvol ) ;

//   VectorXd vol0  = field_to_vctr( sfield_list::Dvol0 ) ;

//   VectorXd Dvol = vol.array() - vol0.array()  ;

//   VectorXd Dp =  Delta_solver.solve( Dvol );

//   // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
//   vctr_to_field( -0.5 * Dp / ( ddt * ddt) , sfield_list::p  ) ;

  
//   return;
// }



//  Laplacian = div of grad, div v as source

void linear::p_equation_divgrad_div_source(const FT dt , const bool ws ) {

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );
  VectorXd p;

  // Possible correction due to w  field .-
  if( ws ) {
    VectorXd w  = field_to_vctr( sfield_list::w );
    VectorXd w0 = field_to_vctr( sfield_list::w0 );

    VectorXd Delta_w = Delta * ( w - w0);
    p =  LL_solver.solve( divUstar + Delta_w / ddt );
  }
  else
    p =  LL_solver.solve( divUstar  );
  
  vctr_to_field( p / ddt ,  sfield_list::p ) ;

  return;
}


//  Laplacian = div of grad, D vol as source

void linear::p_equation_divgrad_Dvol_source(const FT dt , const bool ws ) {

  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  
  // diagnostics on volumes.-

  VectorXd vol  = field_to_vctr( sfield_list::vol ) ;

  VectorXd vol0  = field_to_vctr( sfield_list::vol0 ) ;

  VectorXd Dvol = vol.array() - vol0.array()  ;

  int N = vol.size();

  FT Dvol_sigma =  Dvol.array().square().sum() / N ; // / FT( vol.size() );
  FT Dvol_mean  =  vol.array().sum() / N ; // / FT( vol.size() );

  cout << "Pressure  "
       << " rel Dvol std dev: " << sqrt( Dvol_sigma ) / Dvol_mean 
       << endl;

  // C1: LL Laplacian
  VectorXd Dp  =  LL_solver.solve( Dvol );

  vctr_to_field( Dp / ( ddt * ddt) , sfield_list::p  ) ;


  return;
}



void linear::p_equation_s(const FT dt ) {

  cout << "Solving pressure equation with s term" << endl;
  
  FT ddt = dt;
  if( dt < 1e-10 ) ddt = 1;  // for debugging, mainly

  VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );

  VectorXd s  = field_to_vctr( sfield_list::s );

  VectorXd LNs = LN * s;

  VectorXd p =  Delta_solver.solve( divUstar / ddt + LNs );
  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V

  vctr_to_field( -0.5 * p ,  sfield_list::p ) ;

  
  //VectorXd  p =  LL_solver.solve( divUstar / dt + LNs );
  
  //  vctr_to_field( p ,  sfield_list::p ) ;

  return;
}

