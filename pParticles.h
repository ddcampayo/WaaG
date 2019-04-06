#ifndef _PPARTICLES_H_
#define _PPARTICLES_H_



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>

//#include <CGAL/Triangulation_vertex_base_2.h>

#include <CGAL/point_generators_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/centroid.h>


//#include <CGAL/spatial_sort.h>

using std::vector;
using std::list;
using std::endl;
using std::cout;
using std::cin;

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT FT;

// basic
typedef CGAL::Vector_2<K>          Vector_2;
typedef CGAL::Segment_2<K>         Segment;
typedef CGAL::Triangle_2<K>        Triangle;
typedef CGAL::Point_2<K>           Point;
typedef CGAL::Weighted_point_2<K>  wPoint;
typedef FT         weight;
typedef CGAL::Polygon_2<K> Polygon;

//typedef CGAL::Regular_triangulation_filtered_traits_2<K> TT;

//typedef TT::Weighted_point_2             weightedPoint;
//typedef TT::Weight                       weight;

// vrtx

#include"vrtx_info.h"

typedef My_vertex_base<K>   Vb;

// faces
typedef My_face_base<K>     Fb;

// Triangulation
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>         Tds;

typedef CGAL::Regular_triangulation_2< K , Tds>   Triangulation;

//typedef CGAL::Triangulation_data_structure_2< Vb >         Tds;
//typedef CGAL::Regular_triangulation_2< K >   Triangulation;


typedef Triangulation::Face_handle    Face_handle;
typedef Triangulation::Face           Face;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Vertex         Vertex;
typedef Triangulation::Edge           Edge;
//typedef Triangulation::Facet          Facet;
typedef Triangulation::Locate_type    Locate_type;
//typedef Triangulation::Point          Point;


typedef Triangulation::Finite_vertices_iterator   F_v_it;
typedef Triangulation::Finite_edges_iterator      F_e_it;
typedef Triangulation::Finite_faces_iterator      F_f_it;
typedef Triangulation::Edge_circulator            Edge_circulator;
typedef Triangulation::Vertex_circulator          Vertex_circulator;
typedef Triangulation::Face_circulator            Face_circulator;

void draw(Triangulation& T,  const std::string file_name  ) ;
void draw_diagram(Triangulation& T,  const std::string file_name  ) ;

void create(Triangulation& Tp, const FT& LL) ;
void volumes(Triangulation& T) ;
void number(Triangulation& T);
FT lloyds(Triangulation& T) ;

FT move(Triangulation& Tp, const FT dt , FT& dd0 ) ;
void backup( Triangulation& Tp );
void copy_weights( Triangulation& Tp );
void move_weights( Triangulation& Tp );
void update_full_vel(  Triangulation& Tp );
FT move_from_centroid(Triangulation& T, const FT dt );

void set_vels_rotating(Triangulation& T);
void set_vels_Lamb_Oseen(Triangulation& T) ;
void set_vels_Gresho(Triangulation& T) ;
FT L2_vel_Gresho( Triangulation& T) ;

#endif
