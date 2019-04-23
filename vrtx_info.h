// simple class to define vector or scalar values

template <  class type >
class field {
 public:
  field() : f_()  {}
  type f() const {return f_;}
  type val() const {return f_;}
  type operator()() const {return f_;}

  void set(const type& ff) {f_ = ff;}
  void reset() {f_= type();}

  void operator =(const type& ff) {f_=ff;}
  type operator +(const type& ff) const {return f_ +ff;}
  type operator += (const type& ff) {type f2=  f_ + ff ; f_=f2 ; return f2;}
  type operator -= (const type& ff) {type f2=  f_ - ff ; f_=f2 ; return f2;}
  type operator/=(const FT& ff) {return (f_ = (1.0)/ff * f_ );}
  type operator*=(const FT& ff) {return (f_ =       ff * f_ );}

 private:
  type f_;
};


// vector that sets to zero
class autoVector_2
  : public Vector_2
{
 public:

 autoVector_2() : Vector_2(CGAL::NULL_VECTOR) {}
 autoVector_2(const Vector_2& v) : Vector_2(v) {}
};



#include"fields_enum.h"



/* From "A vertex class with an additionnal handle " CGAL example*/

template < class Gt, class Vb = CGAL::Regular_triangulation_vertex_base_2< Gt > >
class My_vertex_base  : public  Vb
{
  typedef Vb                              Base;
public:

  //  typedef typename Gt::Weight   Weight;
  //  typedef typename Gt::Point_2  Point;
  //  typedef typename Gt::Vector_2 Vector;


  //  typedef typename Vb::Triangulation_data_structure TDS;

  typedef typename Gt::Weighted_point_2 Weighted_point;
  //  typedef  Weighted_point      wPoint;
  typedef typename Vb::Face_handle   Face_handle;
  typedef typename Vb::Vertex_handle Vertex_handle;

  typedef typename Vb::Point              wPoint;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other    Vb2;
    typedef My_vertex_base< Gt ,Vb2>                         Other;
  };

public:
  My_vertex_base() : Base() {}
  My_vertex_base(const wPoint & p) : Base(p) {}
  My_vertex_base(const wPoint & p, Face_handle f) : Base(f,p) {}
  My_vertex_base(Face_handle f) : Base(f) {}

  typedef field<FT>            scalar_field;
  typedef field<autoVector_2>  vector_field;
  typedef field<Point>         point_field;
  typedef field<int>           int_field;
  //  typedef field<weight>        weight_field;

  scalar_field& sfield(const sfield_list::take sf ) {
    switch(sf)
      {
      case sfield_list::p : return p;
      case sfield_list::p0 : return p0;
      case sfield_list::vol0 : return vol0;
      case sfield_list::vol  : return vol;
      case sfield_list::w    : return w;
      case sfield_list::w0   : return w0;//  return this->point().weight();
      case sfield_list::s : return s;
      case sfield_list::I   : return I;
      case sfield_list::I0   : return I0;
      default : return vol;
      }
  }

  vector_field& vfield(const vfield_list::take vf ) {
    switch(vf)
      {
      case vfield_list::U  : return U;
      case vfield_list::U0  : return U0;
      case vfield_list::Ustar  : return Ustar;
      case vfield_list::Dr : return Dr;
      default : return U;
      }
  }

  
  
  ///////////

  // inter-cell data are stored as maps.-

  //   typedef std::map<Vertex_handle,      FT   >   scalar_link;
  //   typedef std::map<Vertex_handle,autoVector_2>  vector_link;

  
private:
//   scalar_link nabla_;   // vertices linked by discrete laplacian
//   scalar_link Delta_;   // vertices linked by discrete laplacian
//   vector_link grad_;    // idem grad (vector)

  Polygon poly_;

public:

  // vector_link& nabla(void) {
  //   return nabla_;
  // };

  // //  void get_nabla(vector_link& nn) {  nn=nabla_;  };

  // scalar_link& Delta() {
  //   return Delta_;
  // };

  Polygon poly() const {return poly_;}
  void set_poly(Polygon& pp)  { poly_ = pp ;}
  
  scalar_field p;          // pressure
  scalar_field p0;
  scalar_field vol;
  scalar_field vol0;

  vector_field U;
  vector_field U0;
  vector_field Ustar;

  point_field  r0;
  int_field    idx;
  vector_field Dr;
  // FT w() { return this->point().weight(); }
  scalar_field w0;
  scalar_field w;

  point_field  centroid;

  scalar_field s;
  scalar_field I;  // second moment of area
  scalar_field I0;

  //  weight_field w;

//  int idx() const {return idx_; }
//  void set_idx(int i)  { idx_ = i ; }
};


template <  class Gt , class Fb = CGAL::Regular_triangulation_face_base_2< Gt > >
  class My_face_base
  : public Fb
{
 public:

 typedef typename Fb::Vertex_handle                   Vertex_handle;
 typedef typename Fb::Face_handle                    Face_handle;

 template < typename TDS2 >
 struct Rebind_TDS {
   typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
   typedef My_face_base<Gt, Fb2>             Other;
 };

 public:
 My_face_base()
 : Fb() {}

 My_face_base(Vertex_handle v0,
			   Vertex_handle v1,
			   Vertex_handle v2)
 : Fb(v0,v1,v2) {}

 My_face_base(Vertex_handle v0,
			   Vertex_handle v1,
			   Vertex_handle v2,
			   Face_handle n0,
			   Face_handle n1,
			   Face_handle n2)
 : Fb(v0,v1,v2,n0,n1,n2) {}

 static int ccw(int i) {return CGAL::Triangulation_cw_ccw_2::ccw(i);}
 static int  cw(int i) {return CGAL::Triangulation_cw_ccw_2::cw(i);}

#ifndef CGAL_NO_DEPRECATED_CODE
 Vertex_handle mirror_vertex(int i) const;
 int mirror_index(int i) const;
#endif



};
