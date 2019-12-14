// Example program for the centroid() function for 2D points, 3D points and 3D triangles.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/centroid.h>
#include <vector>
#include <iostream>
typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  K;
typedef K::Point_2                  Point_2;

int main()
{
  // centroid of 2D points
  std::vector<Point_2> points_2;

  points_2.push_back(Point_2( -0.32500000000000007,                 -0.375 ));
  points_2.push_back(Point_2(-0.32500000000000007,              -0.42499999999999999 ));
  points_2.push_back(Point_2(-0.27500000000000002, -0.42499999999999999 ));
  points_2.push_back(Point_2(  -0.27500000000000008, -0.375 ));
  //  points_2.push_back(Point_2(  -0.32500000000000007, -0.375 ));

  
  Point_2 c2 = CGAL::centroid(points_2.begin(), points_2.end(),CGAL::Dimension_tag<0>());
  std::cout << c2 << std::endl;


  return 0;
}
