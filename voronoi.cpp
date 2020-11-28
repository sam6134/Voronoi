#include <iostream>
#include <list>
#include<string>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>
#include <CGAL/Segment_Delaunay_graph_site_2.h>
#include <CGAL/Polychain_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>

using namespace std;
typedef CGAL::Cartesian <double> Kernel;
struct Gt_inf
: public CGAL::Segment_Delaunay_graph_Linf_traits_2<Kernel,CGAL::Field_with_sqrt_tag> {};
typedef Gt_inf::Site_2 Site_2;
typedef
  CGAL::SegmentDelaunayGraphLinf_2::Bisector_Linf<Gt_inf> Inf_bis;
Inf_bis bisector_linf;
typedef Kernel::Point_2 Point_2;
list<Point_2> pt_list;

void print_error_message(string s)
{
  cerr<<s<<endl;
  return;
}

int main()
{
  Point_2 points[5] = { Point_2(0,0), Point_2(10,0), Point_2(10,10), Point_2(6,5), Point_2(4,1) };

  pt_list.push_back(points[3]);
  pt_list.push_back(points[4]);
  // -----------------------------------------------------------------
  if (pt_list.empty()) {
        print_error_message(("No points are there"));
        return -1;
      }
      if (pt_list.size() != 2) {
        print_error_message(("Exactly two points should be there for L-inf bisector"));
        return -1;
      }
  //-- now calculate the l-inf bisector----------------------------------
    std::list<Point_2>::iterator it;
    it = pt_list.begin();
    Point_2 p = *it;
    Site_2 sp = Site_2::construct_site_2(p);
    ++it;
    Point_2 q = *it;
    Site_2 sq = Site_2::construct_site_2(q);
    Inf_bis::Polychainline pcl = bisector_linf(sp, sq) ;
    Inf_bis::Polychainline::Vertex_const_iterator it1 = pcl.vertices_begin();

    if (pcl.size() == 1) {
      Point_2 firstpt = *it1;
      CGAL::Direction_2<Kernel> incomingDir = pcl.get_incoming();
      CGAL::Direction_2<Kernel> outgoingDir = pcl.get_outgoing();
      cout<<"Ray emnating from "<<firstpt<<endl;
    }

    else if (pcl.size() == 2) {
      Point_2 firstpt = *it1;
      ++it1;
      Point_2 lastpt = *it1;
      cout<<"Segment whose end points are "<<firstpt<<" and "<<lastpt<<endl;

    }

    else {
      Point_2 firstpt = *it1;
      cout<<"Ray emnating from "<<firstpt<<endl;
      Inf_bis::Polychainline::Vertex_const_iterator it2 = it1+1;
      for (; it2!=pcl.vertices_end(); ++it1, ++it2) {
        Point_2 x1 = *it1;
        Point_2 x2 = *it2;
        cout<<"Segment whose end points are "<<x1<<"and "<<x2<<endl;
      }
    }
  return 0;
}