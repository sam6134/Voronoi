#include <iostream>
#include <list>
#include<string>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>
#include <CGAL/Segment_Delaunay_graph_site_2.h>
#include <CGAL/Polychain_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "Linf2D_voronoi_traits_2.h"
#include "L2_voronoi_traits_2.h"

#include <CGAL/envelope_3.h>
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

typedef CGAL::Exact_predicates_exact_constructions_kernel VD_Kernel;
typedef VD_Kernel::FT                                   Number_type;
typedef VD_Kernel::Iso_rectangle_2                      Iso_rectangle_2;
typedef VD_Kernel::Point_2                              VD_Point_2;
typedef VD_Kernel::Segment_2                            VD_Segment_2;
typedef std::vector<VD_Point_2>                         Points;

typedef CGAL::Linf2D_voronoi_traits_2<VD_Kernel>        VD_Traits_3;
typedef VD_Traits_3::Surface_3                          VD_Surface_3;
typedef CGAL::Envelope_diagram_2<VD_Traits_3>           VD_Envelope_diagram_2;

typedef CGAL::L2_voronoi_traits_2<VD_Kernel>            L2_VD_Traits_3;
typedef L2_VD_Traits_3::Surface_3                       L2_VD_Surface_3;
typedef CGAL::Envelope_diagram_2<L2_VD_Traits_3>        L2_VD_Envelope_diagram_2;

typedef Kernel::Point_2 Point_2;
typedef Kernel::Line_2 Line_2;

list<Point_2> pt_list;
list<VD_Point_2>  vd_pt_list;
Iso_rectangle_2 bbox;
Kernel::FT incr_len = 75;
void print_error_message(string s)
{
  cerr<<s<<endl;
  return;
}

int main()
{
  Point_2 points[5] = { Point_2(0,0), Point_2(10,0), Point_2(10,10), Point_2(6,5), Point_2(4,1) };
  double min_x = points[0].x();
  double min_y = points[0].y();
  double max_x = points[0].x();
  double max_y = points[0].y();
  // Calculate a bounding box of the Points ---------------
  for(int i=0;i<5;i++)
  {
    min_x = min(min_x,points[i].x());
    max_x = max(max_x,points[i].x());
    min_y = min(min_y,points[i].y());
    max_x = max(max_y,points[i].y());
  }
  CGAL::Point_2<VD_Kernel> UR(max_x,max_y), LL(min_x, min_y);
  bbox = Iso_rectangle_2(LL,UR);
  // increase the size of the Bounding box to handle extreme cases
  bbox = Iso_rectangle_2(
    bbox.min() + VD_Kernel::Vector_2(-incr_len, -incr_len),
    bbox.max() + VD_Kernel::Vector_2(incr_len, incr_len)
  );
  // ---------------------------------------------------------
  // Test inputs for L_inf Bisectors -------------------------
  pt_list.push_back(points[1]);
  pt_list.push_back(points[2]);
  int Option;
  cout<<"DEBUG: Select a Option"<<endl;
  cin>>Option;
  // L_Inf Bisector -----------------------------------------------------------------
  if(Option==1){
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
      cout<<"Ray emnating from "<<firstpt<<" in direction "<<incomingDir<<" and "<<outgoingDir<<endl;
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
  }
    // -----------------------------------------------------------------------------------------
    // Voronoi diagram for points
    // ------------------------------------------------------------------------------------------------
    else if (Option==2){
    {
      std::list<Point_2>::iterator it;
        for (it = pt_list.begin(); it != pt_list.end(); ++it) {
          vd_pt_list.push_back(VD_Point_2(it->x(), it->y()));
        }
        
        if (vd_pt_list.empty()) {
          print_error_message(("No mark selected"));
          return -1;
        }
    }

    {
      VD_Envelope_diagram_2 *m_envelope_diagram;
      m_envelope_diagram = new VD_Envelope_diagram_2();
    CGAL::lower_envelope_3(
      vd_pt_list.begin(),
      vd_pt_list.end(),
      *m_envelope_diagram
    );

    //computes the bounding box
    VD_Point_2 bottom_left(bbox.min().x(), bbox.min().y());
    VD_Point_2 top_right(bbox.max().x(), bbox.max().y());

    for (VD_Envelope_diagram_2::Vertex_const_iterator vit =
          m_envelope_diagram->vertices_begin();
         vit != m_envelope_diagram->vertices_end();
         vit++
    ) {
      VD_Point_2 vp = VD_Point_2(vit->point());
      if (CGAL::compare(vp.x(), bottom_left.x()) == CGAL::SMALLER)
        bottom_left = VD_Point_2(vp.x(), bottom_left.y());
      if (CGAL::compare(vp.y(), bottom_left.y()) == CGAL::SMALLER)
        bottom_left = VD_Point_2(bottom_left.x(), vp.y());
      if (CGAL::compare(vp.x(), top_right.x()) == CGAL::LARGER)
        top_right = VD_Point_2(vp.x(), top_right.y());
      if (CGAL::compare(vp.y(), top_right.y()) == CGAL::LARGER)
        top_right = VD_Point_2(top_right.x(), vp.y());       
    }

    CGAL::Point_2<VD_Kernel> bl(to_double(bottom_left.x()), to_double(bottom_left.y()));
    CGAL::Point_2<VD_Kernel> tr(to_double(top_right.x()), to_double(top_right.y()));

    Kernel::FT incr_len= 50;

    bbox = Iso_rectangle_2(
      bl + VD_Kernel::Vector_2(-incr_len,-incr_len),
      tr + VD_Kernel::Vector_2(incr_len,incr_len)
    );

    // The edges of voronoi diagram
    cout<<"The Voronoi diagram is"<<endl;
    VD_Envelope_diagram_2::Edge_const_iterator eit;
    for (eit = m_envelope_diagram->edges_begin();
    eit != m_envelope_diagram->edges_end(); eit++) {
      if (eit->curve().is_segment()) {
        // Case when VD is only a segment
        Point_2 p1(
          to_double(eit->curve().segment().source().x()),
          to_double(eit->curve().segment().source().y())
        );
        Point_2 p2(
          to_double(eit->curve().segment().target().x()),
          to_double(eit->curve().segment().target().y())
        );
        cout<<"A segment "<<p1<<" and "<<p2<<endl;
      } else if (eit->curve().is_ray()) {
        Point_2 p(
          to_double(eit->curve().ray().source().x()),
          to_double(eit->curve().ray().source().y())
        );
        CGAL::Direction_2<Kernel> d(
          to_double(eit->curve().ray().direction().dx()),
          to_double(eit->curve().ray().direction().dy())
        );
        cout<<"A Ray emanating from "<<p<<" in direction "<<d<<endl;
      } else if (eit->curve().is_line()) {
        Line_2 l(
          to_double(eit->curve().line().a()),
          to_double(eit->curve().line().b()),
          to_double(eit->curve().line().c())
        );
        cout<<"A Line of form ax+by+c with a,b,c "<<l<<endl;
      }
    }
  }
 }
  return 0;
}