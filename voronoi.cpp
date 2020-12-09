#include <iostream>
#include <list>
#include<string>
#include <CGAL/intersections.h>
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

#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include "CGAL/L2_segment_voronoi_traits_2.h"

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
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Intersect_2 Intersect_2;
typedef Kernel::Line_2 Line_2;

list<Point_2> pt_list;
list<VD_Point_2>  vd_pt_list;

void print_error_message(string s)
{
  cerr<<s<<endl;
  return;
}
// to check whether a point lies in a quadrants as described in the paper
bool PointInRegion(Point_2 p, double y, double x,int a)
{
  // for quad 1
  if(a == 1)
  {
    if(p.x()>=x && p.y()>=y)
    {
      return true;
    }else return false;
  }
  // for quad 2
  if(a == 2)
  {
    if(p.x()<=x && p.y()>=y)
    {
      return true;
    }else return false;
  }
  // for quad 3
  if(a == 3)
  {
    if(p.x()<=x && p.y()<=y)
    {
      return true;
    }else return false;
  }
  // for quad 4
  if(a == 4)
  {
    if(p.x()>=x && p.y()<=y)
    {
      return true;
    }else return false;
  }
  return false;
}

// function to check whether a segment lies partially inside a region
bool intersectsRegion(Segment_2* s, double y, double x, int a, Segment_2* sMod)
{

    // find intx points with both the lines
    Line_2 lin(0,1,-y);
    CGAL::cpp11::result_of<Intersect_2(Segment_2, Line_2)>::type
    result1 = CGAL::intersection(*s, lin);

    Line_2 lin1(1,0,-x);
    CGAL::cpp11::result_of<Intersect_2(Segment_2, Line_2)>::type
    result2 = CGAL::intersection(*s, lin1);
    // both intersections
    if(result1 && result2)
    {
      const Point_2* p1 = boost::get<Point_2 >(&*result1);
      const Point_2* p2 = boost::get<Point_2 >(&*result2);
      (*sMod) = Segment_2(*p1, *p2);
      return true;
    }
    // else if both intxs are not there 

  return false;
}

bool isInRegion(Segment_2* s, double y, double x, int a)
{
  // for quad 1
  if(a == 1)
  {
    if(s->source().x()>=x && s->target().x()>=x && s->source().y()>=y && s->target().y()>=y)
    {
      return true;
    }else return false;
  }

  // for quad 2
  if(a == 2)
  {
    if(s->source().x()<=x && s->target().x()<=x && s->source().y()>=y && s->target().y()>=y)
    {
      return true;
    }else return false;
  }

  // for quad 3
  if(a == 3)
  {
    if(s->source().x()<=x && s->target().x()<=x && s->source().y()<=y && s->target().y()<=y)
    {
      return true;
    }else return false;
  }

  // for quad 4
  if(a == 4)
  {
    if(s->source().x()>=x && s->target().x()>=x && s->source().y()<=y && s->target().y()<=y)
    {
      return true;
    }else return false;
  }
 return false;
}

int main()
{
  Point_2 points[6] = { Point_2(0,0), Point_2(10,4), Point_2(10,10), Point_2(6,5), Point_2(4,1), Point_2(3,7) };
  


  // ---------------------------------------------------------
  // Test inputs for L_inf Bisectors and VD -------------------------
  pt_list.push_back(points[0]);
  pt_list.push_back(points[1]);
  pt_list.push_back(points[2]);
  pt_list.push_back(points[3]);
  pt_list.push_back(points[4]);
  pt_list.push_back(points[5]);
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
    // Voronoi diagram for points in Linf Space
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
 //Voronoi Diagram in L2 space for points
 else if(Option == 3)
 {
        std::list<Point_2>::iterator it;
        for (it = pt_list.begin(); it != pt_list.end(); ++it) {
          vd_pt_list.push_back(VD_Point_2(it->x(), it->y()));
        }
        
        if (vd_pt_list.empty()) {
          print_error_message(("No mark selected"));
          return -1;
        }
      L2_VD_Envelope_diagram_2 *m_envelope_diagram;
      m_envelope_diagram = new L2_VD_Envelope_diagram_2();
    CGAL::lower_envelope_3(
      vd_pt_list.begin(),
      vd_pt_list.end(),
      *m_envelope_diagram
    );
    
    // print the edges
    for (L2_VD_Envelope_diagram_2::Edge_const_iterator eit =
          m_envelope_diagram->edges_begin();
         eit != m_envelope_diagram->edges_end();
         eit++
    ) {
      if (eit->curve().is_segment()) {
        Point_2 p1(
          to_double(eit->curve().segment().source().x()),
          to_double(eit->curve().segment().source().y())
        );
        Point_2 p2(
          to_double(eit->curve().segment().target().x()),
          to_double(eit->curve().segment().target().y())
        );
        cout<<"Segment whose endpoints are "<<p1<<" and "<<p2<<endl;
      } else if (eit->curve().is_ray()) {
        Point_2 p(
          to_double(eit->curve().ray().source().x()),
          to_double(eit->curve().ray().source().y())
        );
        CGAL::Direction_2<Kernel> d(
          to_double(eit->curve().ray().direction().dx()),
          to_double(eit->curve().ray().direction().dy())
        );
        cout<<"Ray emnating from "<<p<<" in direction "<<d<<endl;
      } else if (eit->curve().is_line()) {
        Line_2 l(
          to_double(eit->curve().line().a()),
          to_double(eit->curve().line().b()),
          to_double(eit->curve().line().c())
        );
        cout<<"Line of form ax+by+c with a,b,c as "<<l<<endl;
      }
    } 
  }
  else if(Option == 4)
  {
    list<Segment_2> seg_list;
    // test segments
    for(int i=0;i<6;i+=2)
    {
      seg_list.push_back(Segment_2(points[i],points[i+1]));
    }
    cout<<"The given segments are "<<endl;
    list<Segment_2>::iterator it;
    for(it=seg_list.begin();it != seg_list.end(); it++)
    {
      cout<<"Segment with end points "<< it->source() <<" , "<<it->target()<<endl;
    }
    // Find Region R

    // find ln y coord
    double ln = max(seg_list.begin()->source().y(), seg_list.begin()->target().y());
    for(it=seg_list.begin();it != seg_list.end(); it++)
    {
      ln = min(ln, max(it->source().y(), it->target().y()));
    }


    // find ls y coord
    double ls = min(seg_list.begin()->source().y(), seg_list.begin()->target().y());
    for(it=seg_list.begin();it != seg_list.end(); it++)
    {
      ln = max(ln, min(it->source().y(), it->target().y()));
    }

    // find le x coord
    double le = max(seg_list.begin()->source().x(), seg_list.begin()->target().x());
    for(it=seg_list.begin();it != seg_list.end(); it++)
    {
      le = min(le, max(it->source().x(), it->target().x()));
    }

    // find lw x coord
    double lw = min(seg_list.begin()->source().x(), seg_list.begin()->target().x());
    for(it=seg_list.begin();it != seg_list.end(); it++)
    {
      lw = max(lw, min(it->source().x(), it->target().x()));
    }
    cout<<endl;
    cout<<"The Region R is "<<endl;
    cout<<"ls = "<<ls<<" ln= "<<ln<<endl;
    cout<<"lw = "<<lw<<" le= "<<le<<endl;

    // dummy segment to store modified segment
    Segment_2* sMod = new Segment_2(Point_2(1,1), Point_2(1,1));


    // find segments straddling quadrant-1
    list<Segment_2> Quadrant1;
    for(it = seg_list.begin();it!=seg_list.end();it++)
    {
      if(intersectsRegion(&(*it),ls,lw,1,sMod))
      {
        Quadrant1.push_back(*sMod);
      }
    }

    // print the quad-1 segments
     cout<<"Quadrant-1"<<endl;
    for(it = Quadrant1.begin(); it != Quadrant1.end(); it++)
    {
      cout<<"Segment with end points "<< it->source() << " , "<<it->target()<<endl;
    }


    // find segments straddling quadrant-2
    list<Segment_2> Quadrant2;
    for(it = seg_list.begin();it!=seg_list.end();it++)
    {
      if(intersectsRegion(&(*it),ls,le,2,sMod))
      {
        Quadrant2.push_back(*sMod);
      }
    }

    // print the quad-2 segments
     cout<<"Quadrant-2"<<endl;
    for(it = Quadrant2.begin(); it != Quadrant2.end(); it++)
    {
      cout<<"Segment with end points "<< it->source() << " , "<<it->target()<<endl;
    }

    // find segments straddling quadrant-3
    list<Segment_2> Quadrant3;
    for(it = seg_list.begin();it!=seg_list.end();it++)
    {
      if(intersectsRegion(&(*it),ln,le,3,sMod))
      {
        Quadrant3.push_back(*sMod);
      }
    }

    // print the quad-3 segments
     cout<<"Quadrant-3"<<endl;
    for(it = Quadrant3.begin(); it != Quadrant3.end(); it++)
    {
      cout<<"Segment with end points "<< it->source() << " , "<<it->target()<<endl;
    }

    // find segments straddling quadrant-4
    list<Segment_2> Quadrant4;
    for(it = seg_list.begin();it!=seg_list.end();it++)
    {
      if(intersectsRegion(&(*it),ln,lw,4,sMod))
      {
        Quadrant4.push_back(*sMod);
      }
    }

    // print the quad-4 segments
     cout<<"Quadrant-4"<<endl;
    for(it = Quadrant4.begin(); it != Quadrant4.end(); it++)
    {
      cout<<"Segment with end points "<< it->source() << " , "<<it->target()<<endl;
    }
  }
  return 0;  
 }
  
