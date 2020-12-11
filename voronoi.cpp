#include <iostream>
#include <list>
#include <string>
#include<queue>
#include<deque>

#include <CGAL/intersections.h>
// envelope
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>
//
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


typedef Kernel::Point_2                                 Point_2;
typedef Kernel::Line_2                                  Line_2;
typedef Kernel::Segment_2                               Segment_2;
typedef Kernel::Intersect_2                             Intersect_2;
typedef Kernel::Line_2                                  Line_2;
typedef Kernel::Ray_2                                   Ray_2;
typedef Kernel::Direction_2                             Direction_2;
typedef Kernel::Vector_2                                Vector_2;

list<Point_2> pt_list;
list<VD_Point_2>  vd_pt_list;

struct HullSegment
{
  Segment_2 e1;
  Ray_2 b;
  Segment_2 e2;
};

bool PointComp(Point_2 a, Point_2 b)
{
  return a.x()>b.x();
}
void print_error_message(string s)
{
  cerr<<s<<endl;
  return;
}
void print_ray(Point_2 p, CGAL::Direction_2<Kernel> d){
  cout<<"R "<<p<<" "<<d<<" \n";
}
// Plot the segment p1---p2
void print_segment(Point_2 p1 , Point_2 p2){
  cout<<"S "<<p1 <<" "<< p2 << " \n";
}
// Plot the point p
void print_point(Point_2 p){
  cout<<"P "<<p << " \n";
}
// Plot the line a.x + b.y + c = 0
void print_line(Line_2 l){
  printf("L %f %f %f \n",l.a(),l.b(),l.c());
};

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
    if(a == 1 || a == 3)
    {
      if((s->source().y() - s->target().y())/(s->source().x() - s->target().x())> 0) return false;
    }
    if(a == 2 || a == 4)
    {
      if((s->source().y() - s->target().y())/(s->source().x() - s->target().x())< 0) return false;
    }
    if(a == 1 || a == 2)
    {
      Line_2 lcheck = Line_2(*s);
      if(((lcheck.a()*x + lcheck.b()*y + lcheck.c())/lcheck.b()) > 0) return false;
    }

    if(a == 3 || a == 4)
    {
      Line_2 lcheck = Line_2(*s);
      if(((lcheck.a()*x + lcheck.b()*y + lcheck.c())/lcheck.b()) < 0) return false;
    }

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
  Point_2 points[6] = { Point_2(9,4), Point_2(10,5), Point_2(7,7), Point_2(4,6), Point_2(2,1), Point_2(9,7) };
  


  // ---------------------------------------------------------
  // Test inputs for L_inf Bisectors and VD -------------------------
  pt_list.push_back(points[0]);
  pt_list.push_back(points[1]);
  pt_list.push_back(points[2]);
  pt_list.push_back(points[3]);
  pt_list.push_back(points[4]);
  pt_list.push_back(points[5]);

  //for( auto x: pt_list){
    //print_point(x);
  //}

  int Option;
  cerr<<"DEBUG: Select a Option"<<endl;
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
    cerr<<"The Voronoi diagram is"<<endl;
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
        cerr<<"A segment "<<p1<<" and "<<p2<<endl;
        print_segment(p1,p2);
      } else if (eit->curve().is_ray()) {
        Point_2 p(
          to_double(eit->curve().ray().source().x()),
          to_double(eit->curve().ray().source().y())
        );
        CGAL::Direction_2<Kernel> d(
          to_double(eit->curve().ray().direction().dx()),
          to_double(eit->curve().ray().direction().dy())
        );
        cerr<<"A Ray emanating from "<<p<<" in direction "<<d<<endl;
        print_ray(p,d);
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
        cerr<<"Segment whose endpoints are "<<p1<<" and "<<p2<<endl;
        print_segment(p1,p2);
      } else if (eit->curve().is_ray()) {
        Point_2 p(
          to_double(eit->curve().ray().source().x()),
          to_double(eit->curve().ray().source().y())
        );
        CGAL::Direction_2<Kernel> d(
          to_double(eit->curve().ray().direction().dx()),
          to_double(eit->curve().ray().direction().dy())
        );
        cerr<<"Ray emnating from "<<p<<" in direction "<<d<<endl;
        print_ray(p,d);
      } else if (eit->curve().is_line()) {
        Line_2 l(
          to_double(eit->curve().line().a()),
          to_double(eit->curve().line().b()),
          to_double(eit->curve().line().c())
        );
        cerr<<"Line of form ax+by+c with a,b,c as "<<l<<endl;
        print_line(l);
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
    cerr<<"The given segments are "<<endl;
    list<Segment_2>::iterator it;
    for(it=seg_list.begin();it != seg_list.end(); it++)
    {
      cerr<<"Segment with end points "<< it->source() <<" , "<<it->target()<<endl;
      print_segment(it->source(),it->target());
    }
    // Find Region R

    // find ln y coord
    double ln = max(seg_list.begin()->source().y(), seg_list.begin()->target().y());
    // find ls y coord
    double ls = min(seg_list.begin()->source().y(), seg_list.begin()->target().y());
    // find le x coord
    double le = max(seg_list.begin()->source().x(), seg_list.begin()->target().x());
    // find lw x coord
    double lw = min(seg_list.begin()->source().x(), seg_list.begin()->target().x());
    for(it=seg_list.begin();it != seg_list.end(); it++)
    {
      ln = min(ln, max(it->source().y(), it->target().y()));
      ls = max(ls, min(it->source().y(), it->target().y()));
      le = min(le, max(it->source().x(), it->target().x()));
      lw = max(lw, min(it->source().x(), it->target().x()));
    }
    cerr<<endl;
    cerr<<"The Region R is "<<endl;
    cerr<<"ls = "<<ls<<" ln= "<<ln<<endl;
    cerr<<"lw = "<<lw<<" le= "<<le<<endl;

    // dummy segment to store modified segment
    Segment_2* sMod = new Segment_2(Point_2(1,1), Point_2(1,1));


    list<Segment_2> Quadrant1;
    list<Segment_2> Quadrant2;
    list<Segment_2> Quadrant3;
    list<Segment_2> Quadrant4;


    for(it = seg_list.begin();it!=seg_list.end();it++)
    {

      // find segments straddling quadrant-1
      if(intersectsRegion(&(*it),ls,lw,1,sMod))
      {
        Quadrant1.push_back(*sMod);
      }

      // find segments straddling quadrant-2
      if(intersectsRegion(&(*it),ls,le,2,sMod))
      {
        Quadrant2.push_back(*sMod);
      }

      // find segments straddling quadrant-3
      if(intersectsRegion(&(*it),ln,le,3,sMod))
      {
        Quadrant3.push_back(*sMod);
      }

      // find segments straddling quadrant-4
      if(intersectsRegion(&(*it),ln,lw,4,sMod))
      {
        Quadrant4.push_back(*sMod);
      }
    }

    // print the quad-1 segments
    cerr<<"Quadrant-1"<<endl;
    for(it = Quadrant1.begin(); it != Quadrant1.end(); it++)
    {
      cerr<<"Segment with end points "<< it->source() << " , "<<it->target()<<endl;
      //print_segment(it->source(), it->target());
    }
    // print the quad-2 segments
    cerr<<"Quadrant-2"<<endl;
    for(it = Quadrant2.begin(); it != Quadrant2.end(); it++)
    {
      cerr<<"Segment with end points "<< it->source() << " , "<<it->target()<<endl;
      //print_segment(it->source(), it->target());
    }
    // print the quad-3 segments
    cerr<<"Quadrant-3"<<endl;
    for(it = Quadrant3.begin(); it != Quadrant3.end(); it++)
    {
      cerr<<"Segment with end points "<< it->source() << " , "<<it->target()<<endl;
      //print_segment(it->source(), it->target());

    }
    // print the quad-4 segments
    cerr<<"Quadrant-4"<<endl;
    for(it = Quadrant4.begin(); it != Quadrant4.end(); it++)
    {
      cerr<<"Segment with end points "<< it->source() << " , "<<it->target()<<endl;
      //print_segment(it->source(), it->target());
    }

    /// Determining the Envelope
    {
    typedef CGAL::Cartesian<double>                         Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
    typedef Segment_traits_2::X_monotone_curve_2            Segment_2_monotone;
    typedef CGAL::Arr_curve_data_traits_2<Segment_traits_2,
                                          char>             Traits_2;
    typedef Traits_2::Point_2                               Point_2_monotone;
    typedef Traits_2::X_monotone_curve_2                    Labeled_segment_2;
    typedef CGAL::Envelope_diagram_1<Traits_2>              Diagram_1;

      std::list<Labeled_segment_2>   montone_list;
      char dummy = 'A';
      list<Point_2> Envelope1;
      for(it = Quadrant1.begin();it!=Quadrant1.end();it++)
      {
        if(it->source() == it->target()) 
        {
          Envelope1.push_back(it->source());
          continue;
        }
        montone_list.push_back(Labeled_segment_2 (Segment_2_monotone \
        (Point_2_monotone(it->source().x(), it->source().y()),Point_2_monotone(it->target().x(), it->target().y())),dummy));
        dummy++;
      }
      Diagram_1 min_diag;
      cerr<<"chk1"<<endl;
      upper_envelope_x_monotone_2 (montone_list.begin(), montone_list.end(),
                               min_diag);
      cerr<<"upper env done"<<endl; 
      Diagram_1::Edge_const_handle     e = min_diag.leftmost();
      Diagram_1::Vertex_const_handle   v;
      Diagram_1::Curve_const_iterator  cit;
      while (e != min_diag.rightmost())
      {
        cerr << "Edge:";
        if (! e->is_empty())
        {
          for (cit = e->curves_begin(); cit != e->curves_end(); ++cit)
            cerr << ' ' << cit->data();
        }
        else
          cerr << " [empty]";
        cerr << std::endl;
        v = e->right();
        cerr << "Vertex (" << v->point() << "):";
        Envelope1.push_back(v->point());
        for (cit = v->curves_begin(); cit != v->curves_end(); ++cit)
          cerr << ' ' << cit->data();
        cerr << std::endl;
        e = v->right();
      }
      if(Envelope1.empty()) Envelope1.push_back(Point_2(lw,ls));


      // quad 2
      montone_list.clear();
      dummy = 'A';
      list<Point_2> Envelope2;
      for(it = Quadrant2.begin();it!=Quadrant2.end();it++)
      {
        if(it->source() == it->target()) 
        {
          Envelope2.push_back(it->source());
          continue;
        }
        montone_list.push_back(Labeled_segment_2 (Segment_2_monotone \
        (Point_2_monotone(it->source().x(), it->source().y()),Point_2_monotone(it->target().x(), it->target().y())),dummy));
        dummy++;
      }

      Diagram_1 min_diag2;
      cerr<<"chk1"<<endl;
      upper_envelope_x_monotone_2 (montone_list.begin(), montone_list.end(),
                               min_diag2);
      cerr<<"upper env done"<<endl; 
      e = min_diag2.leftmost();
      
      while (e != min_diag2.rightmost())
      {
        cerr << "Edge:";
        if (! e->is_empty())
        {
          for (cit = e->curves_begin(); cit != e->curves_end(); ++cit)
            cerr << ' ' << cit->data();
        }
        else
          cerr << " [empty]";
        cerr << std::endl;
        v = e->right();
        cerr << "Vertex (" << v->point() << "):";
        Envelope2.push_back(v->point());
        for (cit = v->curves_begin(); cit != v->curves_end(); ++cit)
          cerr << ' ' << cit->data();
        cerr << std::endl;
        e = v->right();
      }
      if(Envelope2.empty()) Envelope2.push_back(Point_2(le,ls));

      // quad 3
      montone_list.clear();
      dummy = 'A';
      list<Point_2> Envelope3;
      for(it = Quadrant3.begin();it!=Quadrant3.end();it++)
      {
        if(it->source() == it->target()) 
        {
          Envelope3.push_back(it->source());
          continue;
        }
        montone_list.push_back(Labeled_segment_2 (Segment_2_monotone \
        (Point_2_monotone(it->source().x(), it->source().y()),Point_2_monotone(it->target().x(), it->target().y())),dummy));
        dummy++;
      }
      Diagram_1 min_diag3;
      cerr<<"chk1"<<endl;
      lower_envelope_x_monotone_2 (montone_list.begin(), montone_list.end(),
                               min_diag3);
      cerr<<"upper env done"<<endl; 
      e = min_diag3.leftmost();
      
      while (e != min_diag3.rightmost())
      {
        cerr << "Edge:";
        if (! e->is_empty())
        {
          for (cit = e->curves_begin(); cit != e->curves_end(); ++cit)
            cerr << ' ' << cit->data();
        }
        else
          cerr << " [empty]";
        cerr << std::endl;
        v = e->right();
        cerr << "Vertex (" << v->point() << "):";
        Envelope3.push_back(v->point());
        for (cit = v->curves_begin(); cit != v->curves_end(); ++cit)
          cerr << ' ' << cit->data();
        cerr << std::endl;
        e = v->right();
      }
      if(Envelope3.empty()) Envelope3.push_back(Point_2(le,ln));



      // quad4
      montone_list.clear();
      dummy = 'A';
      list<Point_2> Envelope4;
      for(it = Quadrant4.begin();it!=Quadrant4.end();it++)
      {
        if(it->source() == it->target()) 
        {
          Envelope4.push_back(it->source());
          continue;
        }
        montone_list.push_back(Labeled_segment_2 (Segment_2_monotone \
        (Point_2_monotone(it->source().x(), it->source().y()),Point_2_monotone(it->target().x(), it->target().y())),dummy));
        dummy++;
      }
      Diagram_1 min_diag4;
      cerr<<"chk1"<<endl;
      lower_envelope_x_monotone_2 (montone_list.begin(), montone_list.end(),
                               min_diag4);
      cerr<<"upper env done"<<endl; 
      e = min_diag4.leftmost();
      
      
      while (e != min_diag4.rightmost())
      {
        cerr << "Edge:";
        if (! e->is_empty())
        {
          for (cit = e->curves_begin(); cit != e->curves_end(); ++cit)
            cerr << ' ' << cit->data();
        }
        else
          cerr << " [empty]";
        cerr << std::endl;
        v = e->right();
        cerr << "Vertex (" << v->point() << "):";
        Envelope4.push_back(v->point());
        for (cit = v->curves_begin(); cit != v->curves_end(); ++cit)
          cerr << ' ' << cit->data();
        cerr << std::endl;
        e = v->right();
      }
      if(Envelope4.empty()) Envelope4.push_back(Point_2(lw,ln));


      // circular queue implementation
    {
      list<HullSegment> FarthestHull;
      list<Point_2>::iterator it;
      for(it=Envelope1.begin();it!=Envelope1.end();it++)
      {
        Point_2 curr_point = *it;
        HullSegment h1;
        h1.b = Ray_2(curr_point, Direction_2(Vector_2(1,-1)));
        if(it == Envelope1.begin())
        {
          h1.e1 = Segment_2(*Envelope4.begin(), curr_point);
        }else{
          h1.e1 = Segment_2(*(std::prev(it)),*it);
        }
        if(*it == *Envelope1.rbegin())
        {
          h1.e2 = Segment_2(curr_point, *Envelope2.begin());
        }
        else{
          h1.e2 = Segment_2(*it,*(std::next(it)));
        }
        FarthestHull.push_back(h1);
      }
      // push envelope 2
      for(it=Envelope2.begin();it!=Envelope2.end();it++)
      {
        Point_2 curr_point = *it;
        HullSegment h1;
        h1.b = Ray_2(curr_point, Direction_2(Vector_2(-1,-1)));
        if(it == Envelope2.begin())
        {
          h1.e1 = Segment_2(*Envelope1.begin(), curr_point);
        }else{
          h1.e1 = Segment_2(*(std::prev(it)),*it);
        }
        if(*it == *Envelope2.rbegin())
        {
          h1.e2 = Segment_2(curr_point, *Envelope3.rbegin());
        }
        else{
          h1.e2 = Segment_2(*it,*(std::next(it)));
        }
        FarthestHull.push_back(h1);
      }

      // push envelope 3
      reverse(Envelope3.begin(), Envelope3.end());
      for(it=Envelope3.begin();it!=Envelope3.end();it++)
      {
        Point_2 curr_point = *it;
        HullSegment h1;
        h1.b = Ray_2(curr_point, Direction_2(Vector_2(1,1)));
        if(it == Envelope3.begin())
        {
          h1.e1 = Segment_2(*Envelope2.rbegin(), curr_point);
        }else{
          h1.e1 = Segment_2(*(std::prev(it)),*it);
        }
        if(*it == *Envelope3.rbegin())
        {
          h1.e2 = Segment_2(curr_point, *Envelope4.rbegin());
        }
        else{
          h1.e2 = Segment_2(*it,*(std::next(it)));
        }
        FarthestHull.push_back(h1);
      }


      // push envelope 4
      reverse(Envelope4.begin(), Envelope4.end());
      for(it=Envelope4.begin();it!=Envelope4.end();it++)
      {
        Point_2 curr_point = *it;
        HullSegment h1;
        h1.b = Ray_2(curr_point, Direction_2(Vector_2(-1,1)));
        if(it == Envelope4.begin())
        {
          h1.e1 = Segment_2(*Envelope3.rbegin(), curr_point);
        }else{
          h1.e1 = Segment_2(*(std::prev(it)),*it);
        }
        if(*it == *Envelope4.rbegin())
        {
          h1.e2 = Segment_2(curr_point, *Envelope1.begin());
        }
        else{
          h1.e2 = Segment_2(*it,*(std::next(it)));
        }
        FarthestHull.push_back(h1);
      }
      
      list<HullSegment>::iterator it1;
      cerr<<"The farthest Hull is"<<endl;

      print_line(Line_2(0,1,-ln));
      print_line(Line_2(0,1,-ls));
      print_line(Line_2(1,0,-lw));
      print_line(Line_2(1,0,-le));

      // if(Envelope1.size()>1){
      // Point_2 prev_point = *Envelope1.begin();
      // for(it=Envelope1.begin();it!=Envelope1.end();it++)
      // {
      //   Point_2 curr_point = *it;
      //   print_segment(prev_point, curr_point);
      //   prev_point = curr_point;
      // }
      // }

      // if(Envelope2.size()>1){
      // Point_2 prev_point = *Envelope2.begin();
      // for(it=Envelope2.begin();it!=Envelope2.end();it++)
      // {
      //   Point_2 curr_point = *it;
      //   print_segment(prev_point, curr_point);
      //   prev_point = curr_point;
      // }
      // }



      // if(Envelope3.size()>1){
      // Point_2 prev_point = *Envelope3.begin();
      // for(it=Envelope3.begin();it!=Envelope3.end();it++)
      // {
      //   Point_2 curr_point = *it;
      //   print_segment(prev_point, curr_point);
      //   prev_point = curr_point;
      // }
      // }


      // if(Envelope4.size()>1){
      // Point_2 prev_point = *Envelope4.begin();
      // for(it=Envelope4.begin();it!=Envelope4.end();it++)
      // {
      //   Point_2 curr_point = *it;
      //   print_segment(prev_point, curr_point);
      //   prev_point = curr_point;
      // }
      // }

      
      for(it1 = FarthestHull.begin(); it1 != FarthestHull.end(); it1 ++)
      {
        print_point(it1->b.source());
        print_segment(it1->e1.source(), it1->e1.target());
        print_segment(it1->e2.source(), it1->e2.target());
      }
      
      
      
      
      vector<HullSegment> circularList;
      

      

    }

    }





    
    
  }

  return 0;  
}

  
