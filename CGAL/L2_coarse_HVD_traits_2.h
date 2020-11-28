//modified by Elena Khramtcova 7.08.2012

// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Envelope_3/include/CGAL/Env_plane_traits_3.h $
// $Id: Env_plane_traits_3.h 51989 2009-09-21 10:55:53Z efif $
//
// Author(s)     : Ophir Setter

#ifndef CGAL_L2_HVD_TRAITS_2_H
#define CGAL_L2_HVD_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Envelope_3/Env_plane_traits_3_functions.h>
#include <CGAL/L2_voronoi_traits_2.h>
#include <vector>

#include <CGAL/envelope_3.h>
#include <CGAL/Polygon_2.h>

namespace CGAL{

template <class Kernel_>
class L2_HVD_traits_2 : public Arr_linear_traits_2<Kernel_>
{
public:
  typedef Kernel_                              Kernel;
  typedef typename Kernel::FT                  FT;
  typedef Arr_linear_traits_2<Kernel>          Base;
  typedef L2_HVD_traits_2<Kernel>          Self;
  typedef typename Base::Multiplicity          Multiplicity;

  typedef typename Base::Point_2               Point_2;
  typedef typename std::vector<Point_2>        Points;
  typedef typename Base::Curve_2               Curve_2;

  typedef typename CGAL::Polygon_2<Kernel>                  Cluster;
  typedef typename CGAL::Polygon_2<Kernel>::Vertex_iterator Cluster_Iterator;

  typedef typename Base::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Kernel::Segment_2           Segment_2;
  typedef typename Kernel::Ray_2               Ray_2;
  typedef typename Kernel::Line_2              Line_2;
  typedef typename Kernel::Direction_2         Direction_2;
  typedef std::pair<Curve_2, Multiplicity>     Intersection_curve;

  typedef typename Base::Left_side_category    Left_side_category;
  typedef typename Base::Bottom_side_category  Bottom_side_category;
  typedef typename Base::Top_side_category     Top_side_category;
  typedef typename Base::Right_side_category   Right_side_category;


  typedef typename CGAL::L2_voronoi_traits_2<Kernel> FVD_Traits_3;
  typedef typename  FVD_Traits_3::Surface_3   FVD_Surface_3;
  typedef typename CGAL::Envelope_diagram_2<FVD_Traits_3>
    FVD_Envelope_diagram_2;
  typedef typename FVD_Envelope_diagram_2::Edge_const_iterator
    Edge_const_iterator;
  typedef typename FVD_Envelope_diagram_2::Surface_const_iterator
    Surface_const_iterator;


  typedef Cluster                  Xy_monotone_surface_3;
  typedef Cluster                  Surface_3;

  typedef typename std::pair<FT, Point_2> Distance_to_Cluster;

  static void farthest_point(
      const Point_2& p, const Cluster& c, Point_2& result) {
    CGAL_assertion(c.size() > 0);
    Cluster_Iterator cit = c.vertices_begin();
    result = *cit;
    FT max_sqr_dist = CGAL::squared_distance(*cit, p);
    for(; cit != c.vertices_end(); cit++) {
      FT cur_sqr_dist  = CGAL::squared_distance(*cit, p);
      if (CGAL::compare(max_sqr_dist, cur_sqr_dist) == CGAL:: SMALLER){
        max_sqr_dist = cur_sqr_dist;
        result = *cit;
      }
    }
    return;
  }




  static FT squared_farthest_distance(const Point_2& p, const Cluster& c) {
    Point_2 res;
    farthest_point(p, c, res);
    return (CGAL::squared_distance(p, res));
  }


static void nearest_point(
          const Point_2& p, const Cluster& c, Point_2& result) {
      CGAL_assertion(c.size() > 0);
          Cluster_Iterator cit = c.vertices_begin();
              result = *cit;
                  FT min_sqr_dist = CGAL::squared_distance(*cit, p);
                      for(; cit != c.vertices_end(); cit++) {
                              FT cur_sqr_dist  = CGAL::squared_distance(*cit, p);
                                    if (CGAL::compare(min_sqr_dist, cur_sqr_dist) == CGAL::LARGER) {
                                              min_sqr_dist = cur_sqr_dist;
                                                      result = *cit;
                                                            }
                                        }
                               return;
                                }

  static FT squared_nearest_distance(const Point_2& p, const Cluster& c) {
        Point_2 res;
            nearest_point(p, c, res);
                return (CGAL::squared_distance(p, res));
                  }






  static Comparison_result compare_z_at_xy (const X_monotone_curve_2& cv,
                                            const Xy_monotone_surface_3& h1,
                                            const Xy_monotone_surface_3& h2,
                                            bool above) {
    CGAL::Comparison_result side = (above == true) ? CGAL::LARGER : CGAL::SMALLER;
    Line_2 l;
    Point_2 p;
    Kernel k;
    if (cv.is_segment()){
      l = cv.segment().supporting_line();
      p = k.construct_midpoint_2_object()(cv.left(), cv.right());
    }else if (cv.is_ray()){
      l = cv.ray().supporting_line();
      p = k.construct_point_on_2_object()(cv.ray(), 1);
    }else{
      l = cv.line();
      p = k.construct_point_on_2_object()(cv.line(), 1);
    }
    Point_2* p_at_h1;
    Point_2* p_at_h2;

    p_at_h1 = new Point_2();
    p_at_h2 = new Point_2();

    farthest_point(p, h1, *p_at_h1);
    farthest_point(p, h2, *p_at_h2);


    if ((l.is_vertical())) {
      // To be "above" the curve, we actually need to have smaller x coordinate,
      // the order of the comparison function here is opposite to the none vertical
      // case.
      side = CGAL::opposite(side);
        CGAL::Comparison_result res = CGAL::compare_x_at_y(*p_at_h1, l);
       if (res == side)
        return CGAL::SMALLER;
      else{
        CGAL_assertion(CGAL::compare_x_at_y(*p_at_h2, l) == CGAL::opposite(res));
        return CGAL::LARGER;
      }
    } else {

      CGAL::Comparison_result res = CGAL::compare_y_at_x(*p_at_h1, l);

      if (l.is_horizontal()) {
        CGAL_assertion(CGAL::compare_y_at_x(*p_at_h2, l) != res);
        if (res == side)
         return CGAL::SMALLER;
        else{
          CGAL_assertion(CGAL::compare_y_at_x(*p_at_h2, l) == CGAL::opposite(res));
          return CGAL::LARGER;
        }
      }

     //Could be a tie: points can be equidistant from some query points not on the line
     //(which are parts of 2d bisector). In that case it returns CGAL::EQUAL

      if (res == side)
        return CGAL::SMALLER;
      res = CGAL::compare_y_at_x(*p_at_h2, l);

      if(res == side)
        return CGAL::LARGER;

      return CGAL::EQUAL;
    }
  }

  class Make_xy_monotone_3 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Surface_3& s,
                                bool /* is_lower */,
                                OutputIterator o) const {
      *o++ = s;
      return o;
    }
  };

  Make_xy_monotone_3 make_xy_monotone_3_object() const {
    return Make_xy_monotone_3();
  }

  class Compare_z_at_xy_3 {
  public:
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return CGAL::compare(squared_farthest_distance(p, h1),
                           squared_farthest_distance(p, h2));
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
          Kernel k;
      Point_2 p;
      if(cv.is_segment()) {
        p = k.construct_midpoint_2_object()(cv.left(), cv.right());
      }
      else {
        if(cv.is_ray()) {
          p = k.construct_point_on_2_object()(cv.ray(), 1);
        }
        else {
          CGAL_assertion(cv.is_line());
          p = k.construct_point_on_2_object()(cv.line(), 1);
        }
      }
      return this->operator()(p, h1, h2);
    }

    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      // should happen only if the points are equal.
      CGAL_assertion(h1 == h2);
      return EQUAL;
    }
  };

  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  {
    return Compare_z_at_xy_3();
  }

  class Compare_z_at_xy_above_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return compare_z_at_xy (cv, h1, h2, true);
    }

  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3();
  }

  class Compare_z_at_xy_below_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return compare_z_at_xy (cv, h1, h2, false);
    }
  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const {
    return Compare_z_at_xy_below_3();
  }

  class Construct_projected_boundary_2 {
  public:
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const {
      return o;
    }
  };

  Construct_projected_boundary_2
  construct_projected_boundary_2_object() const
  {
    return Construct_projected_boundary_2();
  }

  class Construct_projected_intersections_2 {
  public:
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& c1,
                              const Xy_monotone_surface_3& c2,
                              OutputIterator o) const
    {
      Points pts = c1.container();
      pts.insert( pts.end(), c2.vertices_begin(), c2.vertices_end() );
      FVD_Envelope_diagram_2 *FVD_12;
      FVD_12 = new FVD_Envelope_diagram_2();
      CGAL::upper_envelope_3 (pts.begin(), pts.end(), *FVD_12);
      for(Edge_const_iterator
          eit = FVD_12->edges_begin();
          eit != FVD_12->edges_end();
          eit++) {
        Surface_const_iterator sit = (*eit).surfaces_begin();
        CGAL_assertion((*eit).number_of_surfaces() == 2);
        Point_2 p1 = *sit;
        sit++;
        Point_2 p2 = *sit;
        bool p1_in_c1 =
          (std::find(c1.vertices_begin(), c1.vertices_end(), p1)
            != c1.vertices_end());
        bool p2_in_c1 =
          (std::find(c1.vertices_begin(), c1.vertices_end(), p2)
            != c1.vertices_end());
        bool p1_in_c2 =
          (std::find(c2.vertices_begin(), c2.vertices_end(), p1)
            != c2.vertices_end());
        bool p2_in_c2 =
          (std::find(c2.vertices_begin(), c2.vertices_end(), p2)
            != c2.vertices_end());
        if (((p1_in_c1)&&(p2_in_c2))||((p1_in_c2)&&(p2_in_c1))){
          *o++ = CGAL::make_object(Intersection_curve((*eit).curve(), 1));
        } else {
          CGAL_assertion(((p1_in_c1)&&(p2_in_c1))||((p1_in_c2)&&(p2_in_c2)));
        }
      }

      return o;
    }
  };

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const {
    return Construct_projected_intersections_2();
  }

};

} //namespace CGAL

#endif // CGAL_LINF2D_VORONOI_DIAGRAM_TRAITS_2_H
