//modified by Elena Khramtcova 19.07.2012. The names of classes are not changed, even though it is 
//Linf metrics instead of L1, just for simplicity of building the demo.

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

#ifndef CGAL_LINF2D_VORONOI_DIAGRAM_TRAITS_2_H
#define CGAL_LINF2D_VORONOI_DIAGRAM_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/number_utils.h> 
#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Envelope_3/Env_plane_traits_3_functions.h>

namespace CGAL{

template <class Kernel_>
class Linf2D_voronoi_traits_2 : public Arr_linear_traits_2<Kernel_>
{
public:
  typedef Kernel_                              Kernel;
  typedef typename Kernel::FT                  FT;
  typedef Arr_linear_traits_2<Kernel>          Base;
  //L1_voronoi_traits
  typedef Linf2D_voronoi_traits_2<Kernel>          Self;
  typedef typename Base::Multiplicity          Multiplicity;

  typedef typename Base::Point_2               Point_2;
  typedef typename Base::Curve_2               Curve_2;
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
  
  typedef Point_2                  Xy_monotone_surface_3;
  typedef Point_2                  Surface_3;

  // Returns the distance between a two points in Linf metric.
  static FT distance(const Point_2& p1, const Point_2& p2) {
    FT dist = CGAL::max(CGAL::abs(p1.x() - p2.x()), CGAL::abs(p1.y() - p2.y()));
    return dist;
  }

  // Returns two enpoints of the middle (vertical/horizontal) part of the bisector, 
  //the index determines which point is returned. There are no conventions about how these endpoints situated 
  //depending on the index (like there was no in the oroginal traits) 
  static Point_2 mid_seg_endpoint(const Point_2& p1, const Point_2& p2, std::size_t index) {
    const Point_2 *pp1;
    const Point_2 *pp2;
    
    if (index % 2 == 0) {
      pp1 = &p1;
      pp2 = &p2;
    } else {
      pp1 = &p2;
      pp2 = &p1;
    }

    FT delta_x = pp2->x() - pp1->x();
    FT delta_y = pp2->y() - pp1->y();
    
    FT sign_x = CGAL::sign(delta_x);
    FT sign_y = CGAL::sign(delta_y);

    FT abs_x = CGAL::abs(delta_x);
    FT abs_y = CGAL::abs(delta_y);
    
    //midpoint between pp1 and pp2
    Point_2 mid = CGAL::midpoint(*pp1, *pp2);    
    
    //the case when bisector consists of only one line
    if(abs_x == abs_y){
        return mid;
   }

    FT mid_x = mid.x();
    FT mid_y = mid.y();

    //length of the middle part of bisector
    FT len_mid = CGAL::abs(abs_x - abs_y);
    
    CGAL_assertion(sign_x != CGAL::ZERO || sign_y !=CGAL::ZERO);

    //sign_x (resp. sign_y) is not zero, and will be different for different indexes 
    //since the order of input points changes 
    if (abs_x < abs_y){ 
      mid_x += sign_y * 0.5 * len_mid;
    }else{
      mid_y += sign_x * 0.5 * len_mid;
    }
    
    return Point_2(mid_x, mid_y);
  }

  static Comparison_result compare_z_at_xy (const X_monotone_curve_2& cv,
                                            const Xy_monotone_surface_3& h1,
                                            const Xy_monotone_surface_3& h2,
                                            bool above) {
   // std::cout << "Entered compare_z_at_xy(cv, h1, h2), where " << std::endl << "cv = " << cv << std::endl << "h1 = " << h1 << std::endl << "h2 = " << h2 << std::endl;
    CGAL::Comparison_result side = (above == true) ? CGAL::LARGER : CGAL::SMALLER;

    Line_2 l;
    if (cv.is_segment())
      l = cv.segment().supporting_line();
    else if (cv.is_ray())
      l = cv.ray().supporting_line();
    else
      l = cv.line();

      //std::cout << "l = " << l << std::endl;

    if ((l.is_vertical())){ 
      // To be "above" the curve, we actually need to have smaller x coordinate,
      // the order of the comparison function here is opposite to the none vertical
      // case.
      side = CGAL::opposite(side);
      CGAL::Comparison_result res = CGAL::compare_x_at_y(h1, l);
      if (res == side)
        return CGAL::SMALLER;
      else{
        CGAL_assertion(CGAL::compare_x_at_y(h2, l) == CGAL::opposite(res));
        return CGAL::LARGER;
      }
    } else {

      CGAL::Comparison_result res = CGAL::compare_y_at_x(h1, l);
      //std::cout << "res = " << res << "; side = " << side << std::endl;
      
      if (l.is_horizontal()) {
        CGAL_assertion(CGAL::compare_y_at_x(h2, l) != res);
        if (res == side)
         return CGAL::SMALLER;
        else{
          CGAL_assertion(CGAL::compare_y_at_x(h2, l) == CGAL::opposite(res));
          return CGAL::LARGER;
        }
      }
     
     //Could be a tie: points can be equidistant from some query points not on the line 
     //(which are parts of 2d bisector). In that case it returns CGAL::EQUAL

      if (res == side)
        return CGAL::SMALLER;
      res = CGAL::compare_y_at_x(h2, l);

      if(res == side)
        return CGAL::LARGER;

     // std::cout << "compare returns EQUAL: " << "h1 = " << h1 << "; h2 = " << h2 << std::endl; 
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
      return CGAL::compare(distance(p, h1), distance(p, h2));
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      Kernel k;
      Point_2 p;
      if(cv.is_segment())
        p = k.construct_midpoint_2_object()(cv.left(), cv.right());
      else
        if(cv.is_ray())
          p = k.construct_point_on_2_object()(cv.ray(), 1);
        else {
          CGAL_assertion(cv.is_line());
          p = k.construct_point_on_2_object()(cv.line(), 1);
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
      OutputIterator operator()(const Xy_monotone_surface_3& s1,
                                const Xy_monotone_surface_3& s2,
                                OutputIterator o) const {

      FT delta_x = s2.x() - s1.x();
      FT delta_y = s2.y() - s1.y();

      CGAL::Sign s_x = CGAL::sign(delta_x);
      CGAL::Sign s_y = CGAL::sign(delta_y);

      if (s_x == CGAL::ZERO && s_y == CGAL::ZERO)
        return o;

      if(CGAL::abs(delta_x) == CGAL::abs(delta_y)){
        //std::cout << "Points are in the corners of square" << std::endl;
        *o++ = CGAL::make_object(Intersection_curve(CGAL::bisector(s1, s2), 1));
        return o;
      }

      Point_2 p1 = mid_seg_endpoint(s1, s2, 0);
      Point_2 p2 = mid_seg_endpoint(s1, s2, 1);

      CGAL::Comparison_result cmppx = CGAL::compare(p2.x(), p1.x());
      CGAL::Comparison_result cmppy = CGAL::compare(p2.y(), p1.y());

      Point_2 *top, *bottom, *left, *right; 
    
      CGAL_assertion(cmppx == CGAL::EQUAL || cmppy == CGAL::EQUAL);
     
      if (cmppx == CGAL::LARGER) {
        left = &p1;
        right = &p2;
      }else if(cmppx == CGAL::SMALLER){
        left = &p2;
        right = &p1;
      }
      
      if (cmppy == CGAL::LARGER) {
        bottom = &p1;
        top = &p2;
      }else if(cmppy == CGAL::SMALLER){
        bottom = &p2;
        top = &p1;
      }
      
        //std::cout << "Assigned top, bottom, left, right" << std::endl;
        //std::cout << "s1 = " << s1 << "; s2 = " << s2 << std::endl;
        //std::cout << "p1 = " << p1 << "; p2 = " << p2 << std::endl;
      
      // We construct vertical/horyzontal line either way.
      *o++ = CGAL::make_object(Intersection_curve(Segment_2(p1, p2), 1));
      
      // Now construct rays with slope +-1. Two or four rays. If it is only two rays,
      // then the multiplicity of all the curves is 1.
      
      std::size_t mult;
      if(s_x == CGAL::ZERO || s_y == CGAL::ZERO){
        mult = 0;
      //  std::cout << "2d bisector" << std::endl; 
      }else{
        mult = 1;
      }
    
      bool topleft = false;
      bool topright = false;

      if (cmppx != CGAL::EQUAL) {
        
        // horizontal middle segment
        //if the leftmost point from {s1,s2} is above the horizontal line through p1 and p2, 
        //than slope of first and third piece of the bisector  1,
        //otherwise -1. If it is on the line, we add all four rays.
        
        if  (((s_x != CGAL::NEGATIVE) and (CGAL::compare(s2.y(), p1.y()) == CGAL::LARGER))|| ((s_x != CGAL::POSITIVE)and(CGAL::compare(s1.y(), p1.y()) == CGAL::LARGER))){
          topright = true;
        }

        if  (((s_x != CGAL::NEGATIVE) and (CGAL::compare(s1.y(), p1.y()) == CGAL::LARGER))||((s_x != CGAL::POSITIVE) and (CGAL::compare(s2.y(),p1.y()) == CGAL::LARGER))){
          topleft = true;
        }

       if(topleft){
          *o++ = CGAL::make_object(Intersection_curve(Ray_2(*right, Direction_2(1, 1)), mult));
          *o++ = CGAL::make_object(Intersection_curve(Ray_2(*left, Direction_2(-1, -1)), mult));
        
        }

        if(topright){
          *o++ = CGAL::make_object(Intersection_curve(Ray_2(*right, Direction_2(1, -1)), mult));          
          *o++ = CGAL::make_object(Intersection_curve(Ray_2(*left, Direction_2(-1, 1)), mult));
        }
      }
      
      if (cmppy != CGAL::EQUAL) {
                
        // vertical middle segment
        //the same strategy as in horizontal case
        if  (((s_y != CGAL::NEGATIVE) and (CGAL::compare(s2.x(), p1.x()) == CGAL::LARGER))|| ((s_y != CGAL::LARGER)and(CGAL::compare(s1.x(), p1.x()) == CGAL::LARGER))){
                topright = true;
        } 

        if  (((s_y != CGAL::NEGATIVE) and (CGAL::compare(s1.x(), p1.x()) == CGAL::LARGER))||((s_y != CGAL::POSITIVE) and (CGAL::compare(s2.x(),p1.x()) == CGAL::LARGER))){
                topleft = true;
        } 

                                        
        //if the topmost point from {s1,s2} is to the right of the vertical line through p1 and p2, 
        //than slope of first and third piece of the bisector  1,
        //otherwise -1. If it is on the line, we add all four rays.
        if(topright){
          *o++ = CGAL::make_object(Intersection_curve(Ray_2(*top, Direction_2(-1, 1)), mult));          
          *o++ = CGAL::make_object(Intersection_curve(Ray_2(*bottom, Direction_2(1,-1)), mult));
                //std::cout << "2 rays added because of topright" << std::endl;
        }
                                                                                                              
        if (topleft){
          *o++ = CGAL::make_object(Intersection_curve(Ray_2(*top, Direction_2(1, 1)), mult));         
          *o++ = CGAL::make_object(Intersection_curve(Ray_2(*bottom, Direction_2(-1, -1)), mult));
          //std::cout << "2 rays added because of topleft" << std::endl;
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
