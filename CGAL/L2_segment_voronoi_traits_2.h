// created by Aron Fiechter on 2017-07-13.
// This file implements a model of the concept EnvelopeTraits_3.
// It is mostly is a copy of L2_voronoi_traits_2, except that it is for segments
// instead of points.

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

#ifndef CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H
#define CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/number_utils.h>
#include <CGAL/Envelope_3/Envelope_base.h>

/* to convert from Alg to Rat and viceversa */
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Cartesian.h>

/* to compute bisector of two segments */
#include <CGAL/ch_akl_toussaint.h> // convex hull

/* to format strings */
#include <boost/format.hpp>

namespace CGAL {

/* colour */
#define COUT_COLOUR_RESET   "\033[0m"
#define COUT_COLOUR_RED     "\033[31m"      /* Red */
#define COUT_COLOUR_GREEN   "\033[32m"      /* Green */
#define COUT_COLOUR_YELLOW  "\033[33m"      /* Yellow */
#define COUT_COLOUR_BLUE    "\033[34m"      /* Blue */

template <class Conic_traits_2, class Kernel_>
class L2_segment_voronoi_traits_2 : public Conic_traits_2 {

public:
  typedef Kernel_                                     Kernel;
  typedef Conic_traits_2                              C_traits_2;
  typedef L2_segment_voronoi_traits_2<C_traits_2, Kernel> Self;

  typedef typename C_traits_2::Point_2                Point_2;
  typedef typename C_traits_2::Curve_2                Curve_2;
  typedef typename C_traits_2::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename C_traits_2::Multiplicity           Multiplicity;

  typedef typename C_traits_2::Rat_kernel             Rat_kernel;
  typedef typename C_traits_2::Alg_kernel             Alg_kernel;
  typedef typename C_traits_2::Nt_traits              Nt_traits;

  typedef typename Nt_traits::Integer                 BigInt;

  typedef typename Rat_kernel::FT                     Rational;
  typedef typename Rat_kernel::Point_2                Rat_point_2;
  typedef typename Rat_kernel::Segment_2              Rat_segment_2;
  typedef typename Rat_kernel::Line_2                 Rat_line_2;
  typedef typename Rat_kernel::Ray_2                  Rat_ray_2;
  typedef typename Rat_kernel::Vector_2               Rat_vector_2;
  typedef typename Rat_kernel::Direction_2            Rat_direction_2;
  typedef typename CGAL::Polygon_2<Rat_kernel>        Rat_polygon_2;
  typedef typename Rat_polygon_2::Edge_const_iterator Edge_iterator;

  typedef typename Alg_kernel::FT                     Algebraic;
  typedef typename Alg_kernel::Point_2                Alg_point_2;
  typedef typename Alg_kernel::Segment_2              Alg_segment_2;
  typedef typename Alg_kernel::Line_2                 Alg_line_2;
  typedef typename Alg_kernel::Ray_2                  Alg_ray_2;
  typedef typename Alg_kernel::Direction_2            Alg_direction_2;

  typedef Rational                                    RT;
  typedef Algebraic                                   AT;

  /* Converters */
  typedef CGAL::Cartesian_converter<Alg_kernel, Rat_kernel> AK_to_RK;
  typedef CGAL::Cartesian_converter<Rat_kernel, Alg_kernel> RK_to_AK;
  typedef CGAL::Cartesian<double>                     D_kernel;
  typedef typename D_kernel::Point_2                  D_point_2;
  typedef CGAL::Cartesian_converter<Alg_kernel, D_kernel> AK_to_DK;
  typedef CGAL::Cartesian_converter<D_kernel, Rat_kernel> DK_to_RK;

  /* Define segments as surfaces for simplicity. Each segment should be thought
   * of as a surface representing the function distance, with higher values on
   * the z axis mean farthest points */
  typedef Rat_segment_2                               Surface_3;
  typedef Surface_3                                   Xy_monotone_surface_3;

protected:
  typedef std::pair<X_monotone_curve_2, Multiplicity> Intersection_curve;

private:
  typedef typename std::pair<
    std::pair<Rat_line_2, Rat_line_2>,
    std::pair<Rat_line_2, Rat_line_2>
  >                                                   Rat_delimiter_lines;
  typedef typename std::pair<
    std::pair<Alg_line_2, Alg_line_2>,
    std::pair<Alg_line_2, Alg_line_2>
  >                                                   Alg_delimiter_lines;

  /* Enum to store the three type of bisector parts that can be found in the
   * bisector of two segments (excluding the unbounded rays) */
  enum Bisector_type {
    PARABOLIC_ARC,
    SUPP_LINE_BISECTOR,
    ENDPOINT_BISECTOR
  };

  /* Enum to store the four possible quadrants something can be oriented to.
   * Used when rotating a segment. */
  enum Quadrant {
    FIRST_Q,
    SECOND_Q,
    THIRD_Q,
    FOURTH_Q
  };



/* ########################################################################## */
/* ###                           CLASS PARABOLA                           ### */
/* ########################################################################## */

  /* Parabola class used to provide methods to find intersections with lines,
   * to check point location with respect to a parabola, to create parabolic
   * arcs, to get tangent directions at points. */
  class Parabola {

  private:
    /* Fields */

    /* coefficients of equation: rx^2 + sy^2 + txy + ux + vy + w = 0 */
    RT _r; RT _s; RT _t; RT _u; RT _v; RT _w;

    /* generators of parabola, direction of parabola is the same of directrix;
     * also store orientation for when we build arcs */
    Rat_line_2 _directrix;
    Rat_point_2 _focus;
    Orientation _orientation;

  public:

    /* Empty constructor */
    Parabola() {}

    /* Construct using directrix and focus. Details on the computation of the
     * forumla can be found in doc/parabola.pdf */
    Parabola(Rat_line_2 directrix, Rat_point_2 focus)
      : _directrix(directrix), _focus(focus) {

      /* orientation depends on directrix an focus */
      if (directrix.has_on_positive_side(focus)) {
        this->_orientation = CGAL::COUNTERCLOCKWISE;
      }
      else {
        this->_orientation = CGAL::CLOCKWISE;
      }

      /* extract parameters to compute coefficients */
      RT a = directrix.a();
      RT b = directrix.b();
      RT c = directrix.c();
      RT f_x = focus.x();
      RT f_y = focus.y();
      RT TWO = RT(2);

      /* compute coefficients, see details in doc/parabola.pdf */
      RT r = -CGAL::square(b);
      RT s = -CGAL::square(a);
      RT t = TWO * a * b;
      RT u =
        TWO * a * c +
        TWO * CGAL::square(a) * f_x +
        TWO * CGAL::square(b) * f_x
      ;
      RT v =
        TWO * b * c +
        TWO * CGAL::square(a) * f_y +
        TWO * CGAL::square(b) * f_y
      ;
      RT w =
        CGAL::square(c) -
        CGAL::square(a) * CGAL::square(f_x) -
        CGAL::square(a) * CGAL::square(f_y) -
        CGAL::square(b) * CGAL::square(f_x) -
        CGAL::square(b) * CGAL::square(f_y)
      ;

      this->_r = r;
      this->_s = s;
      this->_t = t;
      this->_u = u;
      this->_v = v;
      this->_w = w;

      RK_to_AK to_alg;
      CGAL_assertion(CGAL::square(_t) - 4 * _r * _s == 0);  // curve is parabola
      CGAL_assertion(this->has_on(to_alg(
        CGAL::midpoint(focus, directrix.projection(focus))  // origin is on p
      )));
    }

    /* Getters */
    RT r() { return _r; }
    RT s() { return _s; }
    RT t() { return _t; }
    RT u() { return _u; }
    RT v() { return _v; }
    RT w() { return _w; }
    Rat_line_2 directrix() { return _directrix; }
    Rat_point_2 focus() { return _focus; }
    Orientation orientation() { return _orientation; }

    /* Methods */

    /* Evaluate the equation of the parabola rx^2 + sy^2 + txy + ux + vy + w = 0
     * using the x and y of the point. If the result is 0, the point is on the
     * parabola, if it is positive the point lies on the positive side, if it
     * is negative the point lies on the negative side. */
    Algebraic evaluate(Point_2 point) {
      Algebraic x = point.x();
      Algebraic y = point.y();
      Algebraic result(
        this->r() * CGAL::square(x) +
        this->s() * CGAL::square(y) +
        this->t() * x * y +
        this->u() * x +
        this->v() * y +
        this->w()
      );
      return result;
    }

    /* Check if a given point lies on the parabola by checking if the values of
     * x and y (the point's coordinates) satisfy the equation of the parabola */
    bool has_on(Point_2 point) {
      return CGAL::is_zero(this->evaluate(point));
    }
    /* Check if a given point lies on the positive side of the parabola. The
     * positive side is the one on the left when traveling on the curve in the
     * same direction as the directrix (by construction they have the "same"
     * oriented side) */
    bool has_on_positive_side(Point_2 point) {
      return CGAL::is_positive(this->evaluate(point));
    }
    /* Check if a given point lies on the negative side of the parabola. The
     * negative side is the one on the right when traveling on the curve in the
     * same direction as the directrix (by construction they have the "same"
     * oriented side) */
    bool has_on_negative_side(Point_2 point) {
      return CGAL::is_negative(this->evaluate(point));
    }

    /* Given a point, return the tangent line at that point. The tangent can be
     * found by looking at the line that is perpendicular to the directrix and
     * goes through this point and the line that goes through the focus and this
     * point. The bisector of these two lines is the tangent at this point.
     * Precondition (checked): the point is on the parabola. */
    Alg_line_2 tangent_at_point(Alg_point_2 point) {
      /* check precondition */
      CGAL_precondition_msg(
        this->has_on(point),
        "Given point for tangent has to be on the parabola"
      );

      /* get converter and convert */
      RK_to_AK to_alg;
      Alg_line_2 alg_directrix = to_alg(this->directrix());
      Alg_point_2 focus = to_alg(this->focus());
      Alg_point_2 proj_point = alg_directrix.projection(point);

      /* tangent */
      Alg_line_2 tangent;

      /* based on an orientation test, distinguish three cases */
      switch (CGAL::orientation(proj_point, point, focus)) {

        /* special case: the point is the vertex of the parabola. In this case,
         * the tangent at point has the same direction as the directrix */
        case CGAL::COLLINEAR: {
          tangent = Alg_line_2(point, alg_directrix.direction());
          break;
        }

        /* in this case the tangent is directed towards the positive part of the
         * directrix, so we need to have one line directed from the focus to the
         * point and the other one from the directrix to the point */
        case CGAL::LEFT_TURN: {
          tangent = CGAL::bisector(
            Alg_line_2(focus, point),
            Alg_line_2(proj_point, point)
          );
          break;
        }

        /* in this case the tangent is directed towards the negative part of the
         * directrix, so we need to have one line directed from the point to the
         * focus and the other one from the point to the directrix */
        case CGAL::RIGHT_TURN: {
          tangent = CGAL::bisector(
            Alg_line_2(point, focus),
            Alg_line_2(point, proj_point)
          );
          break;
        }

        /* impossible, throw error */
        default: {
          CGAL_error_msg(
            "Point on parabola, its projection on the directrix and the focus"
            " failed to have one of three orientations."
          );
          break;
        }
      }

      /* invert tangent if necessary */
      if (!generally_same_direction(tangent, alg_directrix.direction())) {
        tangent = tangent.opposite();
      }
      return tangent;
    }

    /* Save into the OutputIterator o the intersection(s) of the parabola with
     * a given line l. The type of o must be Alg_point_2.
     * Return a past the end iterator o. */
    template <class OutputIterator>
    OutputIterator get_intersections(Rat_line_2 line, OutputIterator o) {
      /* equation of line:      ax + by + c = 0
       * equation of parabola:  rx^2 + sy^2 + txy + ux + vy + w = 0
       * we can find intersections by substituting line in parabola; described
       * in detail in docs/parabola.pdf, verified using wxMaxima */
      RT a = line.a();
      RT b = line.b();
      RT c = line.c();

      /* convert line to algebraic, get nt_traits to solve quadratic equation */
      RK_to_AK to_alg;
      Alg_line_2 alg_line = to_alg(line);
      Nt_traits nt_traits;

      /* in this case the intersection is simpler, since we can substitute the x
       * in the parabola equation with just:    x = -c/a
       * We get a quadratic equation in y, which we can solve using CGAL */
      if (b == 0) {
        /* the quadratic equation in y is:
         *                                    s y^2
         *                     + (v - (ct / a)) y
         *      + ((rc^2 / a^2) - (cu / a) + w)
         *                                          = 0
         */
        RT EQ_A = this->s();
        RT EQ_B = this->v() - ((c * this->t()) / a);
        RT EQ_C = ((this->r() * CGAL::square(c)) / CGAL::square(a)) -
          ((c * this->u()) / a) +
          this->w()
        ;

        /* to store the 0, 1, or 2 results, indicating the intersections.
         * For all resulting y, find the corresponding x, and add a point to
         * the OutputIterator o */
        Algebraic  ys[2];
        Algebraic * ys_end;
        int n_ys;
        ys_end = nt_traits.solve_quadratic_equation(EQ_A, EQ_B, EQ_C, ys);
        n_ys = ys_end - ys;

        /* if no intersections return */
        if (n_ys == 0) {
          return o;
        }
        /* else find xs for all ys, add points to iterator */
        else while (--n_ys >= 0) {
          Algebraic current_y = ys[n_ys];
          Algebraic respective_x = alg_line.x_at_y(current_y);
          Alg_point_2 intersection(respective_x, current_y);
          CGAL_assertion(this->has_on(intersection));
          *o++ = intersection;
        }

        return o; // already one past the end, post-incremented when adding
      }
      /* in the general case we substitute the y in the parabola equation with
       * the value:                             y = -c/b + -ax/b
       * We get a quadratic equation in x, which we can solve using CGAL */
      else {
        /* the quadratic equation in x is:
         *                (r + (sa^2 / b^2) - (at / b)) x^2
         *   + ((2acs / b^2) - (ct / b) - (av / b) + u) x
         *              + ((sc^2 / b^2) - (cv / b) + w)
         *                                                  = 0
         */
        RT EQ_A = this->r() +
          (this->s() * CGAL::square(a) / CGAL::square(b)) -
          (a * this->t() / b)
        ;
        RT EQ_B = (2 * a * c * this->s() / CGAL::square(b)) -
          (c * this->t() / b) -
          (a * this->v() / b) +
          this->u()
        ;
        RT EQ_C = (this->s() * CGAL::square(c) / CGAL::square(b)) -
          (c * this->v() / b) +
          this->w()
        ;

        /* to store the 0, 1, or 2 results, indicating the intersections.
         * For all resulting x, find the corresponding y, and add a point to
         * the OutputIterator o */
        Algebraic  xs[2];
        Algebraic * xs_end;
        int n_xs;
        xs_end = nt_traits.solve_quadratic_equation(EQ_A, EQ_B, EQ_C, xs);
        n_xs = xs_end - xs;

        /* if no intersections return */
        if (n_xs == 0) {
          return o;
        }
        /* else find xs for all xs, add points to iterator */
        else while (--n_xs >= 0) {
          Algebraic current_x = xs[n_xs];
          Algebraic respective_y = alg_line.y_at_x(current_x);
          Alg_point_2 intersection(current_x, respective_y);
          CGAL_assertion(this->has_on(intersection));
          *o++ = intersection;
        }

        return o; // already one past the end, post-incremented when adding
      }
    }

    /* Given a point on the parabola and a vector of lines, finds all
     * intersections of the parabola with those lines, then looks for the next
     * point on the parabola (in the same direction as the directrix) after the
     * given start point and returns that point. */
    Alg_point_2 next_intersection(
      Alg_point_2 start,
      std::vector<Rat_line_2> delimiters
    ) {
      /* get intersections */
      std::list<Alg_point_2> intersections;
      for (auto& delimiter : delimiters) {
        this->get_intersections(delimiter, std::back_inserter(intersections));
      }
      CGAL_assertion(intersections.size() > 0);

      /* convert directrix and project start point on it */
      RK_to_AK to_alg;
      Alg_line_2 alg_directrix = to_alg(this->directrix());
      Alg_point_2 alg_start_pt = alg_directrix.projection(start);

      /* project intersections on alg_directrix, filter points that are "before"
       * alg_start_pt on the directrix */
      std::list<Alg_point_2> relevant_proj_intersections;
      for (auto& p : intersections) {
        Alg_point_2 proj_p = alg_directrix.projection(p);
        if (
          alg_directrix.direction()
          ==
          Alg_segment_2(alg_start_pt, proj_p).direction()
        ) {
          relevant_proj_intersections.push_back(proj_p);
        }
      };
      CGAL_assertion(relevant_proj_intersections.size() > 0);

      /* find next intersection on directrix, return corresponding point on the
       * parabola. Just look through the intersections and find it. */
      Alg_point_2 next_on_directrix = closest_point<Alg_kernel>(
        alg_start_pt,
        relevant_proj_intersections
      );
      Alg_point_2 result; bool assigned;
      for (auto& intersection : intersections) {
        if (alg_directrix.projection(intersection) == next_on_directrix) {
          result = intersection;
          assigned = true;
        }
      }

      CGAL_assertion_msg(assigned, "Could not find closest point");
      return result;
    }

    /* Construct a parabolic arc on the parabola from point p1 to point p2.
     * Precondition (checked): p1 and p2 are on the parabola */
    Curve_2 construct_parabolic_arc(Point_2 p1, Point_2 p2) {
      /* check precondition: both points lie on the parabola */
      CGAL_assertion(this->has_on(p1));
      CGAL_assertion(this->has_on(p2));

      /* construct the curve using the parameters and the endpoints */
      Curve_2 arc(_r,_s,_t,_u,_v,_w, this->orientation(), p1, p2);
      CGAL_assertion(arc.is_valid()); // valid arc
      return arc;
    }
  };

/* ########################################################################## */
/* ###                        END OF CLASS PARABOLA                       ### */
/* ########################################################################## */



  /* Returns the squared distance between two points in L2 metric. */
  static Algebraic sqdistance(const Point_2& p1, const Point_2& p2) {
    return CGAL::squared_distance(p1, p2);
  }
  static Rational sqdistance(const Rat_point_2& p1, const Rat_point_2& p2) {
    return CGAL::squared_distance(p1, p2);
  }

  /* Returns the squared distance between a point and a segment in L2 metric. */
  static Algebraic sqdistance(const Point_2& p, const Rat_segment_2& s) {
    /* if segment is degenerate (a point), call other function */
    if (s.is_degenerate()) {
      RK_to_AK to_alg;
      return sqdistance(p, to_alg(s.source()));
    }

    /* find projection of p on supporting line of s */
    Alg_segment_2 alg_s(
      Alg_point_2(s.source().x(), s.source().y()),
      Alg_point_2(s.target().x(), s.target().y())
    );
    Alg_line_2 l = alg_s.supporting_line();
    Alg_point_2 proj = l.projection(p);

    /* if the projection is on s, the distance is d(p,proj) */
    if (alg_s.has_on(proj)) {
      return sqdistance(p, proj);
    }
    /* otherwise, the distance is min(d(p,s1),d(p,s2)), where s1 and s2 are the
    endpoints of the segment s */
    else {
      return CGAL::min(
        sqdistance(p, Alg_point_2(s.source().x(), s.source().y())),
        sqdistance(p, Alg_point_2(s.target().x(), s.target().y()))
      );
    }
  }

  /* Given a point p and a list of points, return the closest point.
   * Precondition (checked): the list of points is not empty */
  template <class K> // Kernel
  static typename K::Point_2 closest_point(
    typename K::Point_2 p,
    std::list<typename K::Point_2> points
  ) {
    CGAL_precondition_msg(points.size() > 0, "List of points cannot be empty.");
    typename K::Point_2 result;
    typename K::FT smaller_sqdistance = -1;
    for (auto& q : points) {
      typename K::FT sqdist_pq = sqdistance(p, q);
      if (smaller_sqdistance < 0 || smaller_sqdistance > sqdist_pq) {
        result = q;
        smaller_sqdistance = sqdist_pq;
      }
    }

    CGAL_assertion(smaller_sqdistance >= 0);
    return result;
  }

  /* Given an Algebraic point and a Quadrant, move the point in the in the
   * general direction of that quadrant (e.g. 2nd quadrant == decrease x,
   * increase y) and return it as algebraic. */
  static Rat_point_2 move_point_alg_to_rat(Alg_point_2 p, Quadrant q) {
    Rat_point_2 rp;
    switch (q) {
      case FIRST_Q: {
        rp = Rat_point_2(
          approximate_algebraic(p.x(), true),   // increase x
          approximate_algebraic(p.y(), true)   // increase y
        );
        break;
      }
      case SECOND_Q: {
        rp = Rat_point_2(
          approximate_algebraic(p.x(), false),  // decrease x
          approximate_algebraic(p.y(), true)   // increase y
        );
        break;
      }
      case THIRD_Q: {
        rp = Rat_point_2(
          approximate_algebraic(p.x(), false),  // decrease x
          approximate_algebraic(p.y(), false)  // decrease y
        );
        break;
      }
      case FOURTH_Q: {
        rp = Rat_point_2(
          approximate_algebraic(p.x(), true),   // increase x
          approximate_algebraic(p.y(), false)  // decrease y
        );
        break;
      }
      default: {
        CGAL_error_msg(
          "This function is not working properly: move_point_alg_to_rat."
        );
        break;
      }
    }
    return rp;
  }

  /* Given a segment, determine the Quadrant where it is oriented "into".
   * Precondition (checked): the segment is not degenerate.
   * Return the Quadrant. */
  template <class K> // Kernel
  static Quadrant general_direction(typename K::Segment_2 s) {
    CGAL_precondition_msg(!s.is_degenerate(), "Segment cannot be degenerate.");

    /* get direction */
    typename K::Direction_2 s_direction = s.direction();
    typename K::RT dx = s_direction.dx(), dy = s_direction.dy();

    /* for each quadrant, always include the "starting" vertical or horizontal
     * direction of the quadrant, but exclude the last direction */
    if (dx > 0 && dy >= 0) {
      return FIRST_Q;
    }
    else if (dx <= 0 && dy > 0) {
      return SECOND_Q;
    }
    else if (dx < 0 && dy <= 0) {
      return THIRD_Q;
    }
    else if (dx >= 0 && dy < 0) {
      return FOURTH_Q;
    }
    else {
      CGAL_error_msg(
        "This function is not working properly: general_direction."
      );
    }
  }

  /* Given an Algebraic number, approximate it to a Rational according to the
   * given flag (ceiling or floor) */
  static RT approximate_algebraic(AT n, bool up) {
    /* precision, better if power of 2 */
    BigInt precision(65536); // use string for more precision

    /* multiply number by precision, take ceil or floor */
    AT big_n = n * precision;
    BigInt approximate_big_n;
    if (up) approximate_big_n = big_n.BigIntValue() + 1;  // ceil
    else approximate_big_n = big_n.BigIntValue() - 1;     // floor

    /* create rational with approximate_big_n and precision */
    RT result(approximate_big_n, precision);
    return result;
  }

  /* Given a segment, rotate it very slightly clockwise or counterclockwise
   * according to the flag. It doesn't matter if the length changes.
   * The supporting line must keep approximately the same orientation.
   * Return the rotated segment, it has rational coordinates. */
  static Rat_segment_2 slightly_rotate_segment(Alg_segment_2 s, bool ccw) {
    /* depending on towards which Quadrant the segment is oriented, rotate it
     * by changing individually each one of its coordinates according to the
     * given flag ccw (if true, rotate counterclockwise) */
    switch (general_direction<Alg_kernel>(s)) {
      case FIRST_Q: {
        return Rat_segment_2(
          Rat_point_2(
            approximate_algebraic(s.source().x(), /* ceil or floor */ ccw),
            approximate_algebraic(s.source().y(), /* ceil or floor */ !ccw)
          ),
          Rat_point_2(
            approximate_algebraic(s.target().x(), /* ceil or floor */ !ccw),
            approximate_algebraic(s.target().y(), /* ceil or floor */ ccw)
          )
        );
        break;
      }
      case SECOND_Q: {
        return Rat_segment_2(
          Rat_point_2(
            approximate_algebraic(s.source().x(), /* ceil or floor */ ccw),
            approximate_algebraic(s.source().y(), /* ceil or floor */ ccw)
          ),
          Rat_point_2(
            approximate_algebraic(s.target().x(), /* ceil or floor */ !ccw),
            approximate_algebraic(s.target().y(), /* ceil or floor */ !ccw)
          )
        );
        break;
      }
      case THIRD_Q: {
        return Rat_segment_2(
          Rat_point_2(
            approximate_algebraic(s.source().x(), /* ceil or floor */ !ccw),
            approximate_algebraic(s.source().y(), /* ceil or floor */ ccw)
          ),
          Rat_point_2(
            approximate_algebraic(s.target().x(), /* ceil or floor */ ccw),
            approximate_algebraic(s.target().y(), /* ceil or floor */ !ccw)
          )
        );
        break;
      }
      case FOURTH_Q: {
        return Rat_segment_2(
          Rat_point_2(
            approximate_algebraic(s.source().x(), /* ceil or floor */ !ccw),
            approximate_algebraic(s.source().y(), /* ceil or floor */ !ccw)
          ),
          Rat_point_2(
            approximate_algebraic(s.target().x(), /* ceil or floor */ ccw),
            approximate_algebraic(s.target().y(), /* ceil or floor */ ccw)
          )
        );
        break;
      }
      default: break; // should never happen
    }
  }

  /* Given a segment, translate it very slightly up or down according to the
   * given flag. It doesn't matter if the length changes.
   * The supporting line must keep approximately the same orientation.
   * Return the translated segment, it has rational coordinates. */
  static Rat_segment_2 slightly_translate_segment(Alg_segment_2 s, bool up) {
    /* depending on towards which Quadrant the segment is oriented, translate it
     * by changing individually each one of its coordinates according to the
     * given flag up (if true, translate towards positive side) */
    switch (general_direction<Alg_kernel>(s)) {
      case FIRST_Q: {
        return Rat_segment_2(
          Rat_point_2(
            approximate_algebraic(s.source().x(), /* ceil or floor */ !up),
            approximate_algebraic(s.source().y(), /* ceil or floor */ up)
          ),
          Rat_point_2(
            approximate_algebraic(s.target().x(), /* ceil or floor */ !up),
            approximate_algebraic(s.target().y(), /* ceil or floor */ up)
          )
        );
        break;
      }
      case SECOND_Q: {
        return Rat_segment_2(
          Rat_point_2(
            approximate_algebraic(s.source().x(), /* ceil or floor */ !up),
            approximate_algebraic(s.source().y(), /* ceil or floor */ !up)
          ),
          Rat_point_2(
            approximate_algebraic(s.target().x(), /* ceil or floor */ !up),
            approximate_algebraic(s.target().y(), /* ceil or floor */ !up)
          )
        );
        break;
      }
      case THIRD_Q: {
        return Rat_segment_2(
          Rat_point_2(
            approximate_algebraic(s.source().x(), /* ceil or floor */ up),
            approximate_algebraic(s.source().y(), /* ceil or floor */ !up)
          ),
          Rat_point_2(
            approximate_algebraic(s.target().x(), /* ceil or floor */ up),
            approximate_algebraic(s.target().y(), /* ceil or floor */ !up)
          )
        );
        break;
      }
      case FOURTH_Q: {
        return Rat_segment_2(
          Rat_point_2(
            approximate_algebraic(s.source().x(), /* ceil or floor */ up),
            approximate_algebraic(s.source().y(), /* ceil or floor */ up)
          ),
          Rat_point_2(
            approximate_algebraic(s.target().x(), /* ceil or floor */ up),
            approximate_algebraic(s.target().y(), /* ceil or floor */ up)
          )
        );
        break;
      }
      default: break; // should never happen
    }
  }

  /* Given a segment with one algebraic endpoint and potentially an underlying
   * supporting line with algebraic coefficients, adjust the algebraic endpoint
   * of the segment up or down to a rational point.
   * `source` indicates if the source of the segment that is algebraic or not
   * (if it is false, it meanst that the target is algebraic)
   * `up` indicates whether the endpoint needs to be moved to the positive (up)
   * or negative (!up) part of the segment.
   * As last argument, pass the exact point that is the endpoint that doesn't
   * have to be approximated.
   * Return a rational segment.
   */
  static Rat_segment_2 adjust_endpoint(
    Alg_segment_2 s,
    bool up,
    bool source,
    Rat_point_2 exact_endpoint
  ) {
    /* determine direction to approximate (in the form of a Quadrant) */
    Quadrant quadrant_direction_to_move_point_towards;
    switch (general_direction<Alg_kernel>(s)) {
      case FIRST_Q: {
        if (up) quadrant_direction_to_move_point_towards = SECOND_Q;
        else quadrant_direction_to_move_point_towards = FOURTH_Q;
        break;
      }
      case SECOND_Q: {
        if (up) quadrant_direction_to_move_point_towards = THIRD_Q;
        else quadrant_direction_to_move_point_towards = FIRST_Q;
        break;
      }
      case THIRD_Q: {
        if (up) quadrant_direction_to_move_point_towards = FOURTH_Q;
        else quadrant_direction_to_move_point_towards = SECOND_Q;
        break;
      }
      case FOURTH_Q: {
        if (up) quadrant_direction_to_move_point_towards = FIRST_Q;
        else quadrant_direction_to_move_point_towards = THIRD_Q;
        break;
      }
    }

    /* approximate the correct endpoint */
    Alg_point_2 endpoint_to_approximate = source ? s.source() : s.target();
    Rat_point_2 approximated_endpoint = move_point_alg_to_rat(
      endpoint_to_approximate, quadrant_direction_to_move_point_towards
    );

    Rat_segment_2 approximated_s;
    if (source) {
      approximated_s = Rat_segment_2(approximated_endpoint, exact_endpoint);
    }
    else {
      approximated_s = Rat_segment_2(exact_endpoint, approximated_endpoint);
    }

    return approximated_s;
  }

  /* Given two segments return true if they are collinear */
  static bool segments_collinear(Rat_segment_2 s1, Rat_segment_2 s2) {
    return
      CGAL::collinear(s1.source(), s1.target(), s2.source())
      &&
      CGAL::collinear(s2.source(), s2.target(), s1.source())
    ;
  }

  /* Given two Curve_2 objects, return true if they are the same */
  static bool same_curves(Curve_2 cv1, Curve_2 cv2) {
    return
      (cv1.r() == cv2.r()) &&
      (cv1.s() == cv2.s()) &&
      (cv1.t() == cv2.t()) &&
      (cv1.u() == cv2.u()) &&
      (cv1.v() == cv2.v()) &&
      (cv1.w() == cv2.w()) &&
      (cv1.orientation() == cv2.orientation()) &&
      (cv1.source() == cv2.source()) &&
      (cv1.target() == cv2.target())
    ;
  }

  /* Construct a point in the middle of the curve cv. This function is copied
   * from Env_sphere_traits_3.h */
  static Point_2 construct_middle_point(const X_monotone_curve_2& cv) {
    Point_2 result;

    /* get the x-value of the middle point */
    Alg_point_2 mid_x = CGAL::midpoint(cv.source(),cv.target());

    /* if cv is vertical, it is just a segment */
    if (cv.is_vertical()) result = Point_2(mid_x);
    /* otherwise take the point with the same x coordinate but on cv */
    else result = Point_2(cv.point_at_x(mid_x));

    return result;
  }

  /* Given a ray, find an endpoint as the intersection with an imaginary
   * boundary (for example a box at 10000 around the origin).
   * Create a segment from the source of the ray and invert it if the flag is
   * set to true.
   * This is used so to be able to store the ray as a Curve_2, that in the
   * traits class we are using, Arr_conic_traits_2, can only be a bounded curve.
   * Return the segment. */
  static Rat_segment_2 make_segment_from_ray(Rat_ray_2 ray, bool invert) {

    //TODO change, can be smarter than brute force

    RT far_l = 10000;
    std::vector<Rat_segment_2> border = {
      Rat_segment_2(Rat_point_2(-far_l, -far_l), Rat_point_2(far_l, -far_l)),
      Rat_segment_2(Rat_point_2(far_l, -far_l), Rat_point_2(far_l, far_l)),
      Rat_segment_2(Rat_point_2(far_l, far_l), Rat_point_2(-far_l, far_l)),
      Rat_segment_2(Rat_point_2(-far_l, far_l), Rat_point_2(-far_l, -far_l))
    };

    Rat_point_2 start_point = ray.source();
    Rat_point_2 end_point;
    bool assigned = false;
    for (auto& seg : border) {
      if (CGAL::do_intersect(ray, seg) && !assigned) {
        assigned = true;
        CGAL_assertion_msg(
          CGAL::assign(end_point, CGAL::intersection(ray, seg)),
          "Could not assign end."
        );
        break;
      }
    }
    CGAL_assertion_msg(assigned, "Could not find ray end_point.");

    /* make segment, invert if needed and return */
    Rat_segment_2 segment(start_point, end_point);
    if (invert) segment = segment.opposite();
    return segment;
  }

  /* Convert the Curve_2 cv into multiple X_monotone_curve_2 using the provided
   * make_x_monotone function. Store the results into the list x_mono_curves.
   * Precondition (checked): cv is a valid curve. */
  static void make_curve_2_into_many_x_monotone_curve_2(Curve_2& cv,
    std::vector<X_monotone_curve_2>& x_mono_curves) {
    /* check precondition */
    CGAL_precondition_msg(
      cv.is_valid(),
      "The given curve cannot be converted to X_monotone parts because it is "
      "not valid"
    );

    /* instantiate traits, we need the provided function */
    C_traits_2 c_traits;
    typename C_traits_2::Make_x_monotone_2 make_x_monotone =
      c_traits.make_x_monotone_2_object();

    /* call the provided function */
    std::vector<CGAL::Object> pre_x_mono_curves;
    make_x_monotone(cv, std::back_inserter(pre_x_mono_curves));

    /* cast all CGAL::Objects into X_monotone_curve_2 and add to list */
    for(size_t i = 0; i < pre_x_mono_curves.size(); i++ ) {
      X_monotone_curve_2 curr;
      bool check = CGAL::assign(curr, pre_x_mono_curves[i]);
      assert(check); CGAL_USE(check);
      x_mono_curves.push_back(curr);
    }
  }

  /* Check if a given segment called edge connects two the other segments s1 and
   * s2 by any of their endpoints.
   * Return true if this is the case, return false if edge is actyally just
   * connecting s1's endpoints (or s2's). We cannot just check for equality
   * because edge could be just the same as one segment but in the other
   * direction. */
  static bool edge_connects_segments(Rat_segment_2 edge, Rat_segment_2 s1,
    Rat_segment_2 s2) {
    /* create a copy of edge but in the other direction, then check equality for
     * both versions of edge */
    Rat_segment_2 rev_edge(edge.target(), edge.source());
    return !(edge == s1 || edge == s2 || rev_edge == s1 || rev_edge == s2);
  }

  /* Given a bisector finds the point that is the furthest intersection
   * (following the direction of the bisector) of the bisector with the four
   * lines saved in delimiters.
   * Return the unbounded ray. */
  static Rat_ray_2 find_unbounded_ray(
    Rat_line_2 bisector,
    Rat_delimiter_lines delimiters
  ) {
    /* get the four intersection points, add them to two lists, sort one by x
     * and the other one by y */
    Rat_point_2 p1; // intersection between bisector and delimiters[0][0]
    Rat_point_2 p2; // intersection between bisector and delimiters[0][1]
    Rat_point_2 p3; // intersection between bisector and delimiters[1][0]
    Rat_point_2 p4; // intersection between bisector and delimiters[1][1]
    CGAL::assign(p1, CGAL::intersection(bisector, delimiters.first.first));
    CGAL::assign(p2, CGAL::intersection(bisector, delimiters.first.second));
    CGAL::assign(p3, CGAL::intersection(bisector, delimiters.second.first));
    CGAL::assign(p4, CGAL::intersection(bisector, delimiters.second.second));
    std::vector<Rat_ray_2> intersections_x = {
      Rat_ray_2(p1, bisector.direction()),
      Rat_ray_2(p2, bisector.direction()),
      Rat_ray_2(p3, bisector.direction()),
      Rat_ray_2(p4, bisector.direction())
    };
    std::vector<Rat_ray_2> intersections_y;
    std::copy(intersections_x.begin(), intersections_x.end(),
              std::back_inserter(intersections_y));
    std::sort(intersections_x.begin(), intersections_x.end(),
      [](Rat_ray_2 a, Rat_ray_2 b) {
      return a.source().x() < b.source().x();
    });
    std::sort(intersections_y.begin(), intersections_y.end(),
      [](Rat_ray_2 a, Rat_ray_2 b) {
      return a.source().y() < b.source().y();
    });

    /* find the farthest point according to the direction of bisector */
    Rat_direction_2 dir = bisector.direction();
    if (dir.dx() == 0) {
      return (dir.dy() > 0) ? intersections_y.back() : intersections_y.front();
    }
    else {
      return (dir.dx() > 0) ? intersections_x.back() : intersections_x.front();
    }
  }


  /* Given a direction and a start point finds the point that is the first
   * intersection (after this start point with the four lines saved in the
   * vector of delimiter lines.
   * Return the found intersection point. */
  static Alg_point_2 find_next_intersection(
    Alg_direction_2 direction,
    Alg_point_2 start_pt,
    std::vector<Rat_line_2> delimiters
  ) {
    /* converter */
    RK_to_AK to_alg;
    /* list to store intersections */
    std::list<Alg_point_2> intersections;

    /* for each delimiter add the intersection with ray, if it's not start_pt */
    Alg_ray_2 ray(start_pt, direction);
    for (auto& delimiter : delimiters) {
      if (CGAL::do_intersect(to_alg(delimiter), ray)) {
        Alg_point_2 intersection;
        CGAL::assign(intersection, CGAL::intersection(to_alg(delimiter), ray));
        if (intersection != start_pt) {
          intersections.push_back(intersection);
        }
      }
    }

    CGAL_assertion_msg(intersections.size() > 0, "There must be intersections");

    /* all intersections are in the correct direction because we used a ray
     * starting from start_pt, so return the closest one */
    return closest_point<Alg_kernel>(start_pt, intersections);
  }

  /* Determine the position of the point p relative to the segments s1 and s2.
   * We do not have the segments though: we have four lines, each orthogonal to
   * one endpoint of one of the two segments. Every line is oriented so to have
   * the inner part of the segment on their negative side (right side).
   * There are three main cases (described by enum Bisector_type):
   * - PARABOLIC_ARC: when p is closer to one segment's inner part and to one of
   *   the other segment's endpoints. In this case, save the supporting_line of
   *   the first segment and the endpoint of the second segment.
   * - SUPP_LINE_BISECTOR: when p is closer to both inner parts of both the two
   *   segments. In this case save the two supporting_lines of the two segments.
   * - ENDPOINT_BISECTOR: when p is closer to two endpoints of the two segments.
   *   In this case save those two endpoints.
   * In all three cases we save in o1 the correct endpoint or supporting_line of
   * s1, and in o2 the same for s2.
   * Precondition (checked): the point p is not on any of the four delimiters
   */
  static Bisector_type find_position(
    Alg_point_2 p,
    Alg_delimiter_lines delimiter_lines,
    Rat_segment_2 s1,
    Rat_segment_2 s2,
    Object& o1,  // to store [directrix1 or focus1]/line1/point1
    Object& o2   // to store [directrix2 or focus2]/line2/point2
  ) {
    /* check precondition */
    CGAL_precondition(!delimiter_lines.first.first.has_on(p));
    CGAL_precondition(!delimiter_lines.first.second.has_on(p));
    CGAL_precondition(!delimiter_lines.second.first.has_on(p));
    CGAL_precondition(!delimiter_lines.second.second.has_on(p));

    /* assume point is not on any delimiter, consider all other cases. To do so,
     * first determine what must be stored in o1, then in o2. Save in two flags
     * information about the case. In the end, determine the case. */
    bool o1_is_line = false;
    bool o2_is_line = false;

    /* determine o1 */
    if (delimiter_lines.first.first.has_on_positive_side(p)) {
      o1 = CGAL::make_object(s1.source());
    }
    else if (delimiter_lines.first.second.has_on_positive_side(p)) {
      o1 = CGAL::make_object(s1.target());
    }
    else {
      CGAL_assertion(
        delimiter_lines.first.first.has_on_negative_side(p)
        &&
        delimiter_lines.first.second.has_on_negative_side(p)
      );
      o1 = CGAL::make_object(s1.supporting_line());
      o1_is_line = true;
    }

    /* determine o2 */
    if (delimiter_lines.second.first.has_on_positive_side(p)) {
      o2 = CGAL::make_object(s2.source());
    }
    else if (delimiter_lines.second.second.has_on_positive_side(p)) {
      o2 = CGAL::make_object(s2.target());
    }
    else {
      CGAL_assertion(
        delimiter_lines.second.first.has_on_negative_side(p)
        &&
        delimiter_lines.second.second.has_on_negative_side(p)
      );
      o2 = CGAL::make_object(s2.supporting_line());
      o2_is_line = true;
    }

    /* determine case using flags */
    if (o1_is_line) {
      return (o2_is_line) ? SUPP_LINE_BISECTOR : PARABOLIC_ARC;
    }
    else {
      return (o2_is_line) ? PARABOLIC_ARC : ENDPOINT_BISECTOR;
    }
  }

  /* Given a line and a direction determine whether the line is oriented in that
   * genreal direction, that is in a range of [-90˚, 90˚] around the Given
   * direction. */
  static bool generally_same_direction(Alg_line_2 line, Alg_direction_2 dir) {
    return line.direction().counterclockwise_in_between(
      dir.vector().perpendicular(CGAL::CLOCKWISE).direction(),
      dir.vector().perpendicular(CGAL::COUNTERCLOCKWISE).direction()
    ); //TODO what if they are perpendicular? which case is it?
  }



public:

  class Make_xy_monotone_3 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Surface_3& s,
                                bool /* is_lower */,
                                OutputIterator o) const {
      /* the surfaces we are considering are distance functions from line
       * segments and are already xy_monotone because there is only one possible
       * distance value for any point on the plane */
      *o++ = s; // just insert the surface in o, return o one past the end
      return o;
    }
  };

  Make_xy_monotone_3 make_xy_monotone_3_object() const {
    return Make_xy_monotone_3();
  }



  class Construct_projected_boundary_2 {
  public:
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const {
      // /* the surfaces we are considering are distance functions of line
      //  * segments and are infinite, so they have no projected boundary */
      // return o; // the iterator remains empty

      /* save boundary for intersection with it */
      RT far_l = 10000;
      std::vector<Rat_segment_2> border = {
        Rat_segment_2(Rat_point_2(-far_l, -far_l), Rat_point_2(far_l, -far_l)),
        Rat_segment_2(Rat_point_2(far_l, -far_l), Rat_point_2(far_l, far_l)),
        Rat_segment_2(Rat_point_2(far_l, far_l), Rat_point_2(-far_l, far_l)),
        Rat_segment_2(Rat_point_2(-far_l, far_l), Rat_point_2(-far_l, -far_l))
      };

      /* The surfaces representing distance functions are infinite, but to make
       * this work with bounded bisectors we need a boundary. We use the same
       * for all segments */
      for (auto& seg : border) {
        X_monotone_curve_2 x_seg = X_monotone_curve_2(seg);
        *o++ = CGAL::make_object(std::make_pair(x_seg, CGAL::ON_NEGATIVE_SIDE));
      }
      return o;
    }
  };

  Construct_projected_boundary_2
  construct_projected_boundary_2_object() const
  {
    return Construct_projected_boundary_2();
  }


/* ########################################################################## */
/* ###            MAIN PART: COMPUTE BISECTOR OF TWO SEGMENTS             ### */
/* ########################################################################## */

  class Construct_projected_intersections_2 {
  public:
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const {

      std::cout << "\n#######################################################"
                << "#######################################################\n"
                << "#\tFinding bisector of s1 = {"
                << s1 << "} and s2 = {" << s2 << "}:\n"
      ;

      /* create converter functors to convert from:
       * - Rational to Algebraic
       * - Algebraic to Cartesian<double>
       * - Cartesian<double> to Rational
       * The last two are used together to convert by approximation from
       * Algebraic to Rational */
      RK_to_AK to_alg;
      AK_to_DK to_dbl;
      DK_to_RK to_rat;

      /* this is a very bad solution, and will need to be changed, for example
       * by using Arr_algebraic_segment_traits_2 that supports both unbounded
       * and bounded curves (also left-/right-unbounded) */
      /* save boundary for intersection with it */
      RT far_l = 10000;
      std::vector<Rat_segment_2> border = {
        Rat_segment_2(Rat_point_2(-far_l, -far_l), Rat_point_2(far_l, -far_l)),
        Rat_segment_2(Rat_point_2(far_l, -far_l), Rat_point_2(far_l, far_l)),
        Rat_segment_2(Rat_point_2(far_l, far_l), Rat_point_2(-far_l, far_l)),
        Rat_segment_2(Rat_point_2(-far_l, far_l), Rat_point_2(-far_l, -far_l))
      };

      /* if the two segments are the same (also if one is just the other but
       * reversed), their distance function is the same, so there is no
       * intersection */
      if (s1 == s2 || s1 == Rat_segment_2(s2.target(), s2.source())) {
        return o;
      }
      /* if one of the segments is degenerate, the bisector is a parabola and
       * two rays, if instead they are both degenerate (that is, they are two
       * two points) the bisector is a line */
      else if (s1.is_degenerate() || s2.is_degenerate()) {
        /* line */
        if (s1.is_degenerate() && s2.is_degenerate()) {
          Rat_line_2 bisector = CGAL::bisector(s1.source(), s2.source());

          /* find endpoints */ //TODO abstract this
          Rat_point_2 start, end;
          bool assigned_start = false, assigned_end = false;
          for (auto& segment : border) {
            if (CGAL::do_intersect(bisector, segment)) {
              if (assigned_start && !assigned_end) {
                assigned_end = true;
                CGAL_assertion_msg(
                  CGAL::assign(end, CGAL::intersection(bisector, segment)),
                  "Could not assing end."
                );
              }
              else if (!assigned_start && !assigned_end) {
                assigned_start = true;
                CGAL_assertion_msg(
                  CGAL::assign(start, CGAL::intersection(bisector, segment)),
                  "Could not assing start."
                );
              }
            }
          }
          Rat_segment_2 seg_bisector(start, end);
          X_monotone_curve_2 x_seg_bisector = X_monotone_curve_2(seg_bisector);
          *o++ = CGAL::make_object(
            Intersection_curve(x_seg_bisector, 0)
          );
        }
        /* parabolic arc and two rays (or in other special cases just lines) */
        else {
          /* determine which one is the non-degenerate segment */
          bool s1_is_degenerate = s1.is_degenerate();
          CGAL_assertion_msg(
            (s1_is_degenerate != s2.is_degenerate()),
            "One segment should be degenerate and the other one not."
          );

          /* get directrix and focus */
          Rat_segment_2 directrix_generator = s1_is_degenerate ? s2 : s1;
          Rat_line_2 directrix = directrix_generator.supporting_line();
          Rat_point_2 focus = s1_is_degenerate ? s1.source() : s2.source();

          /* if segment and point collinear, the bisector is just the bisector
           * of the closest endpoint and the point */
          if (CGAL::collinear(
            focus,
            directrix_generator.source(),
            directrix_generator.target())
          ) {
            Rat_line_2 bisector;

            /* if the point is not on the segment just get the bisector */
            if (!directrix_generator.has_on(focus)) {

              /* get bisector */
              if (
                sqdistance(focus, directrix_generator.source())
                <
                sqdistance(focus, directrix_generator.target())
              ) {
                bisector = CGAL::bisector(focus, directrix_generator.source());
              }
              else {
                bisector = CGAL::bisector(focus, directrix_generator.target());
              }
            }
            /* otherwise the bisector is the line on orthogonal to the segment
             * that goes through the point on it */
            else {
              bisector = directrix.perpendicular(focus);
            }

            /* find endpoints */ //TODO abstract this part
            Rat_point_2 start, end;
            bool assigned_start = false, assigned_end = false;
            for (auto& segment : border) {
              if (CGAL::do_intersect(bisector, segment)) {
                if (assigned_start && !assigned_end) {
                  assigned_end = true;
                  CGAL_assertion_msg(
                    CGAL::assign(end, CGAL::intersection(bisector, segment)),
                    "Could not assing end."
                  );
                }
                else if (!assigned_start && !assigned_end) {
                  assigned_start = true;
                  CGAL_assertion_msg(
                    CGAL::assign(
                      start,
                      CGAL::intersection(bisector, segment)
                    ),
                    "Could not assing start."
                  );
                }
              }
            }
            Rat_segment_2 sg_bisector(start, end);
            X_monotone_curve_2 x_sg_bisector = X_monotone_curve_2(sg_bisector);
            *o++ = CGAL::make_object(
              Intersection_curve(x_sg_bisector, 0)
            );
          }
          /* otherwise, we have the actual parabola and the two rays */
          else {
            /* ... */ //TODO bisector of point and segment
          }

          /* make parabola */
          // *o++ = CGAL::make_object(
          //   Intersection_curve(curve_seg, 0)
          // );
        }
      }
      /* if the two segments lie on a single line, we just have a line as a
       * bisector */
      else if (segments_collinear(s1, s2)) {
        Rat_point_2 p1, p2;
        /* find the closest pair of points */
        if (s1.direction() == s2.direction()) {
          if (
            sqdistance(s1.source(), s2.target())
            <
            sqdistance(s1.target(), s2.source())
          ) {
            p1 = s1.source();
            p2 = s2.target();
          }
          else {
            p1 = s1.target();
            p2 = s2.source();
          }
        }
        else {
          if (
            sqdistance(s1.source(), s2.source())
            <
            sqdistance(s1.target(), s2.target())
          ) {
            p1 = s1.source(); p2 = s2.source();
          }
          else {
            p1 = s1.target(); p2 = s2.target();
          }
        }

        /* bisector is bisector of these two points */
        Rat_line_2 bisector = CGAL::bisector(p1, p2);
        Rat_point_2 start, end;
        bool assigned_start = false, assigned_end = false;
        for (auto& segment : border) {
          if (CGAL::do_intersect(bisector, segment)) {
            if (assigned_start && !assigned_end) {
              assigned_end = true;
              CGAL_assertion_msg(
                CGAL::assign(end, CGAL::intersection(bisector, segment)),
                "Could not assing end."
              );
            }
            else if (!assigned_start && !assigned_end) {
              assigned_start = true;
              CGAL_assertion_msg(
                CGAL::assign(start, CGAL::intersection(bisector, segment)),
                "Could not assing start."
              );
            }
          }
        }
        Rat_segment_2 seg_bisector(start, end);
        X_monotone_curve_2 x_seg_bisector = X_monotone_curve_2(seg_bisector);
        *o++ = CGAL::make_object(
          Intersection_curve(x_seg_bisector, 0)
        );
      }
      /* if the two segments are not the same and not degenerate, compute all
       * parts of their plane bisector */
      else {
        /* first of all, for each segment create the two lines that divide the
         * plane in three areas: one of all points closest to the inner part of
         * the segment, the other two of all points closest to the two endpoints
         * of the segment.
         * The lines are saved with an orientation such that they both have the
         * inner part of the segment on their right side (negative side).
         * Note: a vector constructed using a segment is oriented from source to
         * target of that segment, so to build a line such that the segment lies
         * on the right side of it, we need to use:
         * - the source of the segment and as direction the vector oriented 90
         *   degrees counterclockwise from the segment vector
         * - the target of the segment but as direction the vector oriented 90
         *   degrees clockwise. */
        Rat_delimiter_lines delimiter_lines = {
          {
            Rat_line_2(
              s1.source(),
              Rat_vector_2(s1).perpendicular(CGAL::COUNTERCLOCKWISE)
            ),
            Rat_line_2(
              s1.target(),
              Rat_vector_2(s1).perpendicular(CGAL::CLOCKWISE)
            )
          },
          {
            Rat_line_2(
              s2.source(),
              Rat_vector_2(s2).perpendicular(CGAL::COUNTERCLOCKWISE)
            ),
            Rat_line_2(
              s2.target(),
              Rat_vector_2(s2).perpendicular(CGAL::CLOCKWISE)
            )
          }
        };
        /* also save the lines in alg form for convenience */
        Alg_delimiter_lines alg_delimiter_lines = {
          {
            to_alg(delimiter_lines.first.first),
            to_alg(delimiter_lines.first.second)
          },
          {
            to_alg(delimiter_lines.second.first),
            to_alg(delimiter_lines.second.second)
          }
        };
        /* also save the lines in a vector for convenience */
        std::vector<Rat_line_2> delimiter_lines_vector = {
          delimiter_lines.first.first, delimiter_lines.first.second,
          delimiter_lines.second.first, delimiter_lines.second.second
        };
        /* also save segment endpoints "generating" these lines */
        std::vector<Rat_point_2> segment_endpoints = {
          s1.source(), s1.target(), s2.source(), s2.target()
        };

        /* then compute the 2 or 4 unbounded edges of the bisector.
         * To do this, first compute the convex hull of the endpoints of the
         * segments. The pairs of vertices of the hull that are not of the same
         * segment are the pairs of which the bisector lines contain the
         * unbounded rays that are the unbounded rays of the plane bisector of
         * the two segments */

        /* compute hull of endpoints; the hull points will be stored inside
         * ch_points in counterclockwise order */
        std::list<Rat_point_2> ch_points;
        std::list<Rat_point_2> points = {
          s1.source(), s1.target(), s2.source(), s2.target()
        };
        CGAL::ch_akl_toussaint(
          points.begin(), points.end(),
          std::back_insert_iterator<std::list<Rat_point_2>>(ch_points)
        );

        /* make a polygon out of the hull points, iterate over vertices to find
         * pairs to make rays, directed towards outside of polygon */
        Rat_polygon_2 ch_polygon(ch_points.begin(), ch_points.end());
        CGAL_assertion(ch_polygon.is_convex()); // it is a hull
        CGAL_assertion(ch_polygon.area() >= 0); // it is counterclockwise

        /* list to save the unbounded rays of the bisector */
        std::list<Rat_ray_2> unbounded_ray_list;

        for ( // for all edges
          Edge_iterator eit = ch_polygon.edges_begin();
          eit != ch_polygon.edges_end();
          ++eit
        ) {
          if (edge_connects_segments(*eit, s1, s2)) { // if it's not s1 or s2
            /* create line that bisects the segment, orient it outside */
            Rat_line_2 bisector_line = CGAL::bisector(
              eit->target(), eit->source()
            );

            /* find "farthest" intersection with delimiters, ray starts there */
            Rat_ray_2 unbounded_ray = find_unbounded_ray(
              bisector_line, delimiter_lines
            );

            /* add ray to list of unbounded rays */
            unbounded_ray_list.push_back(unbounded_ray);
          }
        }

        /* if the two segments do NOT intersect, construct the bisector starting
         * from one unbounded edge, finding the correct intersection points
         * using the delimiter_lines.
         * In this case, the ray start points should be only two. */
        if (!CGAL::do_intersect(s1, s2)) { // segments do not intersect
          CGAL_assertion(unbounded_ray_list.size() == 2);

          /* starting from the source of one unbounded ray and finishing at the
           * source of the other, compute the rest of the bisector, consisting
           * of:
           * - parabolic arcs: when we are in the "area of influence" of the
           *   interior of a segment and of one endpoint of the other
           * - segments: when we are in the "area of influence" of the interiors
           *   of the two segments or of two endpoints of the two segments */
          Rat_ray_2 start_ray = unbounded_ray_list.front();
          Rat_ray_2 end_ray = unbounded_ray_list.back();

          Alg_direction_2 curr_direction = to_alg(
            - start_ray.direction()
          );

          /* list to store all parts of the bisector */
          std::list<Curve_2> parts_of_bisector;

          /* add start ray (as a long segment) */
          parts_of_bisector.push_back(make_segment_from_ray(start_ray, true));

          /* call big private function that iteratively constructs the parts of
           * the plane bisector of the segments s1 and s2 starting from a start
           * point going in a given direction and finishing at an end point */
          this->construct_bisector_from_point_to_point(
            s1, s2,                 // the two segments
            std::back_inserter(parts_of_bisector),  // OutputIterator
            start_ray, end_ray,     // construct bisector between two sources
            curr_direction,         // initial direction, updated
            alg_delimiter_lines,    // delimiter lines of s1 and s2
            delimiter_lines_vector  // same but as vector and in rational
          );

          /* add end ray (as a long segment) */
          parts_of_bisector.push_back(make_segment_from_ray(end_ray, false));

          /* convert and add curves to OutputIterator o */
          o = this->convert_and_add_curves(parts_of_bisector, o);

        } // end of segments do not intersect

        /* if instead they do intersect, assert it, then proceed to computing
         * the bisector in this case.
         * In this case, the ray start points should be four, but only if the
         * intersection is not by one or two endpoints (weak intersection) */
        else {
          CGAL_assertion(CGAL::do_intersect(s1, s2)); // they HAVE to intersect

          //TODO add check for touching segments, maybe even before (also add
          // check for three collinear points among the segments)

          /* starting from the source of one unbounded ray and finishing at the
           * source of another, compute the rest of the bisector, consisting of:
           * - parabolic arcs: when we are in the "area of influence" of the
           *   interior of a segment and of one endpoint of the other
           * - segments: when we are in the "area of influence" of the interiors
           *   of the two segments or of two endpoints of the two segments
           * Since in this case we have four rays, do this process twice, and we
           * create two main parts of the bisector that "touch" at the
           * intersection of the two segments. */
          CGAL_assertion(unbounded_ray_list.size() == 4);
          Rat_ray_2 start_ray_one = unbounded_ray_list.front();
          unbounded_ray_list.pop_front();
          Rat_ray_2 end_ray_one = unbounded_ray_list.front();
          unbounded_ray_list.pop_front();
          Rat_ray_2 start_ray_two = unbounded_ray_list.front();
          unbounded_ray_list.pop_front();
          Rat_ray_2 end_ray_two = unbounded_ray_list.front();
          unbounded_ray_list.pop_front();

          CGAL_assertion_msg(
            unbounded_ray_list.empty(),
            "There should only be 4 rays, and all were popped."
          );

          /* get direction for both cases */
          Alg_direction_2 curr_direction_one = to_alg(
            - start_ray_one.direction()
          );
          Alg_direction_2 curr_direction_two = to_alg(
            - start_ray_two.direction()
          );

          /* lists to store all parts of the bisector */
          std::list<Curve_2> parts_of_bisector_one;
          std::list<Curve_2> parts_of_bisector_two;

          /* add start rays (as long segments) */
          parts_of_bisector_one.push_back(make_segment_from_ray(
            start_ray_one, true
          ));
          parts_of_bisector_two.push_back(make_segment_from_ray(
            start_ray_two, true
          ));

          /* call big private function that iteratively constructs the parts of
           * the bisector of s1 and s2 starting from a start point going in a
           * given direction and finishing at an end point. Do this two times,
           * for both bisectors */
          this->construct_bisector_from_point_to_point(
            s1, s2,                     // the two segments
            std::back_inserter(parts_of_bisector_one),  // OutputIterator
            start_ray_one, end_ray_one, // construct bisector between sources
            curr_direction_one,         // initial direction, updated
            alg_delimiter_lines,        // delimiter lines of s1 and s2
            delimiter_lines_vector      // same but as vector and in rational
          );
          this->construct_bisector_from_point_to_point(
            s1, s2,                     // the two segments
            std::back_inserter(parts_of_bisector_two),  // OutputIterator
            start_ray_two, end_ray_two, // construct bisector between sources
            curr_direction_two,         // initial direction, updated
            alg_delimiter_lines,        // delimiter lines of s1 and s2
            delimiter_lines_vector      // same but as vector and in rational
          );

          /* add end rays (as long segments) */
          parts_of_bisector_one.push_back(make_segment_from_ray(
            end_ray_one, false
          ));
          parts_of_bisector_two.push_back(make_segment_from_ray(
            end_ray_two, false
          ));

          /* convert and add curves to OutputIterator o */
          // o = this->convert_and_add_curves(parts_of_bisector_one, o);
          // o = this->convert_and_add_curves(parts_of_bisector_two, o);
          /* do it this way cuz it's cooler */
          o = this->convert_and_add_curves(
            parts_of_bisector_two,
            this->convert_and_add_curves(parts_of_bisector_one, o)
          );

        } // end of segments intersect

      } // end of segments are not the same

      /* return one past the end iterator */
      return o;
    }

  private:

    /* Given a two Curve_2 parabolic arcs and an algebraic segment that
     * connects them, find an approximated segment supported by a line with
     * rational coefficients. The supporting MUST intersect both arcs.
     * The arcs are tangent to the segment to approximate at its endpoints.
     * Precondition (not checked): the arcs and the segment are oriented in the
     * correct way, so that they are directed from the source of prev_arc to the
     * target of next_arc.
     * Return this supporting line.
     */
    Rat_line_2 get_approximated_inner_segment_supporting_line(
      Alg_segment_2& segment,
      Curve_2& prev_arc,
      Curve_2& next_arc
    ) const {
      /* determine if we have to rotate or translate: if the orientation of the
       * previous and of the next arc is the same, we have to translate the
       * algebraic segment up or down, if the orientations are different instead
       * we have to slightly rotate the segment.
       * This is to ensure that the supporting_line line of the approximated
       * Rational segment still intersects both curves.
       *
       * There are two main cases:
       * - the curves have the same orientation. In this case, the approximated
       *   segment needs to be slightly moved towards the "interior" of both
       *   curves, so that it definitely intersects still both of them
       * - the curves have a different orietnat. In this case, the approximated
       *   segment needs to be slightly rotated so that the source is in the
       *   "interior" of the prev_arc, and the target of the next_arc
       *
       * To approximate, we work on the single coordinates of the points, so we
       * need to create four Rational numbers.
       * According to which direction each point needs to be moved, we multiply
       * the Algebraic coordinate by a large multiple of 2 (say, 2^16), then we
       * take the ceiling or floor of this number so to get an Integer, then
       * create a Rational that is this integer divided by the large multiple of
       * 2 (again, for example 2^16).
       */
      Rat_segment_2 approximated_segment;

      /* for checks */
      RK_to_AK to_alg;

      /* rotate */
      if (prev_arc.orientation() != next_arc.orientation()) {
        /* rotate segment counterclockwise */
        if (prev_arc.orientation() == CGAL::CLOCKWISE) {
          approximated_segment = slightly_rotate_segment(segment, true);
        }
        /* rotate segment clockwise */
        else {
          approximated_segment = slightly_rotate_segment(segment, false);
        }

      }
      /* translate "up or down" */
      else {
        /* move "up" (positive side) */
        if (prev_arc.orientation() == CGAL::COUNTERCLOCKWISE) {
          approximated_segment = slightly_translate_segment(segment, true);
          /* check */
          Alg_segment_2 alg_approx_seg = to_alg(approximated_segment);
          CGAL_assertion_msg(
            (segment.supporting_line().has_on_positive_side(
              alg_approx_seg.source()
            )
            &&
            segment.supporting_line().has_on_positive_side(
              alg_approx_seg.target()
            )),
            "The approximated segment should be above initial segment."
          );
        }
        /* move "down" (negative side) */
        else {
          approximated_segment = slightly_translate_segment(segment, false);
          /* check */
          Alg_segment_2 alg_approx_seg = to_alg(approximated_segment);
          CGAL_assertion_msg(
            (segment.supporting_line().has_on_negative_side(
              alg_approx_seg.source()
            )
            &&
            segment.supporting_line().has_on_negative_side(
              alg_approx_seg.target()
            )),
            "The approximated segment should be below initial segment."
          );
        }
      }

      return approximated_segment.supporting_line();
    }

    /* Given three Curve_2 in sequence, update the endpoints to make themselves
     * connected, if they are not already.
     * Keep curr like it is, update the others.
     * Precondition (not checked): the three curved arcs are oriented in the
     * correct way, so that they are directed from the source of prev to the
     * target of next.
     * Precondition (not checked): the source and target of curr lie on prev and
     * on next respectively.
     */
    void update_endpoints(Curve_2& prev, Curve_2& curr, Curve_2& next) const {
      bool one = false, two = false;
      if (prev.target() != curr.source()) {
        prev.set_target(curr.source());
        one = true;
      }
      if (next.source() != curr.target()) {
        next.set_source(curr.target());
        two = true;
      }
      return;
    }

    /* Given a list of Curve_2 and an OutputIterator o, convert those curves to
     * X_monotone_curve_2 and add them to the OutputIterator o */
    template <class OutputIterator>
    OutputIterator convert_and_add_curves(
      std::list<Curve_2>& cvs, OutputIterator o
    ) const {
      /* While iterating, also keep checking if the curves are all connected */
      int piece = 0;
      Alg_point_2 connection = cvs.front().source();
      for (auto& cv : cvs) {
        /* check and update */
        std::cout << piece++ << ": " << cv << '\n';
        CGAL_warning_msg(
          (cv.source() == connection),
          COUT_COLOUR_RED "Warning: " COUT_COLOUR_YELLOW
          "The curve is not connected." COUT_COLOUR_RESET
        );
        connection = cv.target();

        /* convert and add */
        std::vector<X_monotone_curve_2> cv_x_mono_parts;
        make_curve_2_into_many_x_monotone_curve_2(
          cv,
          cv_x_mono_parts
        );
        for (auto& x_cv : cv_x_mono_parts) {
          *o++ = CGAL::make_object(
            Intersection_curve(x_cv, 0)
            //TODO multiplicity? would need to save it in list bisector_parts
          );
        }
      }

      /* return one past the end iterator */
      return o;
    }

    /* Helper function for the construction of the plane bisecotr of two line
     * segments.
     * Given a start point and and end point, iteratively constructs all pieces
     * of the bisector. Other objects passed to the function are the segments
     * themselves, the initial direction where the bisector must continue from
     * start_pt, the iterator where to store the pieces of the bisector as
     * Curve_2 objects, and two different instances of the delimiter lines, that
     * are the lines orthogonal to the segments' endpoints, such that each
     * segment is inside the negative part of his couple of delimiter lines.
     * Returns a one past the end iterator of the list of Curve_2.
     */
    template <class OutputIterator>
    OutputIterator construct_bisector_from_point_to_point(
      Rat_segment_2 s1, Rat_segment_2 s2,
      OutputIterator o,
      Rat_ray_2 start_ray, Rat_ray_2 end_ray,
      Alg_direction_2 curr_direction,
      Alg_delimiter_lines alg_delimiter_lines,
      std::vector<Rat_line_2> delimiter_lines_vector
    ) const {
      /* create converter functors to convert from:
       * - Rational to Algebraic
       * - Algebraic to Cartesian<double>
       * - Cartesian<double> to Rational
       * The last two are used together to convert by approximation from
       * Algebraic to Rational */
      RK_to_AK to_alg;
      AK_to_DK to_dbl;
      DK_to_RK to_rat;

      /* get start and end points, they are the sources of the rays */
      Alg_point_2 start_pt = to_alg(start_ray.source());
      Alg_point_2 end_pt = to_alg(end_ray.source());

      /* determine if the segments intersect and get the intersection point,
       * which is a rational point because the segments are rational (that is,
       * their endpoints are, and their supporting lines have rational
       * coefficients) */
      bool do_intersect = false;
      Rat_point_2 segments_intersection;
      if (CGAL::do_intersect(s1, s2)) {
        do_intersect = true;
        CGAL_assertion_msg(
          CGAL::assign(segments_intersection, CGAL::intersection(s1, s2)),
          "Could not assign segments_intersection."
        );
      }

      /* list to store all Curve_2 bisector part. They all will be converted to
       * X_monotone_curve_2 and inserted into OutputIterator o at the end of
       * this method */
      std::list<Curve_2> bisector_parts;

      /* when computing the straight parts of the bisector between two conic
       * arcs, it is likely that the starting and ending points of these
       * straight parts (segments) are algebraic points. These cannot be added
       * as Curve_2 objects because their supporting conic (a line) could have
       * algebraic coefficients, and this is not supported by the arrangement
       * traits class (Arr_conic_traits_2).
       * As a solution, keep the coefficients of the parabola of the previous
       * arc and wait until the coefficients of the parabola of the next arc are
       * available; then, just use the Curve_2 constructor that does not require
       * the exact endpoints, but the three conic coefficients lists and two
       * approximate endpoints: it computes the actual endpoints as the
       * intersections of the three curves closest to the approximate endpoints.
       * Then, get the computed endpoints of the new curve and update the
       * end point of the previous Curve_2 and the start point of the successive
       * Curve_2
       *
      //  * Keep two parts to approximate for the case in which the segment is the
      //  * part of the bisector that intersects the two segments at their
      //  * intersection (so this happens only when the two segments are
      //  * intersecting, of course).
       */
      bool prev_arc_exists = false;
      Curve_2 prev_arc;
      bool part_to_approximate_exists = false;
      Alg_segment_2 part_to_approximate;
      bool part_to_approximate_1_exists;
      bool part_to_approximate_2_exists;
      Alg_segment_2 part_to_approximate_1;
      Alg_segment_2 part_to_approximate_2;

      /* rename start point */
      Alg_point_2 curr_pt = start_pt;

      /* "walk" through the bisector to find all parts until every piece has
       * been created and added to the OutputIterator o */
      while (curr_pt != end_pt) {
        /* find next intersection with delimiter_lines when going in the
         * direction saved in "curr_direction", then find a middle point
         * between curr_pt and that intersection */
        Alg_point_2 approximate_next_intersection = find_next_intersection(
          curr_direction, curr_pt, delimiter_lines_vector
        );
        Alg_point_2 midpoint = CGAL::midpoint(
          curr_pt,
          approximate_next_intersection
        );

        /* to store the true next intersection and the next direction */
        Alg_point_2 actual_next_intersection;
        Alg_direction_2 next_direction;

        /* determine where this middle point is relative to the two segments
         * s1 and s2, and create the correct piece of the bisector. The
         * objects o1 and o2 that are passed will store in the cases:
         * - PARABOLIC_ARC:       o1 = focus/directrix  o2 = focus/directrix
         * - SUPP_LINE_BISECTOR:  o1 = supp_line1,      o2 = supp_line2
         * - ENDPOINT_BISECTOR:   o1 = endpoint_1,      o2 = endpoint_2   */
        Object o1, o2;
        switch (find_position(midpoint, alg_delimiter_lines, s1, s2, o1, o2)) {

          case PARABOLIC_ARC: {
            std::cout << "PARABOLIC_ARC -> ";
            /* extract directrix and focus */
            Rat_line_2 directrix; Rat_point_2 focus;
            if (CGAL::assign(directrix, o1)) {
              CGAL_assertion(CGAL::assign(focus, o2));
            }
            else {
              CGAL_assertion(CGAL::assign(focus, o1));
              CGAL_assertion(CGAL::assign(directrix, o2));
            }

            /* keep or invert directrix based on curr_direction */
            if (!generally_same_direction(
              to_alg(directrix), curr_direction
            )) {
              directrix = directrix.opposite();
            }

            /* create parabola */
            Parabola supporting_conic(directrix, focus);
            CGAL_assertion(supporting_conic.has_on(curr_pt));

            /* find actual next intersection of parabola */
            actual_next_intersection = supporting_conic.next_intersection(
              curr_pt, delimiter_lines_vector
            );

            /* get tangent with correctly oriented direction */
            Alg_line_2 tangent = supporting_conic.tangent_at_point(
              actual_next_intersection
            );
            next_direction = tangent.direction();

            /* get parabolic arc */
            Curve_2 this_arc = supporting_conic.construct_parabolic_arc(
              curr_pt,
              actual_next_intersection
            );

            /* deal with approximation of previous segment if necessary;
             * distinguish the cases where the segments intersect and the cases
             * in which they do not
             */
            if (
              !do_intersect &&
              part_to_approximate_exists &&
              prev_arc_exists
            ) {
              part_to_approximate_exists = false; // reset flag
              /* get the approximated segment supporting_conic (a line) */
              Rat_line_2 approx_last_segment_line =
                get_approximated_inner_segment_supporting_line(
                part_to_approximate,
                prev_arc,
                this_arc
              );

              /* the prev_arc must be the last one in the list */
              Curve_2 prev_arc_in_list = bisector_parts.back();
              CGAL_assertion_msg(
                same_curves(prev_arc_in_list, prev_arc),
                "prev_arc should be last element in list of parts of bisector."
              );
              bisector_parts.pop_back(); // remove last curve

              /* create new segment curve */
              Curve_2 approx_last_segment_curve(
                0,
               	0,  // supporting conic is a line, so it's linear
               	0,
               	approx_last_segment_line.a(),
               	approx_last_segment_line.b(),
               	approx_last_segment_line.c(),
                CGAL::COLLINEAR,
                part_to_approximate.source(),
               	prev_arc.r(),
               	prev_arc.s(),
               	prev_arc.t(),
               	prev_arc.u(),
               	prev_arc.v(),
               	prev_arc.w(),
                part_to_approximate.target(),
               	this_arc.r(),
               	this_arc.s(),
               	this_arc.t(),
               	this_arc.u(),
               	this_arc.v(),
               	this_arc.w()
              );
              CGAL_assertion_msg(
                approx_last_segment_curve.is_valid(),
                "Created approximated segment curve is not valid"
              );

              /* if needed (because the endpoints might just be the same ones,
               * for example if the approximation was exact), update prev_arc
               * and this_arc end and start point to coincide with start and end
               * of this new approximated segment curve */
              update_endpoints(prev_arc, approx_last_segment_curve, this_arc);

              /* push in list of bisector parts the updated prev_arc and the
               * now approximated segment curve */
              bisector_parts.push_back(prev_arc);
              bisector_parts.push_back(approx_last_segment_curve);
            }
            /* if the segments intersect, we have two parts to approximate, but
             * they are oriented differently and the target of the first and the
             * source of the second are the intersection point of the two
             * segments (saved at the beginning of the method)
             */
            else if (
              do_intersect
              && part_to_approximate_2_exists
              /* part_1 might be first part, in which case there is also no
               * prev_arc */
            ) {
              part_to_approximate_2_exists = false; // reset flag

              /* to save two curves (or just one) */
              Curve_2 part_1, part_2;

              /* we have to approximate part_2 anyway, so let's do it;
               * create approximated part_2 by changing its target */
              Rat_segment_2 approx_seg_part_2 = adjust_endpoint(
                part_to_approximate_2,  // saved before
                this_arc.orientation() == CGAL::COUNTERCLOCKWISE,
                false,  // approximate target of segment
                segments_intersection
              );
              Rat_line_2 supporting_conic_approx_part_2 =
                approx_seg_part_2.supporting_line()
              ;
              part_2 = Curve_2(
                0,
                0,  // supporting conic is a line, so it's linear
                0,
                supporting_conic_approx_part_2.a(),
                supporting_conic_approx_part_2.b(),
                supporting_conic_approx_part_2.c(),
                CGAL::COLLINEAR,
                to_alg(segments_intersection),
                0,
                0,
                0,
                s1.supporting_line().a(),
                s1.supporting_line().b(),   // also s2 would work
                s1.supporting_line().c(),
                actual_next_intersection,
                this_arc.r(),
                this_arc.s(),
                this_arc.t(),
                this_arc.u(),
                this_arc.v(),
                this_arc.w()
              );
              CGAL_assertion_msg(
                part_2.is_valid(),
                "Curve part_2 is not valid."
              );

              /* the source should be actually exact */
              CGAL_assertion_msg(
                part_2.source() == to_alg(segments_intersection),
                "The source of part_2 must be exactly segments_intersection."
              );

              /* update source of this_arc if needed */
              if (this_arc.source() != part_2.target()) {
                this_arc.set_source(part_2.target());
              }

              /* if part_1 was the first internal part of the bisector (i.e.
               * excluding the rays), then we just have to approximate part_2,
               * as we just did;
               * otherwise both need to be approximated */
              if (part_to_approximate_1_exists) {
                part_to_approximate_1_exists = false; // reset flag

                /* there must exist an arc before this one */
                CGAL_assertion_msg(
                  prev_arc_exists,
                  "There must be a prev_arc, to approximate part_1."
                );

                /* the prev_arc must be the last one in the list */
                Curve_2 prev_arc_in_list = bisector_parts.back();
                CGAL_assertion_msg(
                  same_curves(prev_arc_in_list, prev_arc),
                  "prev_arc must be last element in list of parts of bisector."
                );
                bisector_parts.pop_back(); // remove last curve

                /* create approximated part_1 by changing its source */
                Rat_segment_2 approx_seg_part_1 = adjust_endpoint(
                  part_to_approximate_1,  // saved before
                  prev_arc.orientation() == CGAL::COUNTERCLOCKWISE,
                  true,  // approximate source of segment
                  segments_intersection
                );
                Rat_line_2 supporting_conic_approx_part_1 =
                  approx_seg_part_1.supporting_line()
                ;
                part_1 = Curve_2(
                  0,
                  0,  // supporting conic is a line, so it's linear
                  0,
                  supporting_conic_approx_part_1.a(),
                  supporting_conic_approx_part_1.b(),
                  supporting_conic_approx_part_1.c(),
                  CGAL::COLLINEAR,
                  curr_pt,
                  prev_arc.r(),
                  prev_arc.s(),
                  prev_arc.t(),
                  prev_arc.u(),
                  prev_arc.v(),
                  prev_arc.w(),
                  to_alg(segments_intersection),
                  0,
                  0,
                  0,
                  s1.supporting_line().a(),
                  s1.supporting_line().b(),   // also s2 would work
                  s1.supporting_line().c()
                );
                CGAL_assertion_msg(
                  part_1.is_valid(),
                  "Curve part_2 is not valid."
                );

                /* the target should be actually exact */
                CGAL_assertion_msg(
                  part_1.target() == to_alg(segments_intersection),
                  "The target of part_1 must be exactly segments_intersection."
                );

                /* update target of prev_arc if needed */
                if (prev_arc.target() != part_1.source()) {
                  prev_arc.set_target(part_1.source());
                }

                /* push in list of bisector parts the updated prev_arc and the
                 * now approximated part_1 */
                bisector_parts.push_back(prev_arc);
                bisector_parts.push_back(part_1);
              }

              /* push in list of bisector parts the now approximated part_2 */
              bisector_parts.push_back(part_2);
            }

            /* save as Curve_2 in list of bisector parts, save this curve */
            bisector_parts.push_back(this_arc);
            prev_arc = this_arc;
            prev_arc_exists = true;

            break;
          }

          case SUPP_LINE_BISECTOR: {
            std::cout << "SUPP_LINE_BISECTOR -> ";
            /* extract two supporting lines */
            Rat_line_2 supp_line1; Rat_line_2 supp_line2;
            CGAL_assertion(CGAL::assign(supp_line1, o1));
            CGAL_assertion(CGAL::assign(supp_line2, o2));

            /* orient supporting lines according to curr_direction, get
             * bisector, assert that curr_pt is on it */
            if (!generally_same_direction(
              to_alg(supp_line1), curr_direction
            )) {
              supp_line1 = supp_line1.opposite();
            }
            if (!generally_same_direction(
              to_alg(supp_line2), curr_direction
            )) {
              supp_line2 = supp_line2.opposite();
            }
            Alg_line_2 supp_line_bisector = CGAL::bisector(
              to_alg(supp_line1),
              to_alg(supp_line2)
            );
            CGAL_assertion_msg(
              supp_line_bisector.has_on(curr_pt),
              "The point curr_pt should be on the bisector, but it is not."
            );

            /* if the segments to not intersect, just get the next intersection,
             * build the segment and approximate it if needed */
            if (!do_intersect) {
              /* save next_direction, find actual next intersection */
              next_direction = supp_line_bisector.direction();
              actual_next_intersection = find_next_intersection(
                next_direction, curr_pt, delimiter_lines_vector
              );

              /* special cases: if this SUPP_LINE_BISECTOR is the first or the
               * last part of the internal parts of the bisector (that is,
               * excluding the rays), then it can be added easily as it is just
               * an "extension" of that ray, with the same slope (rational).
               * No need to approximate it.
               */
              if (curr_pt == start_pt) {
                Rat_line_2 supporting_conic = start_ray.supporting_line();

                /* direction should be opposite */
                CGAL_assertion_msg(
                  - to_alg(supporting_conic).direction() == next_direction,
                  "Start ray direction should be the opposite as the current."
                );
                supporting_conic = supporting_conic.opposite();

                /* add segment */
                bisector_parts.push_back(Curve_2(
                  0,
                  0,  // supporting conic is a line, so it's linear
                  0,
                  supporting_conic.a(),
                  supporting_conic.b(),
                  supporting_conic.c(),
                  CGAL::COLLINEAR,
                  start_pt,
                  actual_next_intersection
                ));

                /* BREAK here because we already added the bisector part */
                break;
              }
              else if (actual_next_intersection == end_pt) {
                Rat_line_2 supporting_conic = end_ray.supporting_line();

                /* direction should be the same */
                CGAL_assertion_msg(
                  to_alg(supporting_conic).direction() == next_direction,
                  "Start ray direction should be the opposite as the current."
                );

                /* add segment */
                bisector_parts.push_back(Curve_2(
                  0,
                  0,  // supporting conic is a line, so it's linear
                  0,
                  supporting_conic.a(),
                  supporting_conic.b(),
                  supporting_conic.c(),
                  CGAL::COLLINEAR,
                  curr_pt,
                  end_pt
                ));

                /* BREAK here because we already added the bisector part */
                break;
              }

              /* no special case, then just get the segment */
              Alg_segment_2 bisector_part(curr_pt, actual_next_intersection);

              /* save this algebraic segment in the `part_to_approximate`
               * variable; when the next bisector part will be constructed (and
               * it should be an arc) this segment will be approximated and
               * added to the list of pieces of the bisector
               */
              part_to_approximate = bisector_part;
              part_to_approximate_exists = true;

            }
            /* if the segments intersect, we have to create two parts: the first
             * part from curr_pt to the intersection of the two segments, the
             * second part from the intersection of the two segments to a new
             * point to find.
             * At the intersection, the bisector makes a more or less sharp turn
             * to the right (because we give as start and end points the sources
             * of two adjecent rays).
             */
            else {
              /* check that the intersection (saved at the beginning of the
               * function) is on the bisector */
              CGAL_assertion_msg(
                supp_line_bisector.has_on(to_alg(segments_intersection)),
                "The intersection should be on the bisector, but it is not."
              );

              /* before everything we need to find the second bisector of the
               * supporting lines of the two intersecting segments.
               * More precisely, we'd like to have it oriented towards the right
               * relative to the current direction.
               * This is because we have given as end point of this bisector of
               * two intersecting segments the ray that is adjacent to the start
               * ray in counterclockwise direction (the "other" bisector instead
               * connects the other two rays).
               * But why not just connect the two pairs of rays that are
               * opposite to each other? Because that way, the two pieces we are
               * building here on the intersection of the two segments would
               * have a major difference. If the first part has s1 on its right
               * and s2 on its left, then the second part would have the
               * opposite.
               * Instead, by creating the two bisectors so that they "touch" on
               * the intersection of the two segments and then turn away from
               * each other, we ensure that on each bisector, the two pieces we
               * are building here always have s1 on one side and s2 on the
               * other side.
               *
               * To find the new bisector, we invert the supporting line that
               * has a Direction_2 towards the left with respect to the current
               * supp_line_bisector, then we take the bisector of the supporting
               * lines again so to get the new bisector; this works because the
               * CGAL::bisector() function creates the bisector such that the
               * direction is the sum of the normalized directions of the two
               * lines
               */
              Alg_line_2 next_supp_line_bisector;
              if (
                to_alg(supp_line1).direction().counterclockwise_in_between(
                  supp_line_bisector.direction(),
                  - supp_line_bisector.direction()
                )
              ) {
                supp_line1 = supp_line1.opposite();
              }
              else if (
                to_alg(supp_line2).direction().counterclockwise_in_between(
                  supp_line_bisector.direction(),
                  - supp_line_bisector.direction()
                )
              ) {
                supp_line2 = supp_line2.opposite();
              }
              else {
                CGAL_error_msg("The bisector should lie between the two lines");
              }

              next_supp_line_bisector = CGAL::bisector(
                to_alg(supp_line1), to_alg(supp_line2)
              );
              CGAL_assertion_msg(
                next_supp_line_bisector.has_on(to_alg(segments_intersection)),
                "The new bisector should contain the intersection point."
              );
              next_direction = next_supp_line_bisector.direction();

              /* find actual next intersection according to the new direction */
              actual_next_intersection = find_next_intersection(
                next_direction,
                to_alg(segments_intersection),
                delimiter_lines_vector
              );

              /* special cases:
               * the first part, or the second part, or both, are the first
               * and/or last parts of the bisector.
               * if this is the first bisector part (meaning that before this
               * piece there was only the ray) then also the start point is
               * already known
               */
              bool first_to_approximate = true, second_to_approximate = true;
              if (curr_pt == start_pt) {
                first_to_approximate = false;
                Rat_line_2 supporting_conic =
                  start_ray.supporting_line().opposite()
                ;
                bisector_parts.push_back(Curve_2(
                  0,
                  0,  // supporting conic is a line, so it's linear
                  0,
                  supporting_conic.a(),
                  supporting_conic.b(),
                  supporting_conic.c(),
                  CGAL::COLLINEAR,
                  curr_pt,
                  to_alg(segments_intersection)
                ));
              }
              if (actual_next_intersection == end_pt) {
                second_to_approximate = false;
                Rat_line_2 supporting_conic = end_ray.supporting_line();
                bisector_parts.push_back(Curve_2(
                  0,
                  0,  // supporting conic is a line, so it's linear
                  0,
                  supporting_conic.a(),
                  supporting_conic.b(),
                  supporting_conic.c(),
                  CGAL::COLLINEAR,
                  to_alg(segments_intersection),
                  actual_next_intersection
                ));
              }

              /* if the parts were the first and the last parts of the internal
               * parts of the bisector then we don't have to approximate
               * anything;
               * if the second part was the last internal part of the bisector
               * but the first one was not, we have to deal with the
               * approximation here, as there will be no more PARABOLIC_ARC
               * after this where we can deal with the approximation;
               * in the opposite case, where the first part was the first part
               * of the internal parts of the bisector but the second was not,
               * we save just the second part and leave it to be approximated
               * when the next PARABOLIC_ARC is found;
               * if both parts need to be approximated, we can just save them
               * and let the next PARABOLIC_ARC case deal with the approximation
               */
              if (first_to_approximate && second_to_approximate) {
                part_to_approximate_1 = Alg_segment_2(
                  curr_pt, to_alg(segments_intersection)
                );
                part_to_approximate_2 = Alg_segment_2(
                  to_alg(segments_intersection), actual_next_intersection
                );
                part_to_approximate_1_exists = true;
                part_to_approximate_2_exists = true;
              }
              /* in this case deal with approximation of first part already
               * here, since after this the while loop ends */
              else if (first_to_approximate && !second_to_approximate) {

                /* extract part_2 that was already added (extract because we
                 * have to add the parts in order)
                 * and prev_arc (because its target needs to be updated) */
                Curve_2 part_2 = bisector_parts.back();
                bisector_parts.pop_back();
                CGAL_assertion_msg(
                  prev_arc_exists,
                  "prev_arc should exist in the list, otherwise we cannot "
                  "approximate part_1."
                );
                Curve_2 prev_arc_in_list = bisector_parts.back();
                bisector_parts.pop_back();
                CGAL_assertion_msg(
                  same_curves(prev_arc_in_list, prev_arc),
                  "prev_arc must be last element in list of parts of bisector."
                );

                /* create approximated part_1 by adjusting its source */
                Alg_segment_2 seg_part_1(
                  curr_pt,
                  to_alg(segments_intersection)
                );
                Rat_segment_2 approx_seg_part_1 = adjust_endpoint(
                  seg_part_1,
                  prev_arc.orientation() == CGAL::COUNTERCLOCKWISE,
                  true,   // approximate source of segment
                  segments_intersection
                );
                Rat_line_2 supporting_conic_approx_part_1 =
                  approx_seg_part_1.supporting_line()
                ;
                Curve_2 part_1(
                  0,
                  0,  // supporting conic is a line, so it's linear
                  0,
                  supporting_conic_approx_part_1.a(),
                  supporting_conic_approx_part_1.b(),
                  supporting_conic_approx_part_1.c(),
                  CGAL::COLLINEAR,
                  curr_pt,
                  prev_arc.r(),
                  prev_arc.s(),
                  prev_arc.t(),
                  prev_arc.u(),
                  prev_arc.v(),
                  prev_arc.w(),
                  to_alg(segments_intersection),
                  0,
                  0,
                  0,
                  s1.supporting_line().a(),
                  s1.supporting_line().b(),   // also s2 would work
                  s1.supporting_line().c()
                );
                CGAL_assertion_msg(
                  part_1.is_valid(),
                  "Curve part_2 is not valid."
                );

                /* the target should be actually exact */
                CGAL_assertion_msg(
                  part_1.target() == to_alg(segments_intersection),
                  "The target of part_1 must be exactly segments_intersection."
                );

                /* update target of prev_arc if needed */
                if (prev_arc.target() != part_1.source()) {
                  prev_arc.set_target(part_1.source());
                }

                /* push the updated prev_arc, then the approximated first part,
                 * then the second part (which is the last internal part of the
                 * segments bisector) */
                bisector_parts.push_back(prev_arc);
                bisector_parts.push_back(part_1);
                bisector_parts.push_back(part_2);
              }
              /* in this case leave the approximation of the second part to when
               * the next PARABOLIC_ARC is found */
              else if (!first_to_approximate && second_to_approximate) {
                part_to_approximate_2 = Alg_segment_2(
                  to_alg(segments_intersection), actual_next_intersection
                );
                part_to_approximate_2_exists = true;
              }
              /* both parts were already added, meaning that they were the first
               * and last internal parts of the bisector (that is, excluding the
               * unbounded rays) */
              else break; //nothing to do
            }

            break;
          }

          case ENDPOINT_BISECTOR: {
            std::cout << "ENDPOINT_BISECTOR -> ";
            /* extract two endpoints */
            Rat_point_2 endpoint1; Rat_point_2 endpoint2;
            CGAL_assertion(CGAL::assign(endpoint1, o1));
            CGAL_assertion(CGAL::assign(endpoint2, o2));

            /* create bisector, orient it according to curr_direction */
            Rat_line_2 endpoint_bisector = CGAL::bisector(endpoint1, endpoint2);
            if (!generally_same_direction(
              to_alg(endpoint_bisector), curr_direction
            )) {
              endpoint_bisector = endpoint_bisector.opposite();
            }

            /* assert curr_pt is on bisector */
            CGAL_assertion(to_alg(endpoint_bisector).has_on(curr_pt));

            /* save next_direction, find actual next intersection */
            next_direction = to_alg(endpoint_bisector).direction();
            actual_next_intersection = find_next_intersection(
              next_direction, curr_pt, delimiter_lines_vector
            );

            /* check that there is no SUPP_LINE_BISECTOR to approximate */
            CGAL_assertion_msg(
              !part_to_approximate_exists,
              "There cannot be a SUPP_LINE_BISECTOR followed by an "
              "ENDPOINT_BISECTOR that is internal (i.e. not a ray)."
            );
            if (prev_arc_exists) prev_arc_exists = false;

            /* no need to approximate the segment in this case: the supporting
             * conic (a line) has rational coefficients; just save curve with
             * rational coefficients and the two alg points as endpoints,
             * orientation COLLINEAR.
             * save as Curve_2 in list of bisector parts */
            bisector_parts.push_back(Curve_2(
              0,
              0,  // supporting conic is a line, so it's linear
              0,
              endpoint_bisector.a(),
              endpoint_bisector.b(),
              endpoint_bisector.c(),
              CGAL::COLLINEAR,
              curr_pt,
              actual_next_intersection
            ));

            break;
          }

          default: break; // should never happen
        }

        /* update current starting point and current direction of the next piece
         * of the bisector */
        curr_pt = actual_next_intersection;
        curr_direction = next_direction;
      }
      std::cout << "FINISHED\n";

      /* add all Curve_2 to OutputIterator o (conversion to X_monotone_curve_2
       * and adding to the actual list of projected intersections happens in the
       * main function) */
      for (auto cv : bisector_parts) *o++ = cv;

      /* return one past the end iterator */
      return o;
    }
  }; // end of class Construct_projected_intersections_2

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const {
    return Construct_projected_intersections_2();
  }

/* ########################################################################## */
/* ###               END: COMPUTE BISECTOR OF TWO SEGMENTS                ### */
/* ########################################################################## */



  class Compare_z_at_xy_3 {
  public:
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      printf("---> Compare at point\n");
      return CGAL::compare(sqdistance(p, s1), sqdistance(p, s2));
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      printf("---> Compare at cv\n");
      /* compare using the middle point */
      Point_2 p = construct_middle_point(cv);
      return this->operator()(p, s1, s2);
    }

    Comparison_result operator()(const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      printf("\n ---> Compare not intersecting\n");
      /* if the two unbounded surfaces do not intersect, then they must
       * represent the same segment's distance function */
      CGAL_assertion_msg(
        (s1 == s2 || s1 == Rat_segment_2(s2.target(), s2.source())),
        "Distance function surfaces do not intersect but they are not the same"
      );
      return CGAL::EQUAL; // they are literally the same surface
    }
  };

  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  {
    printf("\n#################################\nCreated Compare_z_at_xy_3 obj#################################\n");
    return Compare_z_at_xy_3();
  }


  /* Call helper function with flag set to true */
  class Compare_z_at_xy_above_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      return compare_z_at_xy_3_helper(cv, s1, s2, true);
    }

  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3();
  }


  /* Call helper function with flag set to false */
  class Compare_z_at_xy_below_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      return compare_z_at_xy_3_helper(cv, s1, s2, false);
    }
  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const {
    return Compare_z_at_xy_below_3();
  }

  /* Helper function for Compare_z_at_xy_above_3 and Compare_z_at_xy_below_3 */
  static Comparison_result compare_z_at_xy_3_helper(
    const X_monotone_curve_2& cv,
    const Xy_monotone_surface_3& s1,
    const Xy_monotone_surface_3& s2,
    bool compare_above
  ) {
    Algebraic move_by = 1;

    /* construct a point on the curve cv, assert equidistant from s1 and s2 */
    Alg_point_2 midpoint = construct_middle_point(cv);
    Algebraic difference = sqdistance(midpoint, s1) - sqdistance(midpoint, s2);

    std::cout << std::endl << "##############################" << std::endl;
    std::cout << "# Compare s1[" << s1 << "] and s2[" << s2 << "] "
              << (compare_above ? "above" : "below") << " cv=[ "
              << cv.r() << "x^2 + "
              << cv.s() << "y^2 + "
              << cv.t() << "xy + "
              << cv.u() << "x + "
              << cv.v() << "y + "
              << cv.w() << std::endl
    ;

    std::cout << "# midpoint(" << midpoint << ") , ";

    /* print warning if necessary */
    char message[100];
    sprintf(
      message,
      COUT_COLOUR_RED "Warning: " COUT_COLOUR_YELLOW
      "s1 and s2 are not equidistant, difference: %lf" COUT_COLOUR_RESET,
      CGAL::to_double(difference)
    );
    CGAL_warning_msg((difference == 0), message);
    std::cout << "the difference is: " << difference.toString() << ".";

    /* get converter and convert */
    RK_to_AK to_alg;
    Alg_segment_2 alg_s1 = to_alg(s1);
    Alg_segment_2 alg_s2 = to_alg(s2);

    Alg_point_2 moved_point;
    Algebraic displacement = compare_above ? move_by : -move_by;
    if (cv.is_vertical()) {
      std::cout << " -- CV IS VERTICAL --";
      moved_point = Alg_point_2(midpoint.x() - displacement, midpoint.y());
    }
    // else if (more_vertical_than_horizontal(cv)) {
    //   moved_point = Alg_point_2(midpoint.x() - displacement, midpoint.y());
    // }
    else {
      moved_point = Alg_point_2(midpoint.x(), midpoint.y() + displacement);
    }

    std::cout << " Moved midpoint to pt(" << moved_point << ")" << std::endl;

    if (sqdistance(moved_point, s1) < sqdistance(moved_point, s2)) {
      std::cout << "# Returning CGAL::SMALLER." << std::endl
                << "##############################" << std::endl
      ;
      return CGAL::SMALLER;
    }
    else {
      std::cout << "# Returning CGAL::LARGER." << std::endl
                << "##############################" << std::endl
      ;
      return CGAL::LARGER;
    }

    CGAL_error_msg("This function is not working properly");
    return CGAL::EQUAL;
  }

  // static more_vertical_than_horizontal(X_monotone_curve_2 cv) {
  //   return
  //     CGAL::abs(cv.source().x() - cv.target().x())  // dx
  //     <
  //     CGAL::abs(cv.source().y() - cv.target().y())  // dy
  //   ;
  // }

}; // class L2_segment_voronoi_traits_2
} // namespace CGAL

#endif // CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H
