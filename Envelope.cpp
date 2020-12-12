#include <iostream>
#include <list>
#include <string>

// envelope
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>
//
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Cartesian.h>

using namespace std;

namespace envelope
{
    typedef CGAL::Cartesian<double> Kernel;
    typedef Kernel::Point_2 Point_2;
    typedef Kernel::Segment_2 Segment_2;
    typedef CGAL::Arr_segment_traits_2<Kernel> Segment_traits_2;
    typedef Segment_traits_2::Curve_2 Curve_2;
    typedef CGAL::Envelope_diagram_1<Segment_traits_2> Diagram_1;

    Diagram_1::Edge_const_handle e;
    Diagram_1::Vertex_const_handle v;

    // If envl_type := True Find upper envelope
    // If envl_type := False Find lower envelope
    void find_envelope(list<Segment_2> &Quadrant, int id, bool envl_type, list<Point_2> &envelope)
    {
        std::list<Curve_2> segment_list;
        Diagram_1 min_diag;

        for (auto it = Quadrant.begin(); it != Quadrant.end(); it++)
        {
            if (it->source() == it->target())
            {
                envelope.push_back(it->source());
                continue;
            }
            segment_list.push_back(Curve_2(it->source(),it->target()));
        }
        // COMPUTING THE RESPECTIVE ENVELOPE
        if (envl_type)
        {
            upper_envelope_x_monotone_2(segment_list.begin(), segment_list.end(), min_diag);
            cerr << "COMPUTED UPPER ENVELOPE " << id << endl;
        }
        else
        {
            lower_envelope_x_monotone_2(segment_list.begin(), segment_list.end(), min_diag);
            cerr << "COMPUTED LOWER ENVELOPE " << id << endl;
        }
        e = min_diag.leftmost();
        // OUTPUTING THE ENVELOPE POINTS
        while (e != min_diag.rightmost())
        {
            v = e->right();
            auto pt = v->point();
            envelope.push_back(pt);
            e = v->right();
        }
        min_diag.clear();
        segment_list.clear();
    }
} // namespace envelope

int main()
{
    using namespace envelope;
    // Consrtuct the input segments and label them 'A' ... 'H'.
    std::list<Segment_2> segments;
    segments.push_back(Segment_2(Point_2(0, 1), Point_2(2, 3)));
    segments.push_back(Segment_2(Point_2(1, 2), Point_2(4, 5)));
    segments.push_back(Segment_2(Point_2(1, 5), Point_2(7, 2)));
    segments.push_back(Segment_2(Point_2(4, 2), Point_2(6, 4)));
    segments.push_back(Segment_2(Point_2(8, 3), Point_2(8, 6)));
    segments.push_back(Segment_2(Point_2(9, 2), Point_2(12, 4)));
    segments.push_back(Segment_2(Point_2(10, 2), Point_2(12, 1)));
    segments.push_back(Segment_2(Point_2(11, 0), Point_2(11, 5)));

    list<Point_2> envelope;

    find_envelope(segments, 0, false, envelope);
    for (auto pt : envelope)
    {
        cout << pt << endl;
    }
}