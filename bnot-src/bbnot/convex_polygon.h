#ifndef _CONVEX_POLYGON_H_
#define _CONVEX_POLYGON_H_

// STL
#include <vector>

// CGAL
#include <CGAL/intersections.h>

// Qt
#include <QtOpenGL>

// local
#include "util.h"

template <class Kernel>
class CConvexPolygon
{
public:
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_2    Point;
    typedef typename Kernel::Vector_2   Vector;
    typedef typename Kernel::Ray_2      Ray;
    typedef typename Kernel::Line_2     Line;
    typedef typename Kernel::Segment_2  Segment;
    typedef typename Kernel::Triangle_2 Triangle;
    
private:
    std::vector<Point> m_points;
    
public:
    CConvexPolygon() { }

    CConvexPolygon(const CConvexPolygon& rhs)
    {
        m_points = rhs.m_points;
    }

    CConvexPolygon& operator = (const CConvexPolygon& rhs)
    {
        m_points = rhs.m_points;
        return *this;
    }
    
    void clear()
    {
        m_points.clear();
    }
    
    // GET //

    std::vector<Point>& points() { return m_points; }

    unsigned nb_points() const { return m_points.size();  }
    
    const Point& get_point(unsigned i) const { return m_points[i]; }
    
    const std::vector<Point>& get_points() const { return m_points; }
    
    Segment get_segment(unsigned i) const
    {
        unsigned j = (i + 1) % m_points.size();
        const Point& a = m_points[i];
        const Point& b = m_points[j];
        return Segment(a, b);
    }    
    
    template <class OutputIterator>
    void collect_segments(OutputIterator out) const
    {
        for (unsigned i = 0; i < m_points.size(); ++i)
        {        
            Segment segment = get_segment(i);
            *out++ = segment;
        }
    }
    
    // APPEND //
    
    void append_point(const Point& p)
    {
        m_points.push_back(p);
    }
    
    template <class InputIterator>
    void append(InputIterator first, InputIterator last)
    {
        m_points.insert(m_points.end(), first, last);
    }
    
    // INIT //
    
    void init_regular_polygon(unsigned nb,
                              FT radius, 
                              Point center = Point(0.0, 0.0))
    {
        clear();
        generate_regular_polygon(nb, radius, center, 
                                 std::back_inserter(m_points));
    }

    void init_rectangle(FT width, 
                        FT height,
                        Point center = Point(0.0, 0.0))
    {
        clear();
        m_points.push_back(Point(center.x() - width, center.y() - height));
        m_points.push_back(Point(center.x() + width, center.y() - height));
        m_points.push_back(Point(center.x() + width, center.y() + height));
        m_points.push_back(Point(center.x() - width, center.y() + height));
    }

    // AREA //

    FT compute_area() const
    {
        if (m_points.size() < 3) return 0.0;
        
        FT area = 0.0;
        const Point& p0 = get_point(0);
        for (unsigned i = 2; i < m_points.size(); ++i)
        {
            Segment segment = get_segment(i-1);
            Triangle triangle(p0, segment.source(), segment.target());
            area += triangle.area();
        }
        return area;
    }

    // CENTROID //

    Point compute_centroid() const
    {
        if (m_points.empty()) return Point();
        if (m_points.size() == 1) return get_point(0);
        if (m_points.size() == 2) return CGAL::midpoint(get_point(0), 
                                                        get_point(1));
        
        FT sum_area = 0.0;
        Vector sum_vector = CGAL::NULL_VECTOR;
        const Point& p0 = get_point(0);
        for (unsigned i = 2; i < m_points.size(); ++i)
        {
            Segment segment = get_segment(i-1);
            Triangle triangle(p0, segment.source(), segment.target());
            
            FT area = triangle.area();
            Point barycenter = CGAL::centroid(triangle);
            
            sum_area += area;
            sum_vector = sum_vector + area*(barycenter - CGAL::ORIGIN);
        }
        if (sum_area == 0.0) return get_point(0);
        return CGAL::ORIGIN + (sum_vector / sum_area);
    }

    // VARIANCE //
    
    FT compute_variance(const Point& q) const
    {
        FT variance = 0.0;
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            Segment segment = get_segment(i);
            variance += compute_variance_per_segment(q, segment);
        }
        return variance;
    }
    
    FT compute_variance_per_segment(const Point& x, const Segment& ab) const
    {
        // \int_(x, a, b) |y - x|^2 dy
        // (x, a, b) can be either CCW or CW
        // return signed variance
        
        Point a = ab.source();
        Point b = ab.target();
        Point q = compute_orthogonal_projection(x, ab);
                
        Triangle xqa(x, q, a);
        FT E0 = compute_variance_per_right_triangle(x, q, a); // >=0
        //if (xqa.orientation() == CGAL::COUNTERCLOCKWISE)
        if (xqa.area() >= 0.0)
            E0 = - E0;
        
        Triangle xqb(x, q, b);
        FT E1 = compute_variance_per_right_triangle(x, q, b); // >=0
        //if (xqb.orientation() == CGAL::CLOCKWISE)
        if (xqb.area() < 0.0)
            E1 = - E1;
        
        return (E0 + E1);
    }
    
    FT compute_variance_per_right_triangle(const Point& x, const Point& a, const Point& b) const
    {
        // E = \int_(x, a, b) |y - x|^2 dy
        // where (x, a, b) is 90o at a
        // E >= 0
        FT base   = std::sqrt(CGAL::squared_distance(a, b));
        FT height = std::sqrt(CGAL::squared_distance(a, x));
        return (height*height*height*base)/4 + (height*base*base*base)/12;
    }

    // INSIDE //
    
    bool is_outside(const Point& p) const
    {
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            Line line = get_segment(i).supporting_line();
            if (line.has_on_negative_side(p))
                return true;
        }
        return false;
    }
    
    // CLAMP //
    
    Point clamp(const Point& p) const
    {
        FT min_dist2 = 1e100;
        Point clamped_point = p;
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            Segment segment = get_segment(i);
            Line line = segment.supporting_line();
            if (!line.has_on_negative_side(p)) continue;
                
            Point q = project_to(p, segment);
            FT dist2 = CGAL::squared_distance(p, q);
            if (dist2 < min_dist2)
            {
                clamped_point = q;
                min_dist2 = dist2;
            }
        }
        return clamped_point;
    }
    
    Segment clamp(const Segment& segment) const
    {
        Point ps = segment.source();
        Point pt = segment.target();
        
        bool source_inside = !is_outside(ps);
        bool target_inside = !is_outside(pt);

        if (source_inside && target_inside)
            return segment;

        if (!source_inside && !target_inside)
            return clamp_if_possible(segment);
        
        Point q = intersect(segment);
        if (source_inside) pt = q;
        if (target_inside) ps = q;
        return Segment(ps, pt);
    }
    
    Segment clamp_if_possible(const Segment& segment) const
    {
        // both ps and pt are outside
        Point ps = segment.source();
        Point pt = segment.target();
        
        unsigned ns = find_nearest_line(ps);
        CGAL::Object rs  = CGAL::intersection(segment, get_segment(ns));
        const Point* iqs = CGAL::object_cast<Point>(&rs);
        if (!iqs) return Segment(ps, ps);

        unsigned nt = find_nearest_line(pt);
        CGAL::Object rt  = CGAL::intersection(segment, get_segment(nt));
        const Point* iqt = CGAL::object_cast<Point>(&rt);
        if (!iqt) return Segment(pt, pt);
        
        return Segment(*iqs, *iqt);
    }
    
    // PROJECT //
    
    Point project(const Point& p) const
    {
        FT min_dist2 = 1e100;
        Point projected_point = p;
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            Segment segment = get_segment(i);
            Point q = project_to(p, segment);
            FT dist2 = CGAL::squared_distance(p, q);
            if (dist2 < min_dist2)
            {
                projected_point = q;
                min_dist2 = dist2;
            }
        }
        return projected_point;
    }
    
    Point project_to(const Point& query, const Segment& segment) const
    {
        Point proj = compute_orthogonal_projection(query, segment);
        if (segment.has_on(proj)) return proj;
        
        const Point& ps = segment.source();
        const Point& pt = segment.target();
        FT ds2 = CGAL::squared_distance(proj, ps);
        FT dt2 = CGAL::squared_distance(proj, pt);
        if (ds2 < dt2) return ps;
        return pt;
    }
    
    Point compute_orthogonal_projection(const Point& query, const Segment& segment) const
    {
        const Point& a = segment.source();
        const Point& b = segment.target();
        if (a == b) return a;
        
        FT dab = std::sqrt(CGAL::squared_distance(a, b));
        Vector n = (b - a) / dab;
        FT daq = n * (query - a);
        Point q = a + daq*n;
        
        // debug
        /*
        Line line = segment.supporting_line();
        Point proj = line.projection(query);
        FT diff = std::sqrt(CGAL::squared_distance(proj, q));
        if (diff > 1e-10) 
            std::cout << "Projection: mine " << q << " ; CGAL " << proj << std::endl;
        */
        //
        
        return q;
    }
    
    unsigned find_nearest_line(const Point& p) const
    {
        unsigned nearest = 0;
        FT min_dist2 = 1e100;
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            Segment segment = get_segment(i);
            Point q = compute_orthogonal_projection(p, segment);
            FT dist2 = CGAL::squared_distance(p, q);
            if (dist2 < min_dist2)
            {
                nearest = i;
                min_dist2 = dist2;
            }
        }
        return nearest;
    }
    
    // INTERSECT //
    
    Point intersect(const Segment& segment) const
    {
        // assume that segment has one endpoint inside and the other outside
        // since the polygon is convex, the intersection is unique
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            Segment si = get_segment(i);
            CGAL::Object result = CGAL::intersection(segment, si);
            const Point* iq = CGAL::object_cast<Point>(&result);
            if (iq) return *iq;
        }
        // never come here
        return Point(0.0, 0.0);
    }
    
    // DRAW //
    
    void draw_boundary() const
    {
        glBegin(GL_LINE_LOOP);
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            const Point& pi = get_point(i);
            glVertex2d(pi.x(), pi.y());
        }
        glEnd();
    }    

    void draw_filled() const
    {
        glBegin(GL_POLYGON);
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            const Point& pi = get_point(i);
            glVertex2d(pi.x(), pi.y());
        }
        glEnd();
    }
};

#endif
