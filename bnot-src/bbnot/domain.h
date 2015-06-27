#ifndef _DOMAIN_H_
#define _DOMAIN_H_

// STL
#include <vector>
#include <string>

// CGAL
#include <CGAL/intersections.h>

// Qt
#include <QtOpenGL>

// local
#include "util.h"

template <class Kernel>
class CDomain
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
    std::vector<Point>  m_points; // boundary
    FT m_width, m_height; // bbox
    FT m_area;
    
public:
    CDomain() 
    { 
        clear();
    }
    
    void clear()
    {
        m_points.clear();
        m_height = 0.0;
        m_width = 0.0;
        m_area = 0.0;
    }
    
    // GET / SET //

    Point origin() const { return Point(0.0, 0.0); }
    
    const Point& get_point(unsigned i) const { return m_points[i]; }    
        
    const FT width() const { return m_width; }
    
    const FT height() const { return m_height; }
    
    // BOUNDARY //
    
    Segment get_segment(unsigned i) const
    {
        unsigned j = (i + 1) % m_points.size();
        const Point& a = m_points[i];
        const Point& b = m_points[j];
        return Segment(a, b);
    }    
    
    template <class OutputIterator>
    void add_segments(OutputIterator out)
    {
        for (unsigned i = 0; i < m_points.size(); ++i)
        {        
            Segment segment = get_segment(i);
            *out++ = segment;
        }
    }

    // INIT //
    
    void init_polygon(const FT radius, 
                      const unsigned nb)
    {
        clear();
        m_width  = radius;
        m_height = radius;
        generate_regular_polygon(nb, radius, Point(0.0, 0.0), 
                                 std::back_inserter(m_points));
    }

    void init_rectangle(const FT width, 
                        const FT height)
    {
        clear();
        m_width  = width;
        m_height = height;
        m_points.push_back(Point(-width, -height));
        m_points.push_back(Point( width, -height));
        m_points.push_back(Point( width,  height));
        m_points.push_back(Point(-width,  height));
    }

    void init_area()
    {
        m_area = compute_polygon_area();
    }
    
    // AREA //

    FT get_area() const { return m_area; }

    FT compute_polygon_area() const
    {
        // rect
        return 4.0*width()*height();
    }
    
    // INSIDE / OUTSIDE //
    
    bool is_inside(const Point& p) const
    {
        if (std::abs(p.x()) > width ()) return false;
        if (std::abs(p.y()) > height()) return false;
        return true;
    }
    
    bool is_outside(const Point& p) const
    {
        return !is_inside(p);
    }
    
    bool is_corner(const Point& p) const
    {
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            const Point& pi = get_point(i);
            if (p == pi) return true;
        }
        return false;
    }
    
    // CLAMP //
    
    Point clamp(const Point& p) const
    {
        // rect
        FT x = p.x();
        if (x >  width()) x =  width();
        if (x < -width()) x = -width();
        
        FT y = p.y();
        if (y >  height()) y =  height();
        if (y < -height()) y = -height();

        return Point(x, y);
    }
    
    Segment clamp(const Segment& segment) const
    {
        Point ps = segment.source();
        Point pt = segment.target();
        
        bool source_inside = is_inside(ps);
        bool target_inside = is_inside(pt);

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
        FT dab = std::sqrt(CGAL::squared_distance(a, b));
        Vector n = (b - a) / dab;
        FT daq = n * (query - a);
        Point q = a + daq*n;
        return q;
    }
    
    unsigned find_nearest_line(const Point& p) const
    {
        // rect
        FT x = p.x();
        FT xpos = std::abs(x - width());
        FT xneg = std::abs(x + width());
        FT xbest = std::min(xpos, xneg);
        
        FT y = p.y();
        FT ypos = std::abs(y - height());
        FT yneg = std::abs(y + height());
        FT ybest = std::min(ypos, yneg);
        
        if (xbest < ybest)
        {
            if (xpos < xneg) return 1;
            return 3;
        }
        
        if (ypos < yneg) return 2;
        return 0;
    }
    
    int find_side(const Point& a) const
    {
        int side = -1;
        if (std::abs(a.x() - width ()) < EPS) side = 1;
        if (std::abs(a.x() + width ()) < EPS) side = 3;
        if (std::abs(a.y() - height()) < EPS) side = 2;
        if (std::abs(a.y() + height()) < EPS) side = 0;
        return side;
    }
    
    // assume that an intersection exists
    // since the polygon is convex, the intersection is unique
    Point intersect(const Segment& segment) const
    {
        for (unsigned i = 0; i < m_points.size(); ++i)
        {
            Segment si = get_segment(i);
            CGAL::Object result = CGAL::intersection(segment, si);
            const Point* iq = CGAL::object_cast<Point>(&result);
            if (iq) return *iq;
        }
        // never come here
        std::cout << red << "intersection not found" << white << std::endl;
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
};

#endif
