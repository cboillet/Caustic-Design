#ifndef _PIXEL_H_
#define _PIXEL_H_

// CGAL
#include <CGAL/ch_graham_andrew.h>

// local
#include "convex_polygon.h"

template <class Kernel> 
class CPixel
{
public:
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_2    Point;
    typedef typename Kernel::Vector_2   Vector;
    typedef typename Kernel::Ray_2      Ray;
    typedef typename Kernel::Line_2     Line;
    typedef typename Kernel::Segment_2  Segment;
    typedef typename Kernel::Triangle_2 Triangle;
    typedef CConvexPolygon<Kernel> ConvexPolygon;
    
private:
    FT m_value; // > 0.0
    ConvexPolygon m_shape;
    
public:
    CPixel(const FT value = 1.0)
    {
        m_value = value;
    }

    CPixel(const ConvexPolygon& shape, const FT value = 1.0)
    {
        m_value = value;
        m_shape = shape;
    }
    
    CPixel(const CPixel& pixel)
    {
        m_value = pixel.m_value;
        m_shape = pixel.m_shape;
    }
    
    CPixel& operator = (const CPixel& pixel)
    {
        m_value = pixel.m_value;
        m_shape = pixel.m_shape;
        return *this;
    }
    
    const FT get_value() const { return m_value; }
    void set_value(const FT value) { m_value = value; }
    
    const ConvexPolygon& get_shape() const { return m_shape; }    
    void set_shape(const ConvexPolygon& shape) { m_shape = shape; } 
    
    unsigned nb_corners() const { return m_shape.nb_points(); }
    const Point& get_corner(unsigned i) const { return m_shape.get_point(i); }
    
    FT compute_area() const
    {
        FT area = m_shape.compute_area();
        return m_value * area;
    }
    
    Point compute_centroid() const
    {
        return m_shape.compute_centroid();
    }
    
    FT compute_variance(const Point& p) const
    {
        FT variance = m_shape.compute_variance(p);
        return m_value * variance;
    }
    
    template <class InputIterator> // value_type = Point
    void compute_shape(InputIterator first, InputIterator last)
    {
        m_shape.clear();
        CGAL::ch_graham_andrew(first, last, std::back_inserter(m_shape.points()));
    }
};

#endif
