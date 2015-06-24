#ifndef _PRIMITIVES_H_
#define _PRIMITIVES_H_

#include <vector>

#include "random.h"
#include "convex_polygon.h"

template <class Kernel, class Vbb>
class My_vertex_base : public Vbb
{
public:
    typedef CConvexPolygon<Kernel> ConvexPolygon;
    
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Weight   Weight;
    typedef typename Kernel::Point_2  Point;
    typedef typename Kernel::Vector_2 Vector;
    typedef typename Kernel::Weighted_point_2 Weighted_point;
    
    typedef typename Vbb::Triangulation_data_structure TDS;
    typedef typename TDS::Face_handle   Face_handle;
    typedef typename TDS::Vertex_handle Vertex_handle;

    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename Vbb::template Rebind_TDS<TDS2>::Other Vb2;
        typedef My_vertex_base<Kernel, Vb2> Other;
    };
        
private:
    int m_index;
    ConvexPolygon m_cell;
    
public:
    My_vertex_base() : Vbb() 
    { 
        m_index = -1;
    }
   
    My_vertex_base(const Weighted_point& p) : Vbb(p)
    {
        m_index = -1;
    }

    My_vertex_base(const Weighted_point& p, Face_handle f) : Vbb(p, f)
    {
        m_index = -1;
    }

    ~My_vertex_base()
    {
    }
    
    void set_index(const int x) { m_index = x; }
    const int get_index() const { return m_index; }

    // POSITION / WEIGHT //
    
    const Point& get_position() const 
    { 
        return this->point().point(); 
    }

    void set_position(const Point& p)
    {
        FT w = get_weight();
        Weighted_point wp(p, w);
        this->set_point(wp);
    }
    
    const FT get_weight() const 
    { 
        return this->point().weight(); 
    }
    
    void set_weight(const FT w)
    {
        const Point& p = get_position();
        Weighted_point wp(p, w);
        this->set_point(wp);
    }
    
    // CELL //

    void reset_cell()
    {
        m_cell.clear();
    }
    
    void set_cell(const ConvexPolygon& shape)
    {
        m_cell = shape;
    }
        
    // AREA ; VARIANCE ; CENTROID //
    
    FT compute_area() const
    {
        return m_cell.compute_area();
    }
    
    Point compute_centroid() const
    {
        return m_cell.compute_centroid();
    }
    
    FT compute_variance() const
    {
        return m_cell.compute_variance(get_position());
    }
};

template <class Kernel, class Fbb>
class My_face_base : public Fbb
{
public:
    typedef typename Fbb::Triangulation_data_structure  TDS;
    typedef typename TDS::Vertex_handle Vertex_handle;
    typedef typename TDS::Face_handle   Face_handle;
    
    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename Fbb::template Rebind_TDS<TDS2>::Other Fb2;
        typedef My_face_base<Kernel, Fb2> Other;
    };
    
public:
    My_face_base() 
    : Fbb() 
    { 
    }
    
    My_face_base(Vertex_handle v1,
                 Vertex_handle v2,
                 Vertex_handle v3)
    : Fbb(v1, v2, v3)
    {
    }
    
    My_face_base(Vertex_handle v1,
                 Vertex_handle v2,
                 Vertex_handle v3,
                 Face_handle f1,
                 Face_handle f2,
                 Face_handle f3)
    : Fbb(v1,v2,v3,f1,f2,f3)
    {
    }
};

#endif // _PRIMITIVES_H_
