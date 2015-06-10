#ifndef _RT2_H_
#define _RT2_H_ 1

// local
#include "convex_polygon.h"

#undef min
#undef max

template <class RT>
class CTriangulation : public RT 
{
public:
    typedef CTriangulation<RT> Rt;
    
    typedef typename Rt::Geom_traits    Kernel;
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Weight     Weight;
    typedef typename Kernel::Point_2    Point;
    typedef typename Kernel::Vector_2   Vector;
    typedef typename Kernel::Ray_2      Ray;
    typedef typename Kernel::Line_2     Line;
    typedef typename Kernel::Segment_2  Segment;
    typedef typename Kernel::Triangle_2 Triangle;
    typedef typename Kernel::Weighted_point_2 Weighted_point;
    
    typedef typename Rt::Vertex                   Vertex;
    typedef typename Rt::Vertex_handle            Vertex_handle;
    typedef typename Rt::Vertex_iterator          Vertex_iterator;
    typedef typename Rt::Vertex_circulator        Vertex_circulator;
    typedef typename Rt::Finite_vertices_iterator Finite_vertices_iterator;
    
    typedef typename Rt::Edge                  Edge;
    typedef typename Rt::Edge_iterator         Edge_iterator;
    typedef typename Rt::Edge_circulator       Edge_circulator;
    typedef typename Rt::Finite_edges_iterator Finite_edges_iterator;
    
    typedef typename Rt::Face                  Face;
    typedef typename Rt::Face_handle           Face_handle;
    typedef typename Rt::Face_iterator         Face_iterator;
    typedef typename Rt::Face_circulator       Face_circulator;
    typedef typename Rt::Finite_faces_iterator Finite_faces_iterator;
    
    typedef CConvexPolygon<Kernel> ConvexPolygon;
    
private:
    ConvexPolygon m_boundary;
    
public:
    CTriangulation() 
    {
    }
    
    void set_boundary(FT dx, FT dy)
    {
        m_boundary.clear();
        m_boundary.init_rectangle(dx, dy);
    }
    
    ////////////
    // ACCESS //
    ////////////
    Vertex_handle get_source(const Edge& edge) const
    {
        return edge.first->vertex( Rt::ccw(edge.second) );
    }    
    Vertex_handle get_target(const Edge& edge) const
    {
        return edge.first->vertex( Rt::cw(edge.second) );        
    }
    Vertex_handle get_opposite(const Edge& edge) const
    {
        return edge.first->vertex( edge.second );
    }    
    Edge get_twin(const Edge& edge) const
    {
        Face_handle f = edge.first;
        Vertex_handle v = get_source(edge);
        Face_handle nf = f->neighbor(edge.second);
        return Edge(nf, Rt::ccw(nf->index(v)));
    }
    Edge get_next(const Edge& edge) const
    {
        Face_handle f = edge.first;
        int index = Rt::ccw(edge.second);
        return Edge(f, index);
    }
    Edge get_prev(const Edge& edge) const
    {
        Face_handle f = edge.first;
        int index = Rt::cw(edge.second);
        return Edge(f, index);
    }
    FT get_length(const Edge& edge) const
    {
        Segment segment = get_segment(edge);
        return std::sqrt(segment.squared_length());
    }
    Segment get_segment(const Edge& edge) const
    {
        const Point& ps = get_source(edge)->get_position();
        const Point& pt = get_target(edge)->get_position();
        return Segment(ps, pt);        
    }
    FT get_area(Face_handle face) const
    {
        Triangle triangle = get_triangle(face);
        return triangle.area();
    }
    Triangle get_triangle(Face_handle face) const
    {
        Vertex_handle v0 = face->vertex(0);
        Vertex_handle v1 = face->vertex(1);
        Vertex_handle v2 = face->vertex(2);
        return Triangle(v0->get_position(), v1->get_position(), v2->get_position());
    }
    Vector get_orthogonal_vector(const Edge& edge) const
    {
        const Point& ps = get_source(edge)->get_position();
        const Point& pt = get_target(edge)->get_position();
        Vector vst = pt - ps;
        return Vector(-vst.y(), vst.x());
    }    
    FT get_average_length() const
    {
        unsigned nb = 0;
        FT avg_length = 0.0;
        for (Finite_edges_iterator
             eit  = RT::finite_edges_begin(); 
             eit != RT::finite_edges_end();
             ++eit)
        {
            Edge edge = *eit;
            FT length = get_length(edge);
            avg_length += length;
            nb++;
        }
        return (avg_length / nb);
    }
    
    //////////
    // AREA //
    //////////
    FT compute_area() const
    {
        FT area = 0.0;
        for (Finite_vertices_iterator
             vit = RT::finite_vertices_begin();
             vit != RT::finite_vertices_end();
             vit++)
        {
            Vertex_handle vi = vit;
            area += vi->compute_area();
        }
        return area;
    }
        
    ///////////////////////
    // INSIDE / BOUNDARY //
    ///////////////////////   
    bool is_inside(Face_handle face) const
    {
        if (RT::is_infinite(face)) return false;
        return true;
    }
    bool is_inside(const Edge& edge) const
    {
        Edge twin = get_twin(edge);
        bool left = is_inside(edge.first);
        bool right = is_inside(twin.first);
        return (left || right);
    }
    bool is_boundary(const Edge& edge) const
    {
        Edge twin = get_twin(edge);
        bool left = is_inside(edge.first);
        bool right = is_inside(twin.first);
        return (left != right);
    }
    bool is_boundary(Vertex_handle vertex) const
    {
        if (vertex->is_hidden()) return false;
        Face_circulator fcirc = RT::incident_faces(vertex);
        Face_circulator fend  = fcirc;
        CGAL_For_all(fcirc, fend)
        {
            Face_handle face = fcirc;
            if (!is_inside(face)) 
                return true;
        }
        return false;
    }
    
    //////////
    // DUAL //
    //////////
    Point get_dual(Face_handle face) const
    {
        return RT::dual(face);
    }
    Segment get_dual(const Edge& edge) const
    {
        /*
         * primal edge = (s,t)
         * dual edge = (right, left)
         * primal x dual > 0
         */
        Edge twin = get_twin(edge);
        Face_handle left_face  = edge.first;
        Face_handle right_face = twin.first;
        bool left_inside  = is_inside( left_face);
        bool right_inside = is_inside(right_face);

        if (left_inside && right_inside)
        {
            Point left_cw  = get_dual( left_face);
            Point right_cw = get_dual(right_face);
            return Segment(right_cw, left_cw);
        }
        
        Vector vec90 = get_orthogonal_vector(edge);
        if (!left_inside && !right_inside)
        {
            Point cw = get_edge_cw(edge);
            Line line(cw, vec90);
            return Segment(line.point(-100), line.point(100));
        }
        
        if (left_inside)
        {
            Point cw = get_dual(left_face);
            Ray ray(cw, -vec90);
            return Segment(ray.point(100), cw);
        }   
        Point cw = get_dual(right_face);
        Ray ray(cw, vec90);
        return Segment(cw, ray.point(100));
    }
    Point get_edge_cw(const Edge& edge) const
    {
        Vertex_handle vi = get_source(edge);
        Vertex_handle vj = get_target(edge);
        const Point& pi = vi->get_position();
        const Point& pj = vj->get_position();
        if (pi == pj) return pi;

        const FT wi = vi->get_weight();
        const FT wj = vj->get_weight();
        const FT lij = get_length(edge);
        const FT dij = 0.5*(lij + (wi - wj)/lij);
        Vector vecij = (pj - pi) / lij;
        return pi + dij*vecij;
    }
    
    /////////////
    // BOUNDED //
    /////////////
    void build_polygon(Vertex_handle vi, 
                       std::vector<Point>& points) const
    {
        points = vi->get_dual().get_points();
    }

    bool pre_build_polygon(Vertex_handle vi,
                           std::vector<Point>& points) const
    {
        std::vector<Segment> segments;
        Edge_circulator ecirc = RT::incident_edges(vi);
        Edge_circulator eend  = ecirc;
        CGAL_For_all(ecirc, eend)
        {
            Edge edge = *ecirc;
            if (!is_inside(edge)) continue;
            
            Edge twin = get_twin(edge);
            Segment segment = build_bounded_dual_edge(twin);
            if (segment.is_degenerate()) continue;
            segments.push_back(segment);
        }
        if (segments.empty()) return false;
        
        Point center;
        Vector size;
        m_boundary.compute_bbox(center, size);
        
        Point first_pt = segments.front().source();        
        Point last_pt = first_pt;
        for (unsigned i = 0; i < segments.size(); ++i)
        {
            Segment segment = segments[i];
            Point ps = segment.source();
            Point pt = segment.target();
            
            if (ps != last_pt)
            {
                fill_gap(center, size, last_pt, ps, points);
                points.push_back(ps);
            }
            
            points.push_back(pt);
            last_pt = pt;
        }
        
        if (first_pt != last_pt)
        {
            fill_gap(center, size, last_pt, first_pt, points);
            points.push_back(first_pt);
        }
        return true;
    }

    void fill_gap(const Point& center, const Vector& size,
                  const Point& a, const Point& b, 
                  std::vector<Point>& points) const
    {
        if (a == b) return;
        
        int aside = find_side(center, size, a);
        if (aside == -1) return;
        
        int bside = find_side(center, size, b);
        if (bside == -1) return;
        
        while (aside != bside)
        {
            aside = (aside + 1) % 4;
            Point q = compute_corner(center, size, aside);
            points.push_back(q);
        }
    }
    
    int find_side(const Point& center, const Vector& size, const Point& a) const
    {
        FT pw = center.x() + size.x();
        FT nw = center.x() - size.x();
        FT ph = center.y() + size.y();
        FT nh = center.y() - size.y();

        // rect
        int side = -1;
        if (std::abs(a.x() - pw) < EPS) side = 1;
        if (std::abs(a.x() - nw) < EPS) side = 3;
        if (std::abs(a.y() - ph) < EPS) side = 2;
        if (std::abs(a.y() - nh) < EPS) side = 0;
        return side;
    }
    
    Point compute_corner(const Point& center, const Vector& size, int side) const
    {
        if (side == 0) return Point(center.x() - size.x(), center.y() - size.y());
        if (side == 1) return Point(center.x() + size.x(), center.y() - size.y());
        if (side == 2) return Point(center.x() + size.x(), center.y() + size.y());
        if (side == 3) return Point(center.x() - size.x(), center.y() + size.y());
        return center;
    }
    
    Segment build_bounded_dual_edge(const Edge& edge) const
    {
        Segment dual_segment = get_dual(edge);
        dual_segment = m_boundary.clamp(dual_segment);        
        return dual_segment;
    }
    
    ////////////
    // LOCATE //
    ////////////
    
    Vertex_handle find_nearest_vertex(const Point& query,
                                      Vertex_handle candidate = Vertex_handle()) const
    {
        typename Kernel::Compare_power_distance_2 cmp_power_distance =
        Rt::geom_traits().compare_power_distance_2_object();
        
        Vertex_handle vertex = candidate;
        if (vertex == Vertex_handle()) vertex = Rt::finite_vertex();

        Vertex_handle vclosest;
        do {
            vclosest = vertex;
            Weighted_point wp = vertex->point();
            
            Vertex_circulator vcirc = Rt::incident_vertices(vertex);
            Vertex_circulator vend  = vcirc;
            CGAL_For_all(vcirc, vend)
            {
                Vertex_handle v = vcirc;
                if (this->is_infinite(v)) continue;
                
                if (cmp_power_distance(query, v->point(), wp) == CGAL::SMALLER ) 
                {
                    vertex = v;
                    break;
                }
            }
        } while (vclosest != vertex);
        return vclosest;  
    }

  /*  std::vector<Vertex_handle> find_neighbors(Vertex_handle vi){
        std::vector<Vertex_handle> neighbors;

        Edge_circulator ecirc = incident_edges(vi);
        Edge_circulator eend  = ecirc;
        if (vi->is_hidden()) return neighbors;
        CGAL_For_all(ecirc, eend)
            {
                Edge edge = *ecirc;
                if (!is_inside(edge)) continue;
                std::cout << "we have a neighbor here" << std::endl;
                Vertex_handle vj = get_source(edge);
                if (vj == vi) vj = get_target(edge);
            }

        return neighbors;
    }
*/

};

#endif
