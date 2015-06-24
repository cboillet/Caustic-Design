// 
#include <cmath>

// Qt
#include <QtOpenGL>

// local
#include "ramp.h"
#include "scene.h"
#include "random.h"

////////////
// BASICS //
////////////

void Scene::draw_point(const Point& point) const
{
    glBegin(GL_POINTS);
    glVertex2d(point.x(), point.y());
    glEnd();    
}

void Scene::draw_segment(const Point& s, const Point& t) const
{
	glBegin(GL_LINES);
    glVertex2d(s.x(), s.y());
    glVertex2d(t.x(), t.y());
	glEnd();
}

void Scene::draw_triangle(const Point& a, const Point& b, const Point& c) const
{
	glBegin(GL_TRIANGLES);
    glVertex2d(a.x(), a.y());
    glVertex2d(b.x(), b.y());
    glVertex2d(c.x(), c.y());
	glEnd();
}

void Scene::draw_polygon(const std::vector<Point>& polygon) const
{
    if (polygon.size() < 3) return;
    glBegin(GL_POLYGON);
    for (unsigned i = 0; i < polygon.size(); ++i)
    {
        const Point& p = polygon[i];
        glVertex2d(p.x(), p.y());
    }
    glEnd();
}

void Scene::draw_circle(const Point& center, 
                        const FT scale,
                        const std::vector<Point>& pts) const
{
    glPushMatrix();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glTranslated(center.x(), center.y(), 0.0);
    glScaled(scale, scale, scale);
    draw_polygon(pts);
    glPopMatrix();
}

/////////////
// OBJECTS //
/////////////

void Scene::draw_image() const
{
    if (!m_domain.is_valid()) return;
    m_domain.draw_image();
}

void Scene::draw_image_grid() const
{
    if (!m_domain.is_valid()) return;
    glLineWidth(1.0);
    m_domain.draw_grid();
}

void Scene::draw_domain(const float line_width,
                        const float red,
                        const float green,
                        const float blue) const
{
    if (!m_domain.is_valid()) return;
    glLineWidth(line_width);
    glColor3d(red, green, blue);
    m_domain.draw_boundary();
}

void Scene::draw_sites(const float point_size,
                       const float red,
                       const float green,
                       const float blue) const
{
    glPointSize(point_size);
    glColor3d(red, green, blue);
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        const Point& pi = vi->get_position();
        draw_point(pi);
    }
}

void Scene::draw_barycenter(const float point_size,
                            const float red,
                            const float green,
                            const float blue) const
{
    glPointSize(point_size);
    glColor3d(red, green, blue);
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        const Point& pi = vi->compute_centroid();
        draw_point(pi);
    }
}

void Scene::draw_movement() const
{
    glColor3d(1.0, 1.0, 0.0);
    for (unsigned i = 0; i < m_vertices.size(); i++)
    {
        Vertex_handle vi = m_vertices[i];
        if(vi->is_hidden()) continue;
        const Point& centroid = vi->compute_centroid();
        const Point& site = vi->get_position();
        draw_segment(centroid, site);
    }
}

void Scene::draw_vertices(const float point_size) const
{
    glPointSize(point_size);

    glColor3d(1., 0., 0.);
    for (Finite_vertices_iterator 
         vit = m_rt.finite_vertices_begin(); 
         vit != m_rt.finite_vertices_end(); 
         ++vit)
    {
        Vertex_handle vi = vit;
        const Point& pi = vi->get_position();
        draw_point(pi);
    }
    
    glColor3d(1., 1., 0.);
    for (Hidden_vertices_iterator 
         vit = m_rt.hidden_vertices_begin(); 
         vit != m_rt.hidden_vertices_end(); 
         ++vit)
    {
        Vertex_handle vi = vit;
        const Point& pi = vi->get_position();
        draw_point(pi);
    }    
}

void Scene::draw_faces(const float red, 
                       const float green,
                       const float blue) const
{
    for (Finite_faces_iterator
         fit = m_rt.finite_faces_begin();
         fit != m_rt.finite_faces_end();
         ++fit)
    {
        Face_handle face = fit;
        if (m_rt.is_inside(face))
            glColor3f(red, green, blue);
        else
            glColor3f(1.0-red, 1.0-green, 1.0-blue);
        
        const Point& pa = face->vertex(0)->get_position();
        const Point& pb = face->vertex(1)->get_position();
        const Point& pc = face->vertex(2)->get_position();
        draw_triangle(pa, pb, pc);
    }
}

void Scene::draw_primal(const float line_width,
                        const float red,
                        const float green,
                        const float blue) const
{
    glLineWidth(line_width);
    for (Finite_edges_iterator
         eit = m_rt.finite_edges_begin(); 
         eit != m_rt.finite_edges_end();
         ++eit)
    {
        Edge edge = *eit;
        if (m_rt.is_inside(edge))
            glColor3f(red, green, blue);
        else
            glColor3f(1.0-red, 1.0-green, 1.0-blue);
            
        const Point& pa = m_rt.get_source(edge)->get_position();
        const Point& pb = m_rt.get_target(edge)->get_position();
        draw_segment(pa, pb);
    }
}

void Scene::draw_bounded_dual(const float line_width,
                              const float red,
                              const float green,
                              const float blue) const
{
    glLineWidth(line_width);
    glColor3d(red, green, blue);
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        draw_cell(vi, false);
    }
}

void Scene::draw_dual(const float line_width,
                      const float red,
                      const float green,
                      const float blue) const
{
    glLineWidth(line_width);
    for (Finite_edges_iterator
         eit = m_rt.finite_edges_begin();
         eit != m_rt.finite_edges_end();
         ++eit)
    {
        Edge edge = *eit;
        
        glColor3d(1.0-red, 1.0-green, 1.0-blue);
        Segment segment = m_rt.get_dual(edge);
        draw_segment(segment[0], segment[1]);
        
        glColor3d(red, green, blue);
        segment = m_rt.build_bounded_dual_edge(edge);
        draw_segment(segment[0], segment[1]);
    }
}

void Scene::draw_weights() const
{
    std::vector<Point> points;
    generate_regular_polygon(100, 1.0, Point(0.0, 0.0), std::back_inserter(points));

    FT min_w = 0.0;
    FT max_w = 0.0;
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        FT w = vi->get_weight();
        min_w = std::min(min_w, w);
        max_w = std::max(max_w, w);
    }
    FT range_w = max_w - min_w;
    if (range_w == 0.0) range_w = 1.0;
    
    std::vector<FT> value;
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) 
        {
            value.push_back(0.0);
            continue;
        }
        FT w = vi->get_weight();
        FT t = (w - min_w) / range_w;
        value.push_back(t);
    }
    
    FT avg_len = 0.25*m_rt.get_average_length();
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        FT t = value[i];
        
        QColor c;
        c.setHsv(360*(1.0-t) + 120*t, 255, 225);
        glColor3d(c.redF(), c.greenF(), c.blueF());
        
        const Point& p = vi->get_position();
        draw_circle(p, avg_len, points);
    }
}

void Scene::draw_pixels() const
{
    for (Finite_vertices_iterator 
         vit = m_rt.finite_vertices_begin();
         vit != m_rt.finite_vertices_end();
         ++vit)
    {
        Vertex_handle vertex = vit;
        double r = m_r[vertex->get_index()];
        double g = m_g[vertex->get_index()];
        double b = m_b[vertex->get_index()];
        glColor3d(r, g, b);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        for (unsigned i = 0; i < vertex->nb_pixels(); ++i)
        {
            const Pixel& pixel = vertex->get_pixel(i);
            const ConvexPolygon& shape = pixel.get_shape();
            draw_polygon(shape.get_points());
        }
    }
}

void Scene::draw_capacity() const
{
    FT min_ratio = 1.0;
    FT max_ratio = 0.0;
    std::vector<FT> values;
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        FT Ci = m_capacities[i];
        FT Vi = vi->compute_area();
        FT ratio = Vi / Ci;
        values.push_back(ratio);
        max_ratio = std::max(max_ratio, ratio);
        min_ratio = std::min(min_ratio, ratio);
    }
    
    Ramp ramp;
    unsigned j = 0;
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        FT ratio = values[j++];
        ramp.gl_color(ratio, max_ratio);
        draw_cell(vi, true);
    }
}

void Scene::draw_regularity() const
{
    std::vector<FT> variance;
    compute_variance_vector(variance);
    
    std::vector<FT> regularity;
    compute_regularity_vector(variance, regularity);
    
    FT max_value = 0.0;
    for (unsigned i = 0; i < regularity.size(); ++i)
        max_value = std::max(max_value, regularity[i]);
    if (max_value < 1e-8) max_value = 1.0;
    
    Ramp ramp;
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        
        FT value = regularity[vi->get_index()];
        ramp.gl_color(value, max_value);
        draw_cell(vi, true);
    }
}

void Scene::draw_regular_sites() const
{
    std::vector<FT> variance;
    compute_variance_vector(variance);
    
    std::vector<FT> regularity;
    compute_regularity_vector(variance, regularity);
    
    FT threshold = compute_regularity_threshold();
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        
        FT value = regularity[vi->get_index()];
        if (value > threshold) continue;
        glColor3d(0.7, 0.2, 0.2);
        draw_cell(vi, true);
    }
}

void Scene::draw_variance() const
{
    FT max_variance = 0.0;
    std::vector<FT> values;
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        
        FT variance = vi->compute_variance();
        values.push_back(variance);
        max_variance = std::max(max_variance, variance);
    }
    
    Ramp ramp;
    unsigned j = 0;
    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        
        FT variance = values[j++];
        ramp.gl_color(variance, max_variance);
        draw_cell(vi, true);
    }
}

void Scene::draw_cell(Vertex_handle vertex, bool filled) const
{
    std::vector<Point> polygon;
    m_rt.build_polygon(vertex, polygon);
 
    if (filled) 
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    draw_polygon(polygon);
}

void Scene::draw_capacity_histogram(const unsigned nbins, 
                                    const double xmin,
                                    const double xmax,
                                    const double ymin,
                                    const double ymax) const
{
    std::vector<unsigned> histogram(nbins, 0);
    compute_capacity_histogram(histogram);
    draw_histogram(histogram, xmin, xmax, ymin, ymax);
}
    
void Scene::draw_weight_histogram(const double range,
                                  const unsigned nbins, 
                                  const double xmin,
                                  const double xmax,
                                  const double ymin,
                                  const double ymax) const
{
    std::vector<unsigned> histogram(nbins, 0);
    compute_weight_histogram(range, histogram);
    draw_histogram(histogram, xmin, xmax, ymin, ymax);    
}

void Scene::draw_histogram(const std::vector<unsigned>& histogram, 
                           const double xmin,
                           const double xmax,
                           const double ymin,
                           const double ymax) const
{
    unsigned nbins = histogram.size();
    if (nbins < 2) return;

    // horizontal line
    glColor3d(0., 0., 0.);
    glBegin(GL_LINES);
    glVertex2d(xmin, ymin);
    glVertex2d(xmax, ymin);
    glEnd();
    
    // mean line
    double xstep = (xmax - xmin) / (nbins-1);
    unsigned offset = unsigned(ceil(0.5*nbins));
    glColor3d(1., 0., 1.);
    glBegin(GL_LINES);
    glVertex2d(xmin + offset*xstep, ymin);
    glVertex2d(xmin + offset*xstep, ymax);
    glEnd();

    unsigned max_value = 0;
    for (unsigned i = 0; i < histogram.size(); ++i)
        max_value = std::max(max_value, histogram[i]);
    if (max_value == 0) return;

    glColor3d(0.3, 0.3, 0.3);
    for (unsigned i = 0; i < nbins; ++i)
    {
        double height = double(histogram[i]) / double(max_value);
        double x0 = xmin + i*xstep;
        double y0 = ymin;
        double y1 = ymin + (ymax - ymin)*height;
     
        glBegin(GL_LINES);
        glVertex2d(x0, y0);
        glVertex2d(x0, y1);
        glEnd();
    }
}

void Scene::draw_Xrs(){
    int i;
    glPointSize(1);
    glColor3d(1.0,1.0,0.0);
    //Point p;
    for (i=0; i<lightpts.size(); ++i) {
        draw_point(lightpts[i]);
    }
}

void Scene::draw_Xr(){
    int i;
    glPointSize(2);
    glColor3d(1.0,1.0,0.0);
    //Point p;
    for (i=0; i<lightpt.size(); ++i) {
        std::cout<<"lightpt x"<<lightpt[i].x()<<std::endl;
        std::cout<<"lightpt y"<<lightpt[i].x()<<std::endl;
        draw_point(lightpt[i]);
    }
}

