#ifndef GLWIDGET_H
#define GLWIDGET_H

// Qt
#include <QGLWidget>
#include <QPaintEvent>

class Scene;

class GlViewer : public QGLWidget 
{
    Q_OBJECT
    
private:
    Scene* m_scene;
    
    // toggles
    bool m_view_image;
    bool m_view_image_grid;
    bool m_view_domain;
    bool m_view_points;
    bool m_view_vertices;
    bool m_view_edges;
    bool m_view_faces;
    bool m_view_weights;
    bool m_view_dual;
    bool m_view_pixels;
    bool m_view_capacity;
    bool m_view_variance;
    bool m_view_regularity;
    bool m_view_regular_sites;
    bool m_view_barycenter;
    bool m_view_bounded_dual;
    bool m_view_weight_histogram;
    bool m_view_capacity_histogram;    
    
    // rendering options
    double m_line_thickness;
    double m_point_size;
    double m_vertex_size;
    
    // histogram
    double m_histogram_range;
    unsigned m_histogram_nbins;
    
    // camera
    double m_scale;
    double m_center_x;
    double m_center_y;
    
    // mouse
    QPoint m_mouse_click;
    QPoint m_mouse_move;
    QPoint m_mouse_pick;
    
public:
    GlViewer(QWidget *parent);
    
    void set_scene(Scene* scene);
    
    void set_camera(const double x, const double y, const double s) 
    {
        m_center_x = x;
        m_center_y = y;
        m_scale = s;
    }
    
    // options
    double& line_thickness() { return m_line_thickness; }
    const double& line_thickness() const { return m_line_thickness; }
    
    double& point_size() { return m_point_size; }
    const double& point_size() const { return m_point_size; }
    
    double& vertex_size() { return m_vertex_size; }
    const double& vertex_size() const { return m_vertex_size; }
    
    double& histogram_range() { return m_histogram_range; }
    const double histogram_range() const { return m_histogram_range; }

    unsigned& histogram_nbins() { return m_histogram_nbins; }
    const unsigned histogram_nbins() const { return m_histogram_nbins; }

    // toggles
    void toggle_view_image() { m_view_image = !m_view_image; }
    
    void toggle_view_image_grid() { m_view_image_grid = !m_view_image_grid; }

    void toggle_view_domain() { m_view_domain = !m_view_domain; } 

    void toggle_view_points() { m_view_points = !m_view_points; } 

    void toggle_view_vertices() { m_view_vertices = !m_view_vertices; }
    
    void toggle_view_edges() { m_view_edges = !m_view_edges; }
    
    void toggle_view_faces() { m_view_faces = !m_view_faces; }
    
    void toggle_view_weights() { m_view_weights = !m_view_weights; }
    
    void toggle_view_dual() { m_view_dual = !m_view_dual; }

    void toggle_view_pixels() { m_view_pixels = !m_view_pixels; }
    
    void toggle_view_capacity() { m_view_capacity = !m_view_capacity; }
    
    void toggle_view_variance() { m_view_variance = !m_view_variance; }

    void toggle_view_regularity() { m_view_regularity = !m_view_regularity; }
    
    void toggle_view_regular_sites() { m_view_regular_sites = !m_view_regular_sites; }

    void toggle_view_barycenter() { m_view_barycenter = !m_view_barycenter; }

    void toggle_view_bounded_dual() { m_view_bounded_dual = !m_view_bounded_dual; }
    
    void toggle_view_weight_histogram() { m_view_weight_histogram = !m_view_weight_histogram; }
    
    void toggle_view_capacity_histogram() { m_view_capacity_histogram = !m_view_capacity_histogram; }

protected:
    // GL
    void paintGL();
    void initializeGL();
    void resizeGL(int width, int height);
    
    // mouse
    void wheelEvent(QWheelEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void move_camera(const QPoint& p0, const QPoint& p1);
};

#endif
