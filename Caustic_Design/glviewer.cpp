// Qt
#include <QtGui>

// local
#include "scene.h"
#include "glviewer.h"

GlViewer::GlViewer(QWidget *pParent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), pParent)
{
    m_scene = NULL;
    
    m_view_image = false;
    m_view_image_grid = false;
    m_view_domain = true;
    m_view_points = true;
    m_view_vertices = false;
    m_view_edges = true;
    m_view_faces = false;
    m_view_weights = true;
    m_view_dual = false;
    m_view_capacity = false;
    m_view_variance = false;
    m_view_regularity = false;
    m_view_regular_sites = false;
    m_view_pixels = false;
    m_view_barycenter = false;
    m_view_bounded_dual = true;;
    m_view_weight_histogram = false;
    m_view_capacity_histogram = false;

    m_line_thickness = 2.0;
    m_point_size = 2.0;
    m_vertex_size = 3.0;
    
    m_histogram_range = 1.0;
    m_histogram_nbins = 512;
    
    m_scale = 0.5;
    m_center_x = 0.;
    m_center_y = 0.;
    
    setAutoFillBackground(false);
}

void GlViewer::set_scene(Scene* scene)
{
    m_scene = scene;
}

void GlViewer::resizeGL(int width, int height) 
{
    glViewport(0, 0, width, height);
    double aspect_ratio = double(height) / double(width);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -aspect_ratio, aspect_ratio, -1.0, 1.0);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void GlViewer::initializeGL() 
{
    glClearColor(1., 1., 1., 0.);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_SMOOTH);
}

void GlViewer::paintGL() 
{
    glClear(GL_COLOR_BUFFER_BIT);
    if (!m_scene) return;

    // mesh
    glPushMatrix();
    glLoadIdentity();
    glScaled(m_scale, m_scale, m_scale);
    glTranslated(-m_center_x, -m_center_y, 0.0);

    if (m_view_domain)
        m_scene->draw_domain(m_line_thickness, 0.2, 0.2, 0.2);
    
    if (m_view_image)
        m_scene->draw_image();
    
    if (m_view_variance)
        m_scene->draw_variance();
    
    if (m_view_capacity)
        m_scene->draw_capacity();
    
    if (m_view_regularity)
        m_scene->draw_regularity();
    
    if (m_view_regular_sites)
        m_scene->draw_regular_sites();

    if (m_view_pixels)
        m_scene->draw_pixels();
    
    if (m_view_image_grid)
        m_scene->draw_image_grid();
    
    if (m_view_faces)
        m_scene->draw_faces(1.0, 0.7, 0.7);
    
    if (m_view_edges)
        m_scene->draw_primal(m_line_thickness, 0., 0., 0.5);
    
    if (m_view_dual)
        m_scene->draw_dual(m_line_thickness, 0., 0.5, 0.);

    if (m_view_bounded_dual)
        m_scene->draw_bounded_dual(m_line_thickness, 0.5, 0., 0.);

    if (m_view_weights)
        m_scene->draw_weights();
    
    if (m_view_vertices)
        m_scene->draw_vertices(m_vertex_size);
    
    if (m_view_barycenter)
        m_scene->draw_barycenter(m_point_size, 0.5, 0.0, 0.0);
    
    if (m_view_points)
        m_scene->draw_sites(m_point_size, 0., 0., 0.);

    if(m_view_movement)
        m_scene->draw_movement();

    if(m_view_Xrs)
        m_scene->draw_Xrs();

    if(m_view_Xr)
        m_scene->draw_Xr();

    glPopMatrix();
    
    // histograms
    glPushMatrix();
    glLoadIdentity();
    
    if (m_view_capacity_histogram)
        m_scene->draw_capacity_histogram(histogram_nbins(),
                                         -0.8, 0.8, -0.5, 0.3);
    if (m_view_weight_histogram)
        m_scene->draw_weight_histogram(histogram_range(),
                                       histogram_nbins(),
                                       -0.8, 0.8, -0.5, 0.3);
    
    glPopMatrix();
}

void GlViewer::wheelEvent(QWheelEvent *event) 
{
    if (!m_scene) return;
    m_scale += 0.05 * (event->delta() / 120);
    if (m_scale <= 0.0) m_scale = 0.0;
    updateGL();
}

void GlViewer::mousePressEvent(QMouseEvent *event) 
{
    if (!m_scene) return;
    m_mouse_click = event->pos();
    
    if (event->button() == Qt::RightButton)
    {
        setCursor(QCursor(Qt::ClosedHandCursor));
    }
}

void GlViewer::mouseMoveEvent(QMouseEvent *event)
{
    if(!m_scene) return;
    m_mouse_move = event->pos();
    
    if (event->buttons() == Qt::RightButton)
    {
        move_camera(m_mouse_click, m_mouse_move);
    }
    
    m_mouse_click = m_mouse_move;
    updateGL();
}

void GlViewer::mouseReleaseEvent(QMouseEvent *event) 
{
    if (!m_scene) return;
    m_mouse_move = event->pos();
    
    if (event->button() == Qt::RightButton)
    {
        move_camera(m_mouse_click, m_mouse_move);
    }
    
    m_mouse_click = m_mouse_move;
    setCursor(QCursor(Qt::ArrowCursor));
    updateGL();
}

void GlViewer::move_camera(const QPoint& p0, const QPoint& p1)
{
    m_center_x -= double(p1.x() - p0.x()) / double(width());
    m_center_y += double(p1.y() - p0.y()) / double(height());
}


