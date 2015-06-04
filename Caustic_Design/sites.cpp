#include "scene.h"
#include "random.h"

Scene::Scene(const Scene& sc){
    srand(0);
    m_tau = 1.0;
    m_timer_on = false;
    m_fixed_connectivity = false;
    m_rt = sc.m_rt;
    m_domain = sc.m_domain;
    m_capacities = sc.m_capacities;
    m_tau = sc.m_tau;
    m_ratio = sc.m_ratio;
    m_r = sc.m_r; m_g=sc.m_g; m_b=sc.m_b;
    m_vertices= sc.m_vertices;
}

void Scene::generate_random_sites(const unsigned nb)
{
    if (!m_domain.is_valid()) return;
    std::vector<Point> points;
    double dx = m_domain.get_dx();
    double dy = m_domain.get_dy();
    while (points.size() != nb)
    {
        double x = random_double(-dx, dx);
        double y = random_double(-dy, dy);
        points.push_back( Point(x, y) );
    }
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);
    init_colors(points.size());
}

void Scene::generate_random_sites_based_on_image(const unsigned nb)
{
    if (!m_domain.is_valid()) return;
    std::vector<Point> points;
    double dx = m_domain.get_dx();
    double dy = m_domain.get_dy();
    while (points.size() != nb)
    {
        double x = random_double(-dx, dx);
        double y = random_double(-dy, dy);
        
        Point point(x, y);
        double prob  = random_double(0.0, 1.0);
        double value = m_domain.get_value(point, true) - PIXEL_EPS;
        if (prob < value) points.push_back(point);
    }
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);
    init_colors(points.size());
}

void Scene::generate_regular_grid(const unsigned nx, const unsigned ny)
{
    if (!m_domain.is_valid()) return;
    FT stepx = 2.0 * m_domain.get_dx() / nx;
    FT stepy = 2.0 * m_domain.get_dy() / ny;
    std::vector<Point> points;
    for (unsigned i = 0; i < nx; ++i)
    {
        FT x = (i + 0.5)*stepx - m_domain.get_dx();
        x += EPS;
        for (unsigned j = 0; j < ny; ++j)
        {
            FT y = (j + 0.5)*stepy - m_domain.get_dy();
            y += EPS;
            points.push_back(Point(x, y));
        }
    }
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);
    init_colors(points.size());
}

void Scene::init_colors(const unsigned nb)
{
    m_r.clear();
    m_g.clear();
    m_b.clear();
    for (unsigned i = 0; i < nb; ++i)
    {
        m_r.push_back(random_double(0.0, 1.0));
        m_g.push_back(random_double(0.0, 1.0));
        m_b.push_back(random_double(0.0, 1.0));
    }
}
