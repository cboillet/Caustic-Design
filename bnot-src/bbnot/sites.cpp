#include "scene.h"
#include "random.h"

void Scene::generate_random_sites(const unsigned nb)
{
    std::vector<Point> points;
    double dx = m_domain.width();
    double dy = m_domain.height();
    while (points.size() != nb)
    {
        double x = random_double(-dx, dx);
        double y = random_double(-dy, dy);
        Point p(x, y);
        if (m_domain.is_inside(p))
            points.push_back(p);
    }
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);
}

void Scene::generate_regular_grid(const unsigned nx, const unsigned ny)
{
    double dx = m_domain.width();
    double dy = m_domain.height();
    FT stepx = 2.0 * dx / nx;
    FT stepy = 2.0 * dy / ny;
    std::vector<Point> points;
    for (unsigned i = 0; i < nx; ++i)
    {
        FT x = (i + 0.5)*stepx - dx;
        for (unsigned j = 0; j < ny; ++j)
        {
            FT y = (j + 0.5)*stepy - dy;
            Point p(x, y);
            if (m_domain.is_inside(p))
                points.push_back(p);
        }
    }    
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);
}

void Scene::generate_hextille_grid(const unsigned nb)
{
    double dx = m_domain.width();
    double dy = m_domain.height();
    FT step = 2.0 * dx / nb;
    std::vector<Point> points;
    for (unsigned i = 0; i < nb; ++i)
    {
        FT y = (i + 0.5)*step - dy;
        for (unsigned j = 0; j < nb; ++j)
        {
            FT x = (j + 0.25)*step - dx;
            if (i % 2 == 1) x += 0.5*step;
            Point p(x, y);
            if (m_domain.is_inside(p))
                points.push_back(p);
        }
    }
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);    
}
