#include "scene.h"
#include "random.h"

void Scene::generate_random_sites(const unsigned nb)
{
    std::vector<Point> points;
    for (unsigned i = 0; i < nb; ++i)
    {
        double x = random_double(0.0, 1.0);
        double y = random_double(0.0, 1.0);
        points.push_back(Point(x, y));
    }

    std::vector<FT> noise(points.size(), 0.0);
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights, noise);
    std::cout << "Insert " << points.size() << std::endl;
}

void Scene::generate_regular_grid(const unsigned nx, const unsigned ny)
{
    FT stepx = 1.0 / nx;
    FT stepy = 1.0 / ny;
    std::vector<Point> points;
    for (unsigned i = 0; i < nx; ++i)
    {
        FT x = (i + 0.5)*stepx;
        for (unsigned j = 0; j < ny; ++j)
        {
            FT y = (j + 0.5)*stepy;
            points.push_back(Point(x, y));
        }
    }
    
    std::vector<FT> noise(points.size(), 0.0);
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights, noise);
    std::cout << "Insert " << points.size() << std::endl;
}

void Scene::generate_hextille_grid(const unsigned nb)
{
    FT step = 1.0 / nb;
    std::vector<Point> points;
    for (unsigned i = 0; i < nb; ++i)
    {
        FT y = (i + 0.5)*step;
        for (unsigned j = 0; j < nb; ++j)
        {
            FT x = (j + 0.5)*step;
            if (i % 2 == 1) x += 0.5*step;
            points.push_back(Point(x, y));
        }
    }
    
    std::vector<FT> noise(points.size(), 0.0);
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights, noise);
    std::cout << "Insert " << points.size() << std::endl;
}
