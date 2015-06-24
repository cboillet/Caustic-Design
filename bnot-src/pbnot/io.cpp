#include <QColor>
#include "scene.h"

unsigned Scene::load(const QString& filename)
{
    std::vector<Point> points;
    std::vector<FT> weights;
    std::vector<FT> noise;
    
    if (filename.contains(".dat", Qt::CaseInsensitive))
        load_dat_file(filename, points);
    
    if (filename.contains(".txt", Qt::CaseInsensitive))
        load_txt_file(filename, points);
    
    if (points.empty()) 
    {
        std::cout << red << "try (.dat) or (.txt) file format" << white << std::endl;
        return 0;
    }
    
    weights.resize(points.size(), 0.0);
    noise.resize(points.size(), 0.0);
    
    clear();
    construct_triangulation(points, weights, noise);
    return m_vertices.size();
}

void Scene::load_dat_file(const QString& filename, std::vector<Point>& points) const
{
    std::ifstream ifs(qPrintable(filename));
    Point point;
    while (ifs >> point) points.push_back(point);
    ifs.close();
}

void Scene::load_txt_file(const QString& filename, std::vector<Point>& points) const
{
    std::ifstream ifs(qPrintable(filename));
    
    // kill the first line (contains the number of points)
    int numPoints;
    ifs >> numPoints;
    
    Point point;
    while (ifs >> point) points.push_back(point);
    ifs.close();
}

void Scene::save(const QString& filename) const
{
    std::vector<Point> points;
    collect_visible_points(points);
    
    if (filename.contains(".dat", Qt::CaseInsensitive))
    {
        save_dat(filename, points);
        return;
    }
    
    if (filename.contains(".txt", Qt::CaseInsensitive))
    {
        save_txt(filename, points);
        return;
    }
    
    if (filename.contains(".eps", Qt::CaseInsensitive))
    {
        save_eps(filename);
        return;
    }
    
    std::cout << red << "try (.dat, .txt, .eps) file format" << white << std::endl;
}

void Scene::save_dat(const QString& filename, const std::vector<Point>& points) const
{
    std::ofstream ofs(qPrintable(filename));

    for (unsigned i = 0; i < points.size(); ++i)
    {
        ofs << points[i] << std::endl;
    }
    ofs.close();
}

void Scene::save_txt(const QString& filename, const std::vector<Point>& points) const
{
    std::ofstream ofs(qPrintable(filename));
    ofs << points.size() << std::endl;
    for (unsigned i = 0; i < points.size(); ++i)
    {
        ofs << points[i] << std::endl;
    }
    ofs.close();
}

void Scene::save_eps(const QString& filename) const
{
    double scale  = 512.;
    double radius = 2.75;
    radius /= scale;

    double min = 0;
    double max = scale;

    std::ofstream ofs(qPrintable(filename));
    
    ofs << "%!PS-Adobe-3.1 EPSF-3.0\n";
    ofs << "%%HiResBoundingBox: " 
    << min << " " << min << " " << max << " " << max << std::endl;
    ofs << "%%BoundingBox: " 
    << min << " " << min << " " << max << " " << max << std::endl;
    ofs << "%%CropBox: " 
    << min << " " << min << " " << max << " " << max << "\n";
    
    ofs << "/radius { " << radius << " } def\n";
    ofs << "/p { radius 0 360 arc closepath fill stroke } def\n";
    ofs << "gsave " << scale << " " << scale << " scale\n";
    ofs << "0 0 0 setrgbcolor" << std::endl;

    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        const Point& pi = vi->get_position();
        ofs << pi.x() << " " << pi.y() << " p" << std::endl;
    }

    ofs << "grestore" << std::endl;
    ofs.close();
}
