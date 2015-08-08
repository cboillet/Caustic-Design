#include <QColor>

#include "scene.h"

#include <string>
#include <fstream>
#include <iostream>
#include "tinyxml.h"
#include "console_color.h"

void Scene::load_image(const QString& filename, const int width)
{
    bool ok = m_domain.load(filename, width);
    if (!ok) return;
    
    m_rt.set_boundary(m_domain.get_dx(),
                      m_domain.get_dy());
    std::cout << "Dx vs Dy: " << m_domain.get_dx() << " ; " << m_domain.get_dy() << std::endl;
}

void Scene::load_points(const QString& filename)
{
    if (!m_domain.is_valid()) return;
    
    std::vector<FT> weights;
    std::vector<Point> points;
    if (filename.contains(".dat", Qt::CaseInsensitive))
    {
        load_dat(filename, points);
        weights.resize(points.size(), 0.0);
    }
    
    if (!points.empty())
    {
        clear();
        init_colors(points.size());
        construct_triangulation(points, weights);
        return;
    }
    
    std::cout << red << "try (.dat) file format" << white << std::endl;
    return;
}

void Scene::load_dat(const QString& filename, std::vector<Point>& points) const
{
    std::ifstream ifs(qPrintable(filename));
    Point point;
    while (ifs >> point) points.push_back(point);
    ifs.close();
}

bool Scene::load_interpolation_points(const QString& filename){
    std::ifstream ifs(qPrintable(filename));
    Point point;
    FT x,y;
    FT maxX,minX,maxY,minY;
    std::vector<Point> points;
    while (ifs >> point) {

        points.push_back(point);
    }

    for(int i=0; i<points.size(); i++){
        if (points[i].x()<minX) minX=i;
        if (points[i].y()<minY) minY=i;
        if (points[i].x()>maxX) maxX=i;
        if (points[i].y()>maxY) maxY=i;

    }
    // FT scalingXfactor = (1 / abs(points[maxX].x())) * m_domain.get_dx();
    // FT scalingYfactor = (1 / abs(points[maxY].y())) * m_domain.get_dy();
    // std::cout<<"scalginXfactor"<<scalingXfactor<<std::endl;
     //std::cout<<"scalginYfactor"<<scalingYfactor<<std::endl;


    for(int i=0; i<points.size(); i++){
        x=points[i].x()*SCALING_X;
        y=points[i].y()*SCALING_Y;
        lightpts.push_back(Point(x,y));
//        if(x<=m_domain.get_dx() && y<=m_domain.get_dy()) lightpts.push_back(Point(x,y));
//        else {
//            std::cout<<"x uninserted"<<x<<std::endl;
//            std::cout<<"y uninserted"<<x<<std::endl;
//            return false;
//        }


    }

    ifs.close();
    return true;
}

std::vector<FT> Scene::load_weights(const QString& filename) const
{
    std::vector<FT> weights = std::vector<FT>();

    std::ifstream ifs(qPrintable(filename));
    FT weight;
    while(ifs >> weight) weights.push_back(weight);
    ifs.close();

    return weights;
    //update_weights(weights);
}

void Scene::load_singularities(const QString& filename, std::vector<PointSingularity>& pSing, std::vector<CurveSingularity>& cSing) const
{

    if (!m_domain.is_valid())
    {
        std::cerr << "Domain not valid" << std::endl;
        return;
    }

    FT dx = m_domain.get_dx();
    FT dy = m_domain.get_dy();
    FT width = m_domain.get_width();
    FT height = m_domain.get_height();

    std::cout << "loading singularities from file: " << filename.toStdString() << std::endl;

    pSing.clear();
    cSing.clear();

    // convert qstring to char*
    QByteArray ba = filename.toLatin1();
    const char *file = ba.data();

    // load document
    TiXmlDocument doc(file);
    if(!doc.LoadFile())
    {
        std::cerr << "error loading file" << std::endl;
        return;
    }

    TiXmlHandle handle(&doc);
    TiXmlElement* elem;
    TiXmlHandle root(0);

    // doc name (should be svg)
    elem = handle.FirstChildElement().Element();
    if(!elem)
    {
        std::cerr << "first element not found" << std::endl;
        return;
    }
    root = TiXmlHandle(elem);

    // read all circles
    elem = root.FirstChild("circle").Element();
    TiXmlAttribute* attr;

    do{
        FT x,y;
        for(attr = elem->FirstAttribute(); attr != NULL; attr = attr->Next())
        {
            std::string name = std::string(attr->Name());
            std::string value = std::string(attr->Value());

            if( name == "cx" )
            {
                x = (FT(std::stod(value, NULL)) / width) - dx;
            }else if( name == "cy")
            {
                y = (FT(std::stod(value, NULL)) / height) - dy;
            }
        }

        PointSingularity ps = PointSingularity(Point(x,y), 0.0);
        pSing.push_back(ps);
    }while((elem = elem->NextSiblingElement()) != NULL);

    // todo read paths (into curve-singularities)
}

void Scene::save_points(const QString& filename) const
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

void Scene::save_weights(const QString& filename) const
{
    std::vector<Point> foo = std::vector<Point>();
    std::vector<FT> weights = std::vector<FT>();
    collect_sites(foo, weights);

    std::ofstream ofs(qPrintable(filename));
    ofs.precision(20);

    for(unsigned i = 0; i < weights.size(); i++)
    {
        ofs << weights[i] << std::endl;
    }
    ofs.close();
}

void Scene::save_dat(const QString& filename, const std::vector<Point>& points) const
{
    std::ofstream ofs(qPrintable(filename));
    ofs.precision(20);

    for (unsigned i = 0; i < points.size(); ++i)
    {
        ofs << points[i] << std::endl;
    }
    ofs.close();
}

void Scene::save_txt(const QString& filename, const std::vector<Point>& points) const
{
    std::ofstream ofs(qPrintable(filename));
    ofs.precision(20);

    ofs << points.size() << std::endl;
    for (unsigned i = 0; i < points.size(); ++i)
    {
        ofs << points[i] << std::endl;
    }
    ofs.close();
}

void Scene::save_eps(const QString& filename) const
{
    double dx = m_domain.get_dx();
    double dy = m_domain.get_dy();

    double scale = 512.0;
    double radius = 0.002;

    double wx = dx * scale;
    double wy = dy * scale;

    double min_x = 0;
    double max_x = 2.0 * wx;
    double min_y = 0;
    double max_y = 2.0 * wy;

    std::ofstream ofs(qPrintable(filename));
    ofs.precision(20);

    ofs << "%!PS-Adobe-3.1 EPSF-3.0\n";
    ofs << "%%HiResBoundingBox: " 
    << min_x << " " << min_y << " " << max_x << " " << max_y << std::endl;
    ofs << "%%BoundingBox: " 
    << min_x << " " << min_y << " " << max_x << " " << max_y << std::endl;
    ofs << "%%CropBox: " 
    << min_x << " " << min_y << " " << max_x << " " << max_y << "\n";
    
    ofs << "/radius { " << radius << " } def\n";
    ofs << "/p { radius 0 360 arc closepath fill stroke } def\n";
    ofs << "gsave " << scale << " " << scale << " scale\n";
    ofs << "0 0 0 setrgbcolor" << std::endl;

    for (unsigned i = 0; i < m_vertices.size(); ++i) 
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        const Point& pi = vi->get_position();
        ofs << pi.x() + dx << " " << pi.y() + dy << " p" << std::endl;
    }
    ofs << "grestore" << std::endl;
    ofs.close();
}


void Scene::save_interpolation_dat(const QString &filename) const{
    std::ofstream ofs(qPrintable(filename));
    FT x;
    FT y;

    for(int i=0; i<lightpt.size(); i++){
        x=lightpt[i].x()*SCALING_X;
        y=lightpt[i].y()*SCALING_Y;
        ofs << Point(x,y) << std::endl;
    }

    ofs.close();
}


