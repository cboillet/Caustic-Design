#include "scene.h"
#include "util.h"
#include "timer.h"

bool Scene::is_valid() const
{
    if (!m_domain.is_valid()) return false;
    if (m_vertices.empty()) return false;
    return true;
}

unsigned Scene::count_visible_sites() const
{
    unsigned nb = 0;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        nb++;
    }
    return nb;   
}

void Scene::collect_visible_points(std::vector<Point>& points) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        points.push_back(vi->get_position());
    }
}

void Scene::collect_visible_weights(std::vector<FT>& weights) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        weights.push_back(vi->get_weight());
    }
}

void Scene::collect_sites(std::vector<Point>& points,
                          std::vector<FT>& weights) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        Point pi = vi->get_position();
        points.push_back(pi);
        
        FT wi = 0.0;
        wi = vi->get_weight();
        weights.push_back(wi);
    }
}

void Scene::collect_singularities(std::vector<PointSingularity> &pointSingularities,
                                  std::vector<CurveSingularity> &curveSingularities) const
{
    for (unsigned i = 0; i < m_point_singularities.size(); i++)
    {
        pointSingularities.push_back(m_point_singularities[i]);
    }

    for (unsigned i = 0; i < m_curve_singularities.size(); i++)
    {
        curveSingularities.push_back(m_curve_singularities[i]);
    }
}

void Scene::clear_triangulation()
{
    m_ratio.clear();
    m_vertices.clear();
    m_capacities.clear();
    m_rt.clear();
}

bool Scene::update_triangulation(bool skip)
{
    std::vector<FT> weights;
    std::vector<Point> points;
    collect_sites(points, weights);
    return construct_triangulation(points, weights, skip);
}

bool Scene::construct_triangulation(const std::vector<Point>& points,
                                    const std::vector<FT>& weights,
                                    bool skip)
{
    if (m_timer_on)
    {
        Timer::start_timer(m_timer, COLOR_BLUE, "Triangulation");
        std::cout << std::endl;
    }

    clear_triangulation();
    bool ok = populate_vertices(points, weights);

    if(!ok){
        std::cerr << "Warning, some vertices are hidden" << std::endl;
    }

    if (ok || !skip)
    {
        pre_build_dual_cells();
        assign_pixels();
        assign_singularites();
        pre_compute_area();
        //compute_capacities(m_capacities);
    }
    
    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
    return (ok || !skip);
}

void Scene::update_triangulation_values()
{
    pre_build_dual_cells();
    assign_pixels();
    assign_singularites();
    pre_compute_area();
}

bool Scene::populate_vertices(const std::vector<Point>& points,
                              const std::vector<FT>& weights)
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_YELLOW, "Populate");
    
    unsigned nb = 0;
    unsigned nsites = points.size();
    FT weight_sum = 0.0;
    for (unsigned i = 0; i < nsites; ++i)
    {
        weight_sum += weights[i];
        Vertex_handle vertex = insert_vertex(points[i], weights[i], nb);
        if (vertex == Vertex_handle())
        {
            std::cerr << "did not insert vertex with index " << nb << std::endl;
            continue;
        }
        m_vertices.push_back(vertex);
        nb++;
    }

    std::cout << "weight sum is " << weight_sum << std::endl;
    
    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_YELLOW);

    bool none_hidden = true;
    if (count_visible_sites() != m_vertices.size())
        none_hidden = false;
    
    return none_hidden;
}

Vertex_handle Scene::insert_vertex(const Point& point,
                                   const FT weight,
                                   const unsigned index)
{
    Weighted_point wp(point, weight);
    Vertex_handle vertex = m_rt.insert(wp);

    if (vertex->get_index() != -1)
    {
        std::cerr << "Did not insert point @ (" << point.x() << ", " << point.y() << ") with weight = " << weight << std::endl;
        return Vertex_handle();
    }else{
        //std::cout << "inserted point @ (" << point.x() << ", " << point.y() << ") with weight = " << weight << std::endl;
    }

    vertex->set_index(index);
    return vertex;
}

void Scene::delete_vertex(Vertex_handle vd)
{
    int i= -1;
    i=findIndexVerticeBySite(vd);

    if (i == -1) std::cout << "Not inside the m_vertices, can't be deleted" << std::endl;
    return;
    //m_rt.delete_vertex(vd);
    //int i=findIndexVertice(vd);
    //m_vertices.erase(m_vertices.begin()+i);

}

FT Scene::compute_mean_capacity() const
{
    FT domain_area = compute_value_integral();
    unsigned nb = count_visible_sites();
    return (domain_area / nb);
}

void Scene::compute_capacities(std::vector<FT>& capacities) const
{
    FT C = compute_mean_capacity();
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        FT Ci = 0.0;
        Vertex_handle vi = m_vertices[i];
        if (!vi->is_hidden()) Ci = C;
        capacities.push_back(Ci);
    }
}

void Scene::update_positions(const std::vector<Point>& points, bool clamp, bool hidden)
{
    unsigned j = 0;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (hidden && vi->is_hidden()) continue;
        
        Point pi = points[j++];
        if (clamp) pi = m_domain.clamp(pi);
        vi->set_position(pi);
    }
}

void Scene::update_weights(const std::vector<FT>& weights, bool hidden)
{    
    unsigned j = 0;
    FT mean = compute_mean(weights);
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (hidden && vi->is_hidden()) continue;
        //vi->set_weight(weights[j++]);
        vi->set_weight(weights[j++] - mean);
    }
}

void Scene::update_singularities(const std::vector<PointSingularity> &pointSingularities, const std::vector<CurveSingularity> &curveSingularities)
{
    m_point_singularities = pointSingularities;
    m_curve_singularities = curveSingularities;
}

void Scene::reset_weights()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->set_weight(0.0);
    }
    update_triangulation();
}

FT Scene::compute_value_integral() const
{
	return m_domain.integrate_intensity();
}

void Scene::pre_build_dual_cells()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        if (vertex->is_hidden()) continue;
        
        bool ok = m_rt.pre_build_polygon(vertex, vertex->dual().points());
        /*
        if (!ok) 
            std::cout << "Vertex " << vertex->get_index() 
            << ": pre_build_dual_cell failed" << std::endl;
        */
    }
}

void Scene::pre_compute_area()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->pre_compute_area();
    }

    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->pre_compute_centroid();
    }
}

FT Scene::integrate_singularities()
{
    FT val = 0.0;

    std::vector<PointSingularity>::iterator it;

    for (it = m_point_singularities.begin(); it != m_point_singularities.end(); it++)
    {
        PointSingularity ps = *it;
        val += ps.get_value();
    }

    return val;
}
