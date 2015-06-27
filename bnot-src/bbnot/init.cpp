#include "scene.h"
#include "timer.h"

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

void Scene::clear_triangulation()
{
    m_capacities.clear();
    m_vertices.clear();
    m_rt.clear();
}

void Scene::update_triangulation()
{
    std::vector<FT> weights;
    std::vector<Point> points;
    collect_sites(points, weights);
    construct_triangulation(points, weights);
}

void Scene::construct_triangulation(const std::vector<Point>& points,
                                    const std::vector<FT>& weights)
{
    if (m_timer_on)
    {
        Timer::start_timer(m_timer, COLOR_BLUE, "Triangulation");
        std::cout << std::endl;
    }
    
    clear_triangulation();
    populate_vertices(points, weights);
    m_rt.pre_compute_cells();
    compute_capacities(m_capacities);

    if (m_timer_on)
        Timer::stop_timer(m_timer, COLOR_BLUE);
}

void Scene::populate_vertices(const std::vector<Point>& points,
                              const std::vector<FT>& weights)
{
    if (m_timer_on)
        Timer::start_timer(m_timer, COLOR_YELLOW, "Populate");

    unsigned nb = 0;
    unsigned nsites = points.size();
    for (unsigned i = 0; i < nsites; ++i)
    {
        Vertex_handle vertex = insert_vertex(points[i], weights[i], nb);
        if (vertex == Vertex_handle()) continue;
        m_vertices.push_back(vertex);
        nb++;
    }
    
    if (m_timer_on)
        Timer::stop_timer(m_timer, COLOR_YELLOW);
}

Vertex_handle Scene::insert_vertex(const Point& point,
                                   const FT weight,
                                   const unsigned index)
{
    Weighted_point wp(point, weight);
    Vertex_handle vertex = m_rt.insert(wp);

    if (vertex->get_index() != -1)
        return Vertex_handle();
    
    vertex->set_index(index);
    return vertex;
}

FT Scene::compute_mean_capacity() const
{
    FT domain_area = m_domain.get_area();
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

void Scene::reset_weights()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->set_weight(0.0);
    }
    update_triangulation();
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
        vi->set_weight(weights[j++] - mean);
    }
}

void Scene::project_positions_to_domain()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        
        Point pi = vi->get_position();
        if (m_domain.is_outside(pi) || m_rt.is_boundary(vi))
            pi = m_domain.project(pi);
        vi->set_position(pi);
    }
}
