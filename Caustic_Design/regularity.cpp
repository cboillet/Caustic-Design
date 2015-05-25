#include "scene.h"
#include "random.h"

FT Scene::compute_regularity_threshold() const
{
    FT m = compute_mean_capacity();
    return m_tau*0.25*m*m;
}

void Scene::detect_and_break_regularity(FT max_radius, unsigned max_teleport)
{
    FT threshold = compute_regularity_threshold();
    FT m = compute_mean_capacity();
    max_radius *= m;
    
    std::vector<FT> variance;
    compute_variance_vector(variance);
    
    std::vector<FT> regularity;
    compute_regularity_vector(variance, regularity);
    
    std::vector<Vertex_handle> blocks;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        if (regularity[i] < threshold) 
            blocks.push_back(vi);
    }
    std::cout << blocks.size() << " blocks detected" << std::endl;
    if (blocks.empty()) return;
    
    std::set<Vertex_handle> for_jittering;
    if (max_radius > 0.0)
    {
        for (unsigned i = 0; i < blocks.size(); ++i)
        {
            Vertex_handle vi = blocks[i];
            for_jittering.insert(vi);
            
            Vertex_circulator vcirc = m_rt.incident_vertices(vi);
            Vertex_circulator vend  = vcirc;
            CGAL_For_all(vcirc, vend)
            {
                Vertex_handle vj = vcirc;
                if (m_rt.is_infinite(vj)) continue;
                for_jittering.insert(vj);
            }
        }
        jitter_vertices(for_jittering, max_radius);
    }
    std::cout << for_jittering.size() << " vertices jittered" << std::endl;
    
    std::set<Vertex_handle> to_teleport;
    if (blocks.size() < max_teleport) 
        max_teleport = blocks.size();
    if (blocks.size() == max_teleport) 
    {
        to_teleport.insert(blocks.begin(), blocks.end());
    } else {
        while (to_teleport.size() < max_teleport)
        {
            unsigned index = random_int(0, blocks.size()-1);
            to_teleport.insert(blocks[index]);
        }
    }
    
    double dx = m_domain.get_dx();
    double dy = m_domain.get_dy();
    for (std::set<Vertex_handle>::iterator it = to_teleport.begin(); it != to_teleport.end(); ++it)
    {
        Vertex_handle vertex = *it;
        double x = random_double(-dx, dx);
        double y = random_double(-dy, dy);
        vertex->set_position(Point(x, y));
    }
    std::cout << max_teleport << " vertices teleported" << std::endl;
    
    reset_weights(); // it calls update
}

FT Scene::compute_max_regularity() const
{
    std::vector<FT> variance;
    compute_variance_vector(variance);
    
    std::vector<FT> regularity;
    compute_regularity_vector(variance, regularity);
    
    FT max_value = 0.0;
    for (unsigned i = 0; i < regularity.size(); ++i)
        max_value = std::max(max_value, regularity[i]);
    return max_value;
}

void Scene::compute_variance_vector(std::vector<FT>& variance) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        FT V = 0.0;
        Vertex_handle vertex = m_vertices[i];        
        if (!vertex->is_hidden()) V = vertex->compute_variance();
        variance.push_back(V);
    }
}

void Scene::compute_regularity_vector(const std::vector<FT>& variance,
                                      std::vector<FT>& regularity) const
{
    // regularity = normalized absolute deviation
    FT C = compute_mean_capacity();
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        FT R = 0.0;
        Vertex_handle vertex = m_vertices[i];
        if (!vertex->is_hidden()) R = compute_regularity(vertex, variance);
        regularity.push_back(R/C);
    }
}

FT Scene::compute_regularity(Vertex_handle vi, const std::vector<FT>& variance) const
{
    FT deviation = 0.0;
    unsigned degree = 0;
    FT Vi = variance[vi->get_index()];
    Vertex_circulator vcirc = m_rt.incident_vertices(vi);
    Vertex_circulator vend  = vcirc;
    CGAL_For_all(vcirc, vend)
    {
        Vertex_handle vj = vcirc;
        if (m_rt.is_infinite(vj)) continue;
        FT Vj = variance[vj->get_index()];
        deviation += std::abs(Vi - Vj);
        degree++;
    }
    return (deviation / degree);
}

void Scene::jitter_vertices(const std::set<Vertex_handle>& vertices, const FT max_radius)
{
    for (std::set<Vertex_handle>::const_iterator it = vertices.begin(); it != vertices.end(); ++it)
    {
        Vertex_handle vertex = *it;
        Point p = vertex->get_position();
        p = jitter_point(p, max_radius);
        vertex->set_position(p);
    }
}

Point Scene::jitter_point(const Point& p, const FT max_radius) const
{
    FT angle = random_double(0.0, 2.0*M_PI);
    FT radius = random_double(0.0, max_radius);
    Vector d(radius*cos(angle), radius*sin(angle));
    return m_domain.clamp(p + d);
}

void Scene::count_sites_per_bin(unsigned N) const
{
    std::vector<FT> bin(N);
    FT xmin = -m_domain.get_dx();
    FT xmax =  m_domain.get_dx();
    FT step = (xmax - xmin) / FT(N);
    for (unsigned j = 0; j < N; ++j)
        bin[j] = xmin + (j+1)*step;
        
    unsigned nb = 0;
    std::vector<FT> counter(N, 0.0);
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        
        nb++;
        const Point& pi = vi->get_position();
        for (unsigned j = 0; j < N; ++j)
        {
            if (pi.x() < bin[j])
            {
                counter[j] += 1.0;
                break;
            }
        }
    }
    if (nb != 0)
    {
        for (unsigned j = 0; j < N; ++j)
            counter[j] /= FT(nb);
    }
    
    std::cout << "Counter: ";
    FT total = 0.0;
    for (unsigned j = 0; j < N; ++j)
    {
        std::cout << counter[j] << " ; ";
        total += counter[j];
    }
    std::cout << total << std::endl;
    
    step = FT(m_domain.get_width()) / FT(N);
    for (unsigned k = 0; k < N; ++k)
        bin[k] = (k+1)*step;

    FT total_mass = 0.0;
    std::vector<FT> mass(N, 0.0);
    std::vector<unsigned> nb_pixel(N, 0);
    for (unsigned i = 0; i < m_domain.get_width(); ++i)
    {
        for (unsigned j = 0; j < m_domain.get_height(); ++j)
        {                
            FT mij = m_domain.get_value(i, j);
            total_mass += mij;
            for (unsigned k = 0; k < N; ++k)
            {
                if (FT(i) < bin[k])
                {
                    nb_pixel[k] ++;
                    mass[k] += mij;
                    break;
                }
            }
        }
    }
    for (unsigned k = 0; k < N; ++k)
        mass[k] /= total_mass;
    
    std::cout << "Ref: ";
    total_mass = 0.0;
    for (unsigned k = 0; k < N; ++k)
    {
        std::cout << mass[k] << " ; ";
        total_mass += mass[k];
    }
    std::cout << total_mass << std::endl;

    std::cout << "NbPixel: ";
    for (unsigned k = 0; k < N; ++k)
        std::cout << nb_pixel[k] << " ; ";
    std::cout << std::endl;
}
