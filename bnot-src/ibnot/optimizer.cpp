#include <map>

#include "scene.h"
#include "util.h"
#include "timer.h"
#include "pw_line_search.h"

#include "matrix/sparse_array.h"
#include "matrix/suite_sparse_qr.h"

typedef CLSWeights<Scene, FT> LSWeights;

typedef CLSPositions<Scene, Point, Vector> LSPositions;

FT Scene::optimize_positions_via_lloyd(bool update)
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "Centroid");
    std::vector<Point> points;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];        
        if (vi->is_hidden()) continue;
        Point ci = vi->compute_centroid();
        points.push_back(ci);
    }
    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);

    update_positions(points);
    if (update) update_triangulation();
    
    std::vector<Vector> gradient;
    compute_position_gradient(gradient);
    return compute_norm(gradient);
}

FT Scene::optimize_positions_via_gradient_ascent(FT timestep, bool update)
{
    std::vector<Point> points;
    collect_visible_points(points);

    std::vector<Vector> gradient;
    compute_position_gradient(gradient);
    
    if (timestep <= 0.0)
    {
        double mean_capacity = compute_mean(m_capacities);
        double max_alpha = 1.0 / mean_capacity;
        LSPositions line_search(this, 10, max_alpha);        
        timestep = line_search.run_bt(points, gradient);
    } else {    
        for (unsigned i = 0; i < points.size(); ++i)
        {
            Point  pi = points[i];
            Vector gi = gradient[i];
            points[i] = pi + timestep*gi;
        }
        update_positions(points);
        if (update) update_triangulation();
    }
    
    compute_position_gradient(gradient);
    return compute_norm(gradient);
}

FT Scene::optimize_weights_via_gradient_descent(FT timestep, bool update)
{
    std::vector<FT> gradient;
    compute_weight_gradient(gradient, -1.0);
    
    std::vector<FT> weights;
    collect_visible_weights(weights);
    
    if (timestep <= 0.0)
    {
        LSWeights line_search(this, 10, 2.0);
        timestep = line_search.run_bt(weights, gradient);
    } else {
        for (unsigned i = 0; i < weights.size(); ++i)
        {
            FT wi = weights[i];
            FT gi = gradient[i];
            weights[i] = wi + timestep*gi;
        }
        update_weights(weights);
        if (update) update_triangulation();
    }

    compute_weight_gradient(gradient);
    return compute_norm(gradient);
}

FT Scene::optimize_weights_via_newton(FT timestep, bool update)
{
    std::vector<FT> gradient;
    compute_weight_gradient(gradient, -1.0);
    
    std::vector<FT> direction;
    bool ok = solve_newton_step(gradient, direction);
    if (!ok) return 0.0;   
    
    std::vector<FT> weights;
    collect_visible_weights(weights);
    
    if (timestep <= 0.0)
    {
        LSWeights line_search(this, 20, 2.0);
        timestep = line_search.run_bt(weights, direction);
    } else {
        for (unsigned i = 0; i < weights.size(); ++i)
        {
            FT wi = weights[i];
            FT gi = direction[i];
            weights[i] = wi + timestep*gi;
        }
        update_weights(weights);
        if (update) update_triangulation();
    }
    
    compute_weight_gradient(gradient);
    return compute_norm(gradient);
}

bool Scene::solve_newton_step(const std::vector<FT>& b, std::vector<FT>& x)
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "LinearSolver");

    unsigned nb = 0;
    std::map<unsigned, unsigned> indices;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;        
        indices[vi->get_index()] = nb++;
    }
    
    SparseMatrix L(nb, nb);
    build_laplacian(0.5, indices, L);
    
    bool ok = solve_linear_system(L, x, b);
    if (!ok) 
    {
        std::cout << red << "linear solver failed" << white << std::endl;
        return false;
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
    return true;
}

void Scene::build_laplacian(const FT scale,
                            const std::map<unsigned, unsigned>& indices,
                            SparseMatrix& A) const
{
    unsigned nb = A.numRows();
    for (unsigned k = 0; k < m_vertices.size(); ++k)
    {
        Vertex_handle vi = m_vertices[k];
        if (vi->is_hidden()) continue;
        unsigned i = indices.find(vi->get_index())->second;
        
        double diagi = 0.0;
        SparseArray rowi(nb);        
        Edge_circulator ecirc = m_rt.incident_edges(vi);
        Edge_circulator eend  = ecirc;
        CGAL_For_all(ecirc, eend)
        {
            Edge edge = *ecirc;
            if (!m_rt.is_inside(edge)) continue;
            
            Vertex_handle vj = m_rt.get_source(edge);
            if (vj == vi) vj = m_rt.get_target(edge);
            
            unsigned j = vj->get_index();
            j = indices.find(j)->second;
            
            double coef = scale * get_ratio(edge);
            if (std::abs(coef) < EPS) continue;
            
            rowi.setValue(j, -coef);
            diagi += coef;
        }
        
        rowi.setValue(i, diagi);
        A.setRow(i, rowi);
    }
}

bool Scene::solve_linear_system(const SparseMatrix& A,
                                std::vector<double>& x,
                                const std::vector<double>& b) const
{
    SuiteSparseQRFactorizer solver;
    bool ok = solver.factorize(A);
    if (!ok) return false;
    
    ok = solver.solve(b, x);
    return ok;
}

unsigned Scene::optimize_weights_via_gradient_descent_until_converge(FT timestep, 
                                                                     FT threshold,
                                                                     unsigned update,
                                                                     unsigned max_iters)
{
    for (unsigned i = 0; i < max_iters; ++i)
    {
        bool flag = (update == 0 || (i+1) % update == 0);
        FT norm = optimize_weights_via_gradient_descent(timestep, flag);
        if (norm < threshold) return i;
    }
    return max_iters;
}

unsigned Scene::optimize_weights_via_newton_until_converge(FT timestep, 
                                                           FT threshold,
                                                           unsigned update,
                                                           unsigned max_iters)
{
    for (unsigned i = 0; i < max_iters; ++i)
    {
        bool flag = (update == 0 || (i+1) % update == 0);
        FT norm = optimize_weights_via_newton(timestep, flag);
        if (norm < threshold) return i;
    }
    return max_iters;
}

unsigned Scene::optimize_all(FT wstep, FT xstep, unsigned max_newton_iters,
                             FT epsilon, unsigned max_iters,
                             std::ostream& out)
{
    bool global_connectivity = m_fixed_connectivity;
    unsigned nb0 = count_visible_sites();
    
    FT xthreshold = compute_position_threshold(epsilon);
    FT wthreshold = compute_weight_threshold(epsilon);
    
    out << "NbSites = " << nb0 << std::endl;
    out << "Threshold: " << xthreshold << " ; " << wthreshold << std::endl;

    m_fixed_connectivity = false;
    FT coarse_xthreshold = 2.0*xthreshold;
    FT coarse_wthreshold = 2.0*wthreshold;
    
    unsigned iters = 0;    
    unsigned nb_assign = 0;
    while (iters < max_iters)
    {
        iters++;
        reset_weights();        
        nb_assign += optimize_weights_via_newton_until_converge(wstep, coarse_wthreshold, 0, max_newton_iters);        
        FT norm = optimize_positions_via_lloyd(true);
        nb_assign++;
        out << "(Coarse) Norm: " << norm << std::endl;        
        if (norm <= coarse_xthreshold) break;
    }
    
    out << "Partial: " << iters << " iters" << std::endl;
    m_fixed_connectivity = global_connectivity;
    if (iters == max_iters) return iters;
    
    m_fixed_connectivity = false;
    FT fine_xthreshold = xthreshold;
    FT fine_wthreshold = wthreshold;

    while (iters < max_iters)
    {
        iters++;
        unsigned nb1 = count_visible_sites();
        if (nb1 != nb0) reset_weights();
        nb_assign += optimize_weights_via_newton_until_converge(wstep, fine_wthreshold, 0, max_newton_iters);
        FT norm = optimize_positions_via_gradient_ascent(xstep, true);
        nb_assign++;
        out << "(Fine) Norm: " << norm << std::endl;
        if (norm <= fine_xthreshold) break;
    }
    optimize_weights_via_newton_until_converge(wstep, 0.1*fine_wthreshold, 0, max_newton_iters);
    
    std::cout << "NbAssign: " << nb_assign << std::endl;
    
    m_fixed_connectivity = global_connectivity;
    return iters;
}
