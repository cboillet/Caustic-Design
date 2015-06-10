#include "optimal_transport.h"
#include "scene.h"

OptimalTransport::OptimalTransport(Scene*m_scene, Scene*source_scene, MainWindow* win):
    m_scene(m_scene),
    source_scene(source_scene),
    win(win)
{
}

void OptimalTransport::runOptimalTransport()
{

    if (!prepare_data())
    {
        std::cerr << "Invalid data for optimal transport" << std::endl;
        return;
    }

    std::cout << "running lfbgs" << std::endl;
    int n, i, ret = 0;
    n = source_weights.size();
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(n);
    lbfgs_parameter_t param;
    if (x == NULL)
    {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return;
    }
    /* Initialize the variables. Initial weight vector is 0 for first try (TODO: make previous guess if multiscale-approach used) */
    for (i = 0;i < n;i ++)
    {
        x[i] = 0.0;
    }
    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/
    /*
        Start the L-BFGS optimization; this will invoke the callback functions
        evaluate() and progress() when necessary.
     */
    ret = lbfgs(n, x, &fx, evaluate, progress, this, &param);

    evaluate_results(ret, x, n);

    /* Report the result. */
    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    lbfgs_free(x);

    std::cout << "done" << std::endl;

}

void OptimalTransport::evaluate_results(int ret, lbfgsfloatval_t *x, int n){

    std::string resultString = "unknown";
    switch(ret)
    {
    case LBFGS_SUCCESS:
        resultString = "LBFGS_SUCCESS";
        break;
    case LBFGS_ALREADY_MINIMIZED:
        resultString = "LBFGS_ALREADY_MINIMIZED";
        break;
    case LBFGSERR_UNKNOWNERROR:
        resultString = "LBFGSERR_UNKNOWNERROR";
        break;
    case LBFGSERR_LOGICERROR:
        resultString = "LBFGSERR_LOGICERROR";
        break;
    case LBFGSERR_OUTOFMEMORY:
        resultString = "LBFGSERR_OUTOFMEMORY";
        break;
    case LBFGSERR_CANCELED:
        resultString = "LBFGSERR_CANCELED";
        break;
    case LBFGSERR_INVALID_N:
        resultString = "LBFGSERR_INVALID_N";
        break;
    case LBFGSERR_INVALID_N_SSE:
        resultString = "LBFGSERR_INVALID_N_SSE";
        break;
    case LBFGSERR_INVALID_X_SSE:
        resultString = "LBFGSERR_INVALID_X_SSE";
        break;
    case LBFGSERR_INVALID_EPSILON:
        resultString = "LBFGSERR_INVALID_EPSILON";
        break;
    case LBFGSERR_INVALID_TESTPERIOD:
        resultString = "LBFGSERR_INVALID_TESTPERIOD";
        break;
    case LBFGSERR_INVALID_DELTA:
        resultString = "LBFGSERR_INVALID_DELTA";
        break;
    case LBFGSERR_INVALID_LINESEARCH:
        resultString = "LBFGSERR_INVALID_LINESEARCH";
        break;
    case LBFGSERR_INVALID_MINSTEP:
        resultString = "LBFGSERR_INVALID_MINSTEP";
        break;
    case LBFGSERR_INVALID_MAXSTEP:
        resultString = "LBFGSERR_INVALID_MAXSTEP";
        break;
    case LBFGSERR_INVALID_FTOL:
        resultString = "LBFGSERR_INVALID_FTOL";
        break;
    case LBFGSERR_INVALID_WOLFE:
        resultString = "LBFGSERR_INVALID_WOLFE";
        break;
    case LBFGSERR_INVALID_GTOL:
        resultString = "LBFGSERR_INVALID_GTOL";
        break;
    case LBFGSERR_INVALID_XTOL :
        resultString = "LBFGSERR_INVALID_XTOL";
        break;
    case LBFGSERR_INVALID_MAXLINESEARCH:
        resultString = "LBFGSERR_INVALID_MAXLINESEARCH";
        break;
    case LBFGSERR_INVALID_ORTHANTWISE:
        resultString = "LBFGSERR_INVALID_ORTHANTWISE";
        break;
    case LBFGSERR_INVALID_ORTHANTWISE_START:
        resultString = "LBFGSERR_INVALID_ORTHANTWISE_START";
        break;
    case LBFGSERR_INVALID_ORTHANTWISE_END:
        resultString = "LBFGSERR_INVALID_ORTHANTWISE_END";
        break;
    case LBFGSERR_OUTOFINTERVAL:
        resultString = "LBFGSERR_OUTOFINTERVAL";
        break;
    case LBFGSERR_INCORRECT_TMINMAX:
        resultString = "LBFGSERR_INCORRECT_TMINMAX";
        break;
    case LBFGSERR_ROUNDING_ERROR:
        resultString = "LBFGSERR_ROUNDING_ERROR";
        break;
    case LBFGSERR_MINIMUMSTEP:
        resultString = "LBFGSERR_MINIMUMSTEP";
        break;
    case LBFGSERR_MAXIMUMSTEP:
        resultString = "LBFGSERR_MAXIMUMSTEP";
        break;
    case LBFGSERR_MAXIMUMLINESEARCH:
        resultString = "LBFGSERR_MAXIMUMLINESEARCH";
        break;
    case LBFGSERR_MAXIMUMITERATION:
        resultString = "LBFGSERR_MAXIMUMITERATION";
        break;
    case LBFGSERR_WIDTHTOOSMALL:
        resultString = "LBFGSERR_WIDTHTOOSMALL";
        break;
    case LBFGSERR_INVALIDPARAMETERS:
        resultString = "LBFGSERR_INVALIDPARAMETERS";
        break;
    case LBFGSERR_INCREASEGRADIENT:
        resultString = "LBFGSERR_INCREASEGRADIENT";
        break;
    }


    std::vector<FT> weights = std::vector<FT>(n);
    for (int i=0; i<n; i++)
    {
        weights[i] = x[i];
        if(weights[i] != 0.0)
        {
            std::cout << "weights[" << i << "] = " << weights[i] << std::endl;
        }
    }

    std::cout << "Finished solving, exit status: " << resultString << std::endl;

    std::vector<Point> points = std::vector<Point>();
    std::vector<FT> w = std::vector<FT>();
    //source_scene->collect_visible_points(points);
    //source_scene->collect_sites(points, w);
    m_scene->update_positions(source_points);
    m_scene->update_weights(weights);
    m_scene->update_triangulation();
}

/*
 * This is the main callback method for liblbfgs.
 * Within here, the new value of the weights (in *x) are given.
 * The gradient needs to be computed from it as well as the new value of the function (fx)
 */
lbfgsfloatval_t OptimalTransport::evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
{

    std::cout << "Eval.. step = " << step << std::endl;
    std::vector<FT> weights = std::vector<FT>(n);
    int i;
    for(i=0; i<n; i++)
    {
        weights[i] = x[i];
    }

    source_scene->update_positions(source_points, false);
    source_scene->update_weights(weights, false);
    source_scene->update_triangulation();

    // TODO: following if clause is only for debugging. should be removed later
    std::vector<Vertex_handle> target_vertices = source_scene->getVertices();
    int nv = target_vertices.size();
    if(n != nv)
    {
        m_scene->update_positions(source_points);
        m_scene->update_weights(weights);
        m_scene->update_triangulation();
        win->update();
        std::cerr << "wrong amount of vertices" << std::endl;
        return 0;
    }

    // fx = f(w)
    lbfgsfloatval_t fx = 0.0;

    // f(w) = --- the convex function to be minimized
    //std::cout << "n = " << n << ", source_vertices.size() = " << nv << std::endl;
    //std::cout << std::flush;
    for(int i=0; i<n; i++)
    {
        FT integration_term = target_vertices[i]->compute_wasserstein( m_points[i], x[i] );
        //FT integration_term = m_vertices[i]->compute_wasserstein( m_points[i], x[i] );
        fx += ( x[i]* source_capacities[i] - integration_term );
    }

    // df/dwi = --- The derivate of the convex function
    for (i=0; i<n; i++)
    {
        g[i] = source_capacities[i] - m_vertices[i]->compute_area() / integrated_m_intensity;
    }

    return fx;
}

int OptimalTransport::progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
{

    return 1;

    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

/*
 * Checks if input is valid and makes all calculations for data that doesn't change during optimization.
 */
bool OptimalTransport::prepare_data()
{
    // --- ensure scenes are available
    if(!m_scene) return false;
    if(!source_scene) return false;
    std::cout << "scenes available.. ";

    // --- retrieve points, weights, vertices
    m_points.clear();
    std::vector<FT> scene_weights = std::vector<FT>();
    m_scene->collect_sites(m_points, scene_weights);
    // TODO: remove following lines
    /*scene_weights[16] = 0.01;
    scene_weights[60] = 0.01;
    scene_weights[118] = 0.01;
    m_scene->update_weights(scene_weights);
    m_scene->update_triangulation();
    return false;
    // --- until here*/

    source_points.clear();
    source_weights.clear();
    source_scene->collect_sites(source_points, source_weights);

    m_vertices = m_scene->getVertices();
    source_vertices = source_scene->getVertices();

    integrated_m_intensity = m_scene->getDomain().integrate_intensity();
    integrated_source_intensity = source_scene->getDomain().integrate_intensity();

    for (int i=0; i< m_vertices.size(); i++)
    {
        source_capacities.push_back(source_vertices[i]->compute_area() / integrated_source_intensity);
    }

    // --- ensure they are of same dimension
    //if(source_points.size() != m_points.size()) return false;
    std::cout << "same point amount.. ";
    //if(source_weights.size() != scene_weights.size()) return false;
    std::cout << "same weight amount.. ";
    if(source_vertices.size() != m_vertices.size())
    {
        std::cout << "error.. source_vertices.size = " << source_vertices.size() << " != " << m_vertices.size() << " = m_vertices.size";
        return false;
    }
    std::cout << "same vertex amount.. ";

    std::cout << std::endl;
    // --- no issue found
    return true;
}
