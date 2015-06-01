#include "optimal_transport.h"

OptimalTransport::OptimalTransport(Scene* sc, Scene* tc):
    m_scene(sc),
    target_scene(tc)
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
    n = target_weights.size();
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(n);
    lbfgs_parameter_t param;
    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return;
    }
    /* Initialize the variables. Initial weight vector is 0 for first try (TODO: make previous guess if multiscale-approach used) */
    for (i = 0;i < n;i ++) {
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
    /* Report the result. */
    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    lbfgs_free(x);

    std::cout << "done" << std::endl;

/*
    std::cout << "running optimal transport" << std::endl;

    std::vector<FT> weights;
    m_scene->collect_visible_weights(weights);
    std::vector<FT>::iterator iterator;

    int i=0;

    for(iterator = weights.begin();
        iterator != weights.end();
        iterator++)
    {
        if(i == 100)
        {
            weights[i] = 0.005;
        }
        i++;
    }


    m_scene->update_weights(weights);
    m_scene->update_triangulation();
*/
}

lbfgsfloatval_t OptimalTransport::evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
{

    if(!target_scene){
        std::cerr << "target scene not available!" << std::endl;
        return 0;
    }
    std::vector<Vertex_handle> vertices = target_scene->getVertices();

    std::cout << "Eval.. n= " << n << std::endl;
    std::vector<FT> weights = std::vector<FT>(n);
    int i;
    for(i=0; i<n; i++){
        weights[i] = x[i];
    }

    target_scene->update_weights(weights);
    target_scene->update_triangulation();

    std::vector<FT> weight_gradients;
    target_scene->compute_weight_gradient(weight_gradients, 1.0);

    for (i=0; i<n; i++){
        g[i] = weight_gradients[i];
    }

    for (i=0; i<n; i++){
        FT area = vertices[i]->compute_area();
    }

    fx = 0;

    return 0;
    /*
    int i;
    lbfgsfloatval_t fx = 0.0;
    for (i = 0;i < n;i += 2) {
        lbfgsfloatval_t t1 = 1.0 - x[i];
        lbfgsfloatval_t t2 = 10.0 * (x[i+1] - x[i] * x[i]);
        g[i+1] = 20.0 * t2;
        g[i] = -2.0 * (x[i] * g[i+1] + t1);
        fx += t1 * t1 + t2 * t2;
    }
    return fx;*/
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
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}


bool OptimalTransport::prepare_data()
{
    // --- ensure scenes are available
    if(!m_scene) return false;
    if(!target_scene) return false;
    std::cout << "scenes available.. ";

    // --- retrieve points, weights, vertices
    source_points.clear();
    std::vector<FT> scene_weights = std::vector<FT>();
    m_scene->collect_sites(source_points, scene_weights);

    target_points.clear();
    target_weights.clear();
    target_scene->collect_sites(target_points, target_weights);

    source_vertices = m_scene->getVertices();
    target_vertices = target_scene->getVertices();

    // --- ensure they are of same dimension
    //if(target_points.size() != source_points.size()) return false;
    std::cout << "same point amount.. ";
    //if(target_weights.size() != scene_weights.size()) return false;
    std::cout << "same weight amount.. ";
    if(target_vertices.size() != source_vertices.size())
    {
        std::cout << "error.. target_vertices.size = " << target_vertices.size() << " != " << source_vertices.size() << " = source.vertices.size";
        return false;
    }
    std::cout << "same vertex amount.. ";

    // --- every cell has same probability (normalized irradiance per cell)
    // TODO: is this a valid cast?
    probability_per_cell = 1.0 / (double)target_points.size();

    std::cout << std::endl;
    // --- no issue found
    return true;
}
