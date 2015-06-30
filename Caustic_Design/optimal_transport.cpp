#include "optimal_transport.h"
#include "config.h"
#include "scene.h"

OptimalTransport::OptimalTransport(Scene*m_scene, Scene*source_scene, MainWindow* win, GlViewer* source_viewer):
    m_scene(m_scene),
    source_scene(source_scene),
    win(win),
    source_viewer(source_viewer)
{
    level_max = LEVEL_MAX;
    site_amount = SITE_AMOUNT;
}

void OptimalTransport::runOptimalTransport()
{

    if (!prepare_data())
    {
        std::cerr << "Invalid data for optimal transport" << std::endl;
        return;
    }

    for (current_level = (level_max-1); current_level >= 0; current_level--)
    {
        unsigned site_amount = get_level_sites(current_level);

        std::cout << "running lfbgs for level " << current_level << std::endl;
        int n, i, ret = 0;
        n = site_amount;
        lbfgsfloatval_t fx;
        lbfgsfloatval_t *x = lbfgs_malloc(n);
        lbfgs_parameter_t param;
        if (x == NULL)
        {
            printf("ERROR: Failed to allocate a memory block for variables.\n");
            return;
        }
        prepare_level_data(x, n);


        /* Initialize the variables. Initial weight vector is 0 for first try (TODO: make previous guess if multiscale-approach used) */
        /*for (i = 0;i < n;i ++)
        {
            x[i] = 0.0;
        }*/
        /* Initialize the parameters for the L-BFGS optimization. */
        lbfgs_parameter_init(&param);
        param.linesearch = LBFGS_LINESEARCH_MORETHUENTE;
        //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
        //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
        param.m = 10;
        param.max_linesearch = 40;
        //param.ftol = 0.000000000000001;
        param.epsilon = 0;
        //param.delta = 0.0001;
        /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/
        /*
            Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.
         */
        ret = lbfgs(n, x, &fx, evaluate, progress, this, &param);

        if(!evaluate_results(ret, x, n))
        {
            std::cerr << "error in level " << current_level << ", stop." << std::endl;
            lbfgs_free(x);
            return;
        }


        /* Report the result. */
        printf("L-BFGS optimization on level %d terminated with status code = %d\n", current_level, ret);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        lbfgs_free(x);


        std::cout << "done" << std::endl;
    }

    clean();

}


bool OptimalTransport::evaluate_results(int ret, lbfgsfloatval_t *x, int n)
{

    std::string resultString = get_result_string(ret);

    std::vector<FT> weights = std::vector<FT>(n);
    for (int i=0; i<n; i++)
    {
        weights[i] = x[i];
    }

    std::cout << "cleaning up" << std::endl;

    std::cout << "Finished solving, exit status: " << resultString << std::endl;

    // assign value to original scene (optional, can be removed!)
    scaled_scenes[current_level]->construct_triangulation(source_points, weights);

    return ret >= 0;
}

std::string OptimalTransport::get_result_string(int ret)
{
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

    return resultString;
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
    FT min_weight = 1000;
    FT max_weight = -1000;
    int i;
    for(i=0; i<n; i++)
    {
        weights[i] = x[i];

        if(weights[i] < min_weight)
            min_weight = weights[i];

        if(weights[i] > max_weight)
            max_weight = weights[i];
    }

    std::cout << "min-weight = " << min_weight << ", max-weight = " << max_weight << std::endl;

    // --- update the triangulation with the old points and the new weights
    //source_scene->update_weights(weights, false);
    //source_scene->update_triangulation();
    scaled_scenes[current_level]->construct_triangulation(source_points, weights);
    current_source_vertices = scaled_scenes[current_level]->getVertices();
    // --- update UI (can be removed for improved performance)
#ifdef LIVE_DEMO
    win->update();
#endif

    // fx = f(w)
    FT fx = 0.0;

    // f(w) = --- the convex function to be minimized
    FT integral_sum = 0.0;
    FT source_sum = 0.0;
    for(int i=0; i<n; i++)
    {

        FT integration_term = current_source_vertices[i]->compute_wasserstein( x[i], integrated_m_intensity );
        integral_sum += integration_term;
        source_sum += x[i] * initial_source_capacity;
        fx += (x[i] * (initial_source_capacity) - integration_term);

        g[i] = ( current_source_vertices[i]->compute_area() / integrated_m_intensity ) - initial_source_capacity;

        //std::cout << "fx += " << x[i] << " * " << initial_source_capacities[i] << " - " << integration_term << std::endl;
    }

    std::cout << "integrated sum = " << integral_sum << std::endl;
    std::cout << "source_sum = " << source_sum << std::endl;

    //source_scene->compute_capacities(areas);
    // df/dwi = --- The derivative of the convex function
    //for (i=0; i<n; i++)
    //{
        //g[i] = 2*(( current_source_vertices[i]->compute_area() / integrated_m_intensity ) - initial_source_capacities[i]) * (x[i] * initial_source_capacities[i] - integration_term );
        //std::cout << "g[" << i << "] = " << g[i] << " = " << initial_source_capacities[i] << " - " << areas[i] << std::endl;
    //}

    std::cout << "evaluate done, fx = " << fx << std::endl;

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

    bool hidden_vertices = false;

    double epsilon = current_level == 0 ? 3e-4 : 3e-3;
    bool norm = gnorm < epsilon;
    int hidden_vertices_amount = 0;

    std::vector<Vertex_handle> vertices = scaled_scenes[current_level]->getVertices();

    for (int i = 0; i<vertices.size(); i++)
    {
        if(vertices[i]->is_hidden())
        {
            hidden_vertices_amount++;
        }
    }

    hidden_vertices = (hidden_vertices_amount != 0);

    bool will_stop = norm & !hidden_vertices;

    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("  hidden vertices: %d, norm reached: %s, will stop: %s", hidden_vertices_amount, norm ? "true" : "false", will_stop ? "true" : "false");
    printf("\n");



    // returning 0 will continue, returning sth else will stop
    if(will_stop)
        return 1;

    return 0;
}

/*
 * Checks if input is valid and makes all calculations for data that doesn't change during optimization.
 */
bool OptimalTransport::prepare_data()
{
    // --- ensure scenes are available
    if(!m_scene)
    {
        std::cerr << "m_scene not available" << std::endl;
        return false;
    }
    if(!source_scene)
    {
        std::cerr << "source_scene not available" << std::endl;
        return false;
    }

    std::cout << "scenes available.. ";

    VoronoiCreator voronoi_creator;// = VoronoiCreator();
    QString source_image = source_scene->getDomain().get_filename();

    scaled_scenes = new Scene*[level_max];

    for(unsigned i=0; i<level_max; i++)
    {

        scaled_scenes[i] = i == 0 ? source_scene : new Scene();
        unsigned scene_sites = get_level_sites(i);

        scaled_scenes[i]->load_image(source_image);

#ifdef LIVE_DEMO
        source_viewer->set_scene(scaled_scenes[i]);
#endif

        voronoi_creator.generate_voronoi(scaled_scenes[i], scene_sites, EPSILON, source_viewer);
        //voronoi_creator.init_points(scene_sites, scaled_scenes[i]);

        /*for(int j=0; j<10; j++)
        {
            voronoi_creator.apply_lloyd_optimization(scaled_scenes[i]);
#ifdef LIVE_DEMO
            win->update();
#endif
        }*/

        std::cout << "site-amount for level " << i << ": " << scene_sites << std::endl;
    }


    // --- retrieve points, weights, vertices
    m_points.clear();
    std::vector<FT> scene_weights = std::vector<FT>();
    m_scene->collect_sites(m_points, scene_weights);

    m_vertices = m_scene->getVertices();
    current_source_vertices = source_scene->getVertices();

    source_scene->collect_sites(source_points, source_weights);

    // --- integrate the intensities (areas)
    integrated_m_intensity = m_scene->getDomain().integrate_intensity();
    integrated_m_intensity += m_scene->integrate_singularities();
    // source does not have singularities
    integrated_source_intensity = source_scene->getDomain().integrate_intensity();

    return true;
/*
    // --- ensure they are of same dimension
    if(source_points.size() != m_points.size())
    {
        std::cerr << "source_points and m_points are not of same length" << std::endl;
        return false;
    }
    if(source_weights.size() != scene_weights.size())
    {
        std::cerr << "source_points and m_points are not of same length" << std::endl;
        return false;
    }
    if(current_source_vertices.size() != m_vertices.size())
    {
        std::cerr << "error.. source_vertices.size = " << current_source_vertices.size() << " != " << m_vertices.size() << " = m_vertices.size" << std::endl;
        return false;
    }
    std::cout << "same vertex amount.. ";

    std::cout << std::endl;
    // --- no issue found
    return true;*/
}

void OptimalTransport::prepare_level_data(lbfgsfloatval_t *initial_weights, unsigned n)
{

    source_viewer->set_scene(scaled_scenes[current_level]);
    source_points.clear();
    source_weights.clear();
    scaled_scenes[current_level]->collect_sites(source_points, source_weights);

    current_source_vertices = scaled_scenes[current_level]->getVertices();


    // --- pre-compute the not-changing value of capacities (== probabilities)
    initial_source_capacity = FT(1.0 / ((FT)current_source_vertices.size()));

    // --- load the target image into the source scene.
    // --- reinsert the points, update triangulation
    // --- points need to be re-inserted and re-read to have valid points
    QString filename = m_scene->getDomain().get_filename();
    std::cout << "loading image " << filename.toStdString() << "...";
    scaled_scenes[current_level]->load_image(filename);
    std::cout << "re-inserting points...";
    source_weights.clear();
    source_points.clear();
    // load the singularities of the target image
    std::vector<PointSingularity> ps;
    std::vector<CurveSingularity> cs;
    m_scene->collect_singularities(ps, cs);
    scaled_scenes[current_level]->update_singularities(ps, cs);
    scaled_scenes[current_level]->update_positions(source_points);
    scaled_scenes[current_level]->update_triangulation();
    scaled_scenes[current_level]->collect_sites(source_points, source_weights);

    // copy points to ensure we don't only have pointers -- possibly can be removed
    for(int i=0; i<source_points.size(); i++)
    {
        source_points[i] = Point(source_points[i]);
    }


    unsigned previous_level = current_level+1;
    for (unsigned i=0; i<n; i++)
    {
        if(current_level == (level_max-1))
        {
            initial_weights[i] = 0.0;
        }
        else
        {
            initial_weights[i] = get_initial_weight(source_points[i], scaled_scenes[previous_level]);
        }
    }
}

// point is the point we need to get the weight for (site of current cell)
// scene is the higher-level scene (with less granularity) where we have already calculated the weights
FT OptimalTransport::get_initial_weight(Point point, Scene *scene)
{
    std::vector<Vertex_handle> vertices = scene->getVertices();

    // find vertex with minimum distance to point
    // TODO check for max distance
    FT minimum_distance = 10000000;
    Vertex_handle minimum_vertex;
    for(std::vector<Vertex_handle>::iterator it = vertices.begin();
        it != vertices.end();
        it++
        )
    {
        Vertex_handle current = *it;
        FT distance = CGAL::squared_distance(point, current->get_position());
        if(distance < minimum_distance){
            minimum_distance = distance;
            minimum_vertex = current;
        }
    }

    // return weight of the vertex we found
    return minimum_vertex->get_weight();
}

unsigned OptimalTransport::get_level_sites(unsigned level)
{
    return site_amount / pow(5, level);
    //return m_scene->getVertices().size() / pow(5, level);
}

void OptimalTransport::clean()
{
    for(unsigned i=1; i<level_max; i++){
        delete scaled_scenes[i];
    }

    delete scaled_scenes;
}
