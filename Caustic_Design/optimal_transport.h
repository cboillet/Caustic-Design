#ifndef OPTIMAL_TRANSPORT_H
#define OPTIMAL_TRANSPORT_H

// System
#include "lbfgs.h"
#include <string>

// Local
#include "scene.h"
#include "window.h"
#include "domain.h"
#include "random.h"

class OptimalTransport
{
public:
    OptimalTransport(Scene*m_scene, Scene*source_scene, MainWindow* win);
    void runOptimalTransport();
    void evaluate_results(int ret, lbfgsfloatval_t *x, int n);

    Scene* m_scene;
    Scene*source_scene;
    MainWindow* win;


protected:
    // L_BFGS

    /*
     *  Callback interface to provide objective function and gradient evaluations.
     *
     *  The lbfgs() function call this function to obtain the values of objective function and its gradients when needed.
     *  A client program must implement this function to evaluate the values of the objective function and its gradients,
     *  given current values of variables.
     */
    static lbfgsfloatval_t evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<OptimalTransport*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
            const lbfgsfloatval_t *x,
            lbfgsfloatval_t *g,
            const int n,
            const lbfgsfloatval_t step
            );

    /*
     * Callback interface to receive the progress of the optimization process.
     *
     * The lbfgs() function call this function for each iteration.
     * Implementing this function, a client program can store or display the current progress of the optimization process.
     */
    static int progress(
        void *instance,
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
        return reinterpret_cast<OptimalTransport*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        );


private:

    std::vector<Vertex_handle> current_source_vertices;
    std::vector<FT> source_weights;
    std::vector<Point> source_points;

    std::vector<Vertex_handle> m_vertices;
    std::vector<Point> m_points;
    std::vector<FT> initial_source_capacities;

    FT integrated_m_intensity;
    FT integrated_source_intensity;

    bool prepare_data();
};

#endif // OPTIMAL_TRANSPORT_H
