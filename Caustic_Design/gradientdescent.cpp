#include "gradientdescent.h"

GradientDescent::GradientDescent(OptimalTransport* ot):
    ot(ot)
{
}

void GradientDescent::run()
{

    lbfgsfloatval_t* weights = new lbfgsfloatval_t[SITE_AMOUNT];
    lbfgsfloatval_t* g = new lbfgsfloatval_t[SITE_AMOUNT];

    for (uint i=0; i<SITE_AMOUNT; i++)
    {
        g[i] = 0.0;
        weights[i] = 0.0;
    }

    ot->prepare_data();
    ot->prepare_level_data(weights, SITE_AMOUNT);

    int iteration = 0;
    int continue_optimization = 0;

    do{
        iteration++;
        lbfgsfloatval_t val = ot->evaluate(weights, g, SITE_AMOUNT, 1.0);

        for (uint i=0; i<SITE_AMOUNT; i++)
        {
            weights[i] -= (0.5*g[i]);
        }

        lbfgsfloatval_t gnorm = calc_norm(g, SITE_AMOUNT);
        lbfgsfloatval_t xnorm = calc_norm(weights, SITE_AMOUNT);

        continue_optimization = ot->progress(weights, g, val, xnorm, gnorm, 1.0, SITE_AMOUNT, iteration, 1);
    }while(continue_optimization == 0);



    delete[] weights;
    delete[] g;
}


lbfgsfloatval_t GradientDescent::calc_norm(lbfgsfloatval_t *array, uint n)
{

    lbfgsfloatval_t val = 0.0;

    for (uint i=0; i<n; i++)
    {
        val += (array[i]*array[i]);
    }

    return sqrt(val);
}
