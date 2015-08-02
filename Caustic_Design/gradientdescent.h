#ifndef GRADIENTDESCENT_H
#define GRADIENTDESCENT_H

#include "optimal_transport.h"
#include "config.h"
#include "lbfgs.h"
#include "math.h"

class GradientDescent
{
public:
    GradientDescent(OptimalTransport* ot);
    void run(lbfgsfloatval_t* weights);
    lbfgsfloatval_t calc_norm(lbfgsfloatval_t* array, uint n);

private:
    OptimalTransport* ot;
};

#endif // GRADIENTDESCENT_H
