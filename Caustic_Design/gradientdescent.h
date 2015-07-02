#ifndef GRADIENTDESCENT_H
#define GRADIENTDESCENT_H

#include "optimal_transport.h"

class GradientDescent
{
public:
    GradientDescent(OptimalTransport* ot);
    void run();

private:
    OptimalTransport* ot;
};

#endif // GRADIENTDESCENT_H
