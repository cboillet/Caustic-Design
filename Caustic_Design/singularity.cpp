#include "singularity.h"

// -------------------
//  Point Singularity
// -------------------

PSingularity::PSingularity():
    m_value(0.0)
{
}

PSingularity::PSingularity(FT value):
    m_value(value)
{
}

PSingularity::PSingularity(Point & position, FT value):
    m_position(position),
    m_value(value)
{
}


// -------------------
//  Curve Singularity
// -------------------

CurveSingularity::CurveSingularity()
{
}
