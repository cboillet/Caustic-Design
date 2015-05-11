#ifndef _UTIL_H_
#define _UTIL_H_

#define MINN -1e100
#define MAXN  1e100
#define EPS   1e-12

#include <cmath>
#include <vector>

template <class T>
T compute_mean(const std::vector<T>& data)
{
    T mean = 0.0;
    for (unsigned i = 0; i < data.size(); ++i)
        mean += data[i];
    return (mean / data.size());
}

template <class T>
void remove_mean(std::vector<T>& data)
{
    T mean = compute_mean(data);
    for (unsigned i = 0; i < data.size(); ++i)
        data[i] -= mean;
}

template <class T>
double compute_norm(std::vector<T>& data)
{
    double norm2 = 0.0;
    for (unsigned i = 0; i < data.size(); ++i)
        norm2 += data[i]*data[i];
    return std::sqrt(norm2);
}

template <class Point, class OutputIterator>
void generate_regular_polygon(unsigned nb, 
                              const double radius,
                              const Point center,
                              OutputIterator out)
{
    if (nb < 3) nb = 3;
    double step = 2*M_PI / nb;
    double start = 0.25*M_PI;
    for (unsigned i = 0; i < nb; ++i)
    {
        double angle = start + i*step;
        Point q(radius*cos(angle) + center.x(), 
                radius*sin(angle) + center.y());
        *out++ = q;
    }
}

#endif
