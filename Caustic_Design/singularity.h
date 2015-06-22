#ifndef SINGULARITY_H
#define SINGULARITY_H

template <class Kernel>
class PSingularity
{
public:
    typedef typename Kernel::FT       FT;
    typedef typename Kernel::Point_2  Point;

    PSingularity();
    PSingularity(FT value);
    PSingularity(Point &position, FT value);

    const FT get_value(){ return m_value; }
    const Point& get_position(){ return m_position; }

private:
    FT m_value;
    Point m_position;
};

class CurveSingularity
{
public:
    CurveSingularity();
};


#endif // SINGULARITY_H
