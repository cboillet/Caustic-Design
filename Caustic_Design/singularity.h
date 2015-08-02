#ifndef SINGULARITY_H
#define SINGULARITY_H

template <class Kernel>
class PSingularity
{
public:
    typedef typename Kernel::FT       FT;
    typedef typename Kernel::Point_2  Point;

private:
    FT m_value;
    Point m_position;

public:
    PSingularity():
        m_value(0.0)
    {}

    PSingularity(const FT value):
        m_value(value)
    {}

    PSingularity(const Point &position, const FT value):
        m_position(position),
        m_value(value)
    {}

    void set_value(const FT value) { m_value = value; }

    const FT get_value(){ return m_value; }
    const Point& get_position(){ return m_position; }

};

class CurveSingularity
{
public:
    CurveSingularity(){}
};


#endif // SINGULARITY_H
