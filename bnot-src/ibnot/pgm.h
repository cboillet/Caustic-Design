#ifndef _PGM_H_
#define _PGM_H_

#include <vector>
#include <string>
#include <fstream>

class PGMImage
{
private:
    unsigned m_w, m_h;
    double m_max_value;
    std::vector<double> m_value;
    double m_max_original;
    std::vector<double> m_original;

public:
    PGMImage()
    {
        m_w = 0;
        m_h = 0;
        m_max_value = 0.0;
        m_max_original = 0.0;
    }
    
    unsigned width () const { return m_w; }
    
    unsigned height() const { return m_h; }
    
    bool isNull() const { return m_value.empty(); }
    
    double max_value() const { return m_max_value; }
    
    unsigned get_index(unsigned i, unsigned j) const { return (j*width() + i); }

    double pixel(unsigned i, unsigned j) const
    {
        unsigned index = get_index(i, j);
        return m_value[index];
    }
    
    bool load(const std::string& filename)
    {
        m_original.clear();
        std::ifstream input(filename.c_str());

        std::string header;
        input >> header;
        input >> m_w >> m_h;        
        input >> m_max_original;
        for (unsigned i = 0; i < m_w*m_h; ++i)
        {
            double value;
            input >> value;
            m_original.push_back(value);
        }
        input.close();

        tonemap(0.5);
        return true;
    }
    
    void save(const std::string& filename) const
    {
        std::ofstream output(filename.c_str());
        output << "P2" << std::endl;
        output << m_w << " " << m_h << std::endl;
        output << "255" << std::endl;
        for (unsigned i = 0; i < m_w*m_h; ++i)
        {
            int value = static_cast<int>(255*m_value[i]);
            output << value << std::endl;
        }
        output.close();
    }
    
    void tonemap(double key)
    {
        std::cout << " tonemapp with key = "<< key << std::endl;        
        m_value = m_original;
        m_max_value = m_max_original;
        
        double maxv = 0.0;
        double maxo = 0.0;
        for (unsigned i = 0; i < m_value.size(); ++i)
        {
            maxv = std::max(maxv, m_value[i]);
            maxo = std::max(maxo, m_original[i]);
        }
        std::cout << "MaxV: " << maxv << " ; " << m_max_value << std::endl;
        std::cout << "MaxO: " << maxo << " ; " << m_max_original << std::endl;
    }
};

#endif
