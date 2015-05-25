#ifndef _GRID_H_
#define _GRID_H_

template <class EnrichedSegment>
class CGrid
{
private:
    unsigned m_width, m_height;
    std::vector< std::vector<EnrichedSegment> > m_data;
    
public:
    CGrid(unsigned nx, unsigned ny)
    {
        m_width  = nx;
        m_height = ny;
        m_data.resize(nx*ny);
    }

    unsigned get_width() const
    {
        return m_width;
    }
    
    unsigned get_height() const
    {
        return m_height;
    }
    
    unsigned get_index(unsigned i, unsigned j) const
    {
        return i + m_width*j;
    }
        
    bool is_empty(unsigned i, unsigned j) const
    {
        unsigned index = get_index(i, j);
        return m_data[index].empty();
    }
    
    const std::vector<EnrichedSegment>& 
    get_enriched_segments(unsigned i, unsigned j) const
    {
        unsigned index = get_index(i, j);
        return m_data[index];
    }
    
    void append(unsigned i, unsigned j, 
                const EnrichedSegment& enriched_segment)
    {
        unsigned index = get_index(i, j);        
        m_data[index].push_back(enriched_segment);
    }
};

#endif
