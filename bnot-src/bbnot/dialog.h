#ifndef _DIAL_OPT_H_
#define _DIAL_OPT_H_

#include <cmath>
#include "ui_dialog.h"

class Dialog : public QDialog, private Ui::Dialog
{
    Q_OBJECT

public:
    Dialog(QWidget *parent = 0)
    {
        setupUi(this);
    }
    
    void set_all_ranges()
    {
        tau_spin->setRange(0.0, 1000.0);
        verbose_spinbox->setRange(0, 1);
        stepX_spin->setRange(-1.0, 1000.);
        stepW_spin->setRange(-1.0, 1000.);
        epsilon_spin->setRange(0.0, 10.0);
        frequency_spin->setRange(0, 10000);
        max_iters_spin->setRange(1, 10000);
        hist_nbins_spin->setRange(2, 4096);
        hist_range_spin->setRange(0.001, 1.0);
    }
    
    double get_line_thickness() const { return thickness_spinbox->value(); }
    void set_line_thickness(const double t) { thickness_spinbox->setValue(t); }
    
    double get_point_size() const { return point_size_spinbox->value(); }
    void set_point_size(const double t) { point_size_spinbox->setValue(t); }
    
    double get_vertex_size() const { return vertex_size_spinbox->value(); }
    void set_vertex_size(const double t) { vertex_size_spinbox->setValue(t); }    
    
    int get_verbose() const { return verbose_spinbox->value(); }
    void set_verbose(const int verbose) { verbose_spinbox->setValue(verbose); }    
    
    double get_stepX() const { return stepX_spin->value(); }
    void set_stepX(const double step) { stepX_spin->setValue(step); }

    double get_stepW() const { return stepW_spin->value(); }
    void set_stepW(const double step) { stepW_spin->setValue(step); }

    double get_epsilon() const { return epsilon_spin->value(); }
    void set_epsilon(const double epsilon) { epsilon_spin->setValue(epsilon); }

    int get_frequency() const { return frequency_spin->value(); }
    void set_frequency(const int freq) { frequency_spin->setValue(freq); }

    int get_max_iters() const { return max_iters_spin->value(); }
    void set_max_iters(const int max_iters) { max_iters_spin->setValue(max_iters); }

    int get_hist_nbins() const { return hist_nbins_spin->value(); }
    void set_hist_nbins(const int nbins) { hist_nbins_spin->setValue(nbins); }
    
    double get_hist_range() const { return hist_range_spin->value(); }
    void set_hist_range(const double range) { hist_range_spin->setValue(range); }
    
    double get_tau() const { return tau_spin->value(); }
    void set_tau(const double tau) { tau_spin->setValue(tau); }
};

#endif
