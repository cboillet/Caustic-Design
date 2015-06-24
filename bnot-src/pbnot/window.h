#ifndef _WINDOW_
#define _WINDOW_

// Qt
#include <QWidget>
#include <QString>

// local
#include "scene.h"
#include "ui_pbnot.h"

class MainWindow : public QMainWindow, public Ui_MainWindow
{
	Q_OBJECT
    
private:
    Scene m_scene;
    
    double m_stepX;
    double m_stepW;
    double m_epsilon;
    unsigned m_verbose;
    unsigned m_frequency;
    unsigned m_max_iters;
    std::vector<double> m_timer;

    unsigned maxNumRecentFiles;
    QAction* recentFilesSeparator;
    QVector<QAction*> recentFileActs;
    
public:
	MainWindow();
	~MainWindow();
    
    double& stepX() { return m_stepX; }
    const double stepX() const { return m_stepX; }

    double& stepW() { return m_stepW; }
    const double stepW() const { return m_stepW; }

    double& epsilon() { return m_epsilon; }
    const double epsilon() const { return m_epsilon; }

    unsigned& frequency() { return m_frequency; }
    const unsigned frequency() const { return m_frequency; } 

    unsigned& max_iters() { return m_max_iters; }
    const unsigned max_iters() const { return m_max_iters; } 
    
    const unsigned verbose() const { return m_verbose; }
    unsigned& verbose() { return m_verbose; }
    
protected slots:
    // drag and drop
    void dropEvent(QDropEvent *event);
    void closeEvent(QCloseEvent *event);
    void dragEnterEvent(QDragEnterEvent *event);
    
    // recent files
    void openRecentFile_aux();
    void updateRecentFileActions();
    void addToRecentFiles(QString fileName);    
	void addRecentFiles(QMenu* menu, QAction* insertBefore = 0);    
	unsigned int maxNumberOfRecentFiles() const {return maxNumRecentFiles;}

    // io
    void update();
    void open(const QString& file);
    void save(const QString& filename) const;

    // file
    void on_actionClear_triggered();
    void on_actionSaveEPS_triggered();
    void on_actionSnapshot_triggered();
    void on_actionOpenPoints_triggered();
    void on_actionSavePoints_triggered();
        
    // view
    void on_actionViewDomain_toggled();
    void on_actionViewPoints_toggled();    
    void on_actionViewVertices_toggled();
    void on_actionViewEdges_toggled();
    void on_actionViewWeights_toggled();
    void on_actionViewDual_toggled();
    void on_actionViewCapacity_toggled();
    void on_actionViewVariance_toggled();
    void on_actionViewBarycenter_toggled();
    void on_actionViewRegularity_toggled();
    void on_actionViewRegularSites_toggled();
    void on_actionViewWeightHistogram_toggled();
    void on_actionViewCapacityHistogram_toggled();

    // data
    void on_actionGenerateGrid_triggered();
    void on_actionGenerateHextille_triggered();
    void on_actionGenerateRandomPoints_triggered();
    
    // algorithm
    void on_actionToggleTimer_toggled();
    void on_actionPrintAlpha_triggered();
    void on_actionSetParameters_triggered();
    void on_actionToggleFixedConnectivity_toggled();
        
    void on_actionResetWeights_triggered();
    void on_actionFullOptimization_triggered();
    void on_actionBreak_Regularity_triggered();
    
    void on_actionOptimizeWeightsGD_triggered();
    void on_actionOptimizeWeightsNewton_triggered();
    void on_actionOptimizeWeightsGDUntil_triggered();
    void on_actionOptimizeWeightsNewtonUntil_triggered();
    
    void on_actionOptimizePointsGD_triggered();
    void on_actionOptimizePointsLloyd_triggered();
    
signals:
    void openRecentFile(QString filename);
};

#endif // _WINDOW_
