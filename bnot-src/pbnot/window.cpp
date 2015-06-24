// STL
#include <fstream>
#include <sstream>

// Qt
#include <QtGui>
#include <QDialog>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QClipboard>

// local
#include "window.h"
#include "dialog.h"
#include "timer.h"

MainWindow::MainWindow() :
QMainWindow(), Ui_MainWindow(),
maxNumRecentFiles(15), recentFileActs(15)
{
    setupUi(this);
    viewer->set_scene(&m_scene);

    m_verbose = 0;
    m_stepX = 0.0;
    m_stepW = 0.0;
    m_epsilon = 1.0;
    m_frequency = 0;
    m_max_iters = 500;
    
    // accepts drop events
    setAcceptDrops(true);
    addRecentFiles(menuFile);
    connect(this, SIGNAL(openRecentFile(QString)), 
            this, SLOT(open(QString)));
}

MainWindow::~MainWindow()
{
}

void MainWindow::addToRecentFiles(QString fileName)
{
    QSettings settings;
    QStringList files = settings.value("recentFileList").toStringList();
    files.removeAll(fileName);
    files.prepend(fileName);
    while (files.size() > int(maxNumRecentFiles))
        files.removeLast();
    settings.setValue("recentFileList", files);
    updateRecentFileActions();
}

void MainWindow::dropEvent(QDropEvent *event)
{
    Q_FOREACH(QUrl url, event->mimeData()->urls())
    {
        QString filename = url.toLocalFile();
        if (!filename.isEmpty())
        {
            QTextStream(stderr) << QString("dropEvent(\"%1\")\n").arg(filename);
            open(filename);
        }
    }
    event->acceptProposedAction();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    event->accept();
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
    if (event->mimeData()->hasFormat("text/uri-list"))
        event->acceptProposedAction();
}

void MainWindow::openRecentFile_aux()
{
    QAction* action = qobject_cast<QAction*>(sender());
    if (action)
        emit openRecentFile(action->data().toString());
}

void MainWindow::updateRecentFileActions()
{
    QSettings settings;
    QStringList files = settings.value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), int(maxNumberOfRecentFiles()));
    for (int i = 0; i < numRecentFiles; ++i) {
        QString strippedName = QFileInfo(files[i]).fileName();
        QString text = tr("&%1 %2").arg(i).arg(strippedName);
        recentFileActs[i]->setText(text);
        recentFileActs[i]->setData(files[i]);
        recentFileActs[i]->setVisible(true);
    }
    for (unsigned j = numRecentFiles; j < maxNumberOfRecentFiles(); ++j)
        recentFileActs[j]->setVisible(false);

    recentFilesSeparator->setVisible(numRecentFiles > 0);
}

void MainWindow::addRecentFiles(QMenu* menu, QAction* insertBeforeAction)
{
    if (insertBeforeAction)
        recentFilesSeparator = menu->insertSeparator(insertBeforeAction);
    else
        recentFilesSeparator = menu->addSeparator();
    recentFilesSeparator->setVisible(false);

    for (unsigned int i = 0; i < maxNumberOfRecentFiles(); ++i) {
        recentFileActs[i] = new QAction(this);
        recentFileActs[i]->setVisible(false);
        connect(recentFileActs[i], SIGNAL(triggered()), this, SLOT(openRecentFile_aux()));
        if (insertBeforeAction)
            menu->insertAction(insertBeforeAction, recentFileActs[i]);
        else
            menu->addAction(recentFileActs[i]);
    }
    updateRecentFileActions();
}

void MainWindow::open(const QString& filename)
{
    std::cerr << "read ...";
    QApplication::setOverrideCursor(Qt::WaitCursor);
    unsigned nb = m_scene.load(filename);
    QApplication::restoreOverrideCursor();
    addToRecentFiles(filename);
    update();
    std::cerr << "done (" << nb << " points)" << std::endl;
}

void MainWindow::save(const QString& filename) const
{
    std::cerr << "save ...";
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene.save(filename);
    QApplication::restoreOverrideCursor();
    std::cerr << "done" << std::endl;
}

void MainWindow::update()
{
    viewer->repaint();
}

void MainWindow::on_actionClear_triggered()
{
    std::cerr << "Clear scene...";
    m_scene.clear();
    update();
    std::cerr << "done" << std::endl;
}

void MainWindow::on_actionOpenPoints_triggered()
{
    QString fileName = 
    QFileDialog::getOpenFileName(this, tr("Open pointset"), ".dat");
    if (fileName.isEmpty()) return;
    open(fileName);
}

void MainWindow::on_actionSavePoints_triggered()
{
    QString filename = 
    QFileDialog::getSaveFileName(this, tr("Save pointset"), ".dat");
    if (filename.isEmpty()) return;
    save(filename);
}

void MainWindow::on_actionSaveEPS_triggered()
{
    QString filename = 
    QFileDialog::getSaveFileName(this, tr("Save eps"), ".eps");
    if (filename.isEmpty()) return;
    save(filename);
}

void MainWindow::on_actionSnapshot_triggered()
{
    std::cout << "snapshot...";
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QClipboard *qb = QApplication::clipboard();
    viewer->makeCurrent();
    viewer->raise();
    QImage snapshot = viewer->grabFrameBuffer(true);
    qb->setImage(snapshot);
    QApplication::restoreOverrideCursor();
    std::cout << "done" << std::endl;
}

///////////////
// ALGORITHM //
///////////////

void MainWindow::on_actionToggleTimer_toggled()
{
    m_scene.toggle_timer();
}

void MainWindow::on_actionPrintAlpha_triggered()
{
    m_scene.compute_alpha();
}

void MainWindow::on_actionToggleFixedConnectivity_toggled()
{
    m_scene.toggle_connectivity();
}

void MainWindow::on_actionSetParameters_triggered()
{
    Dialog dlg;
    dlg.set_all_ranges();
    
    dlg.set_stepX(m_stepX);
    dlg.set_stepW(m_stepW);
    dlg.set_epsilon(m_epsilon);
    dlg.set_verbose(m_verbose);
    dlg.set_frequency(m_frequency);
    dlg.set_max_iters(m_max_iters);
    dlg.set_point_size(viewer->point_size());
    dlg.set_vertex_size(viewer->vertex_size());
    dlg.set_line_thickness(viewer->line_thickness());
    dlg.set_hist_range(viewer->histogram_range());
    dlg.set_hist_nbins(viewer->histogram_nbins());
    dlg.set_tau(m_scene.get_tau());
    
    if (dlg.exec() == QDialog::Accepted)
    {
        m_verbose = dlg.get_verbose();
        m_stepX = dlg.get_stepX();
        m_stepW = dlg.get_stepW();
        m_epsilon = dlg.get_epsilon();
        m_frequency = dlg.get_frequency();
        m_max_iters = dlg.get_max_iters();
        viewer->point_size()  = dlg.get_point_size();
        viewer->vertex_size() = dlg.get_vertex_size();
        viewer->line_thickness() = dlg.get_line_thickness();
        viewer->histogram_range() = dlg.get_hist_range();
        viewer->histogram_nbins() = dlg.get_hist_nbins();
        m_scene.set_tau(dlg.get_tau());
    }
    
    update();
}


void MainWindow::on_actionResetWeights_triggered()
{
    Timer::start_timer(m_timer, COLOR_GREEN, "reset W");

    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene.reset_weights();
    QApplication::restoreOverrideCursor();

    Timer::stop_timer(m_timer, COLOR_GREEN);
    update();
}

void MainWindow::on_actionOptimizePointsLloyd_triggered()
{
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene.compute_wcvt_energy();
        n0 = m_scene.count_visible_sites();
    }

    Timer::start_timer(m_timer, COLOR_YELLOW, "opt X via Lloyd");
    std::cout << std::endl;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    FT norm = m_scene.optimize_positions_via_lloyd(true);
    QApplication::restoreOverrideCursor();

    Timer::stop_timer(m_timer, COLOR_YELLOW);
    update();

    if (m_verbose > 0)
    {
        n1 = m_scene.count_visible_sites();
        E1 = m_scene.compute_wcvt_energy();
        std::cout << "XNorm: " << norm << std::endl;
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ") " << std::endl;
        if (n0 != n1) std::cout << "Visible: " << red << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizePointsGD_triggered()
{
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene.compute_wcvt_energy();
        n0 = m_scene.count_visible_sites();
    }

    Timer::start_timer(m_timer, COLOR_YELLOW, "opt X via GD");
    std::cout << std::endl;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    FT norm = m_scene.optimize_positions_via_gradient_ascent(stepX(), true);
    QApplication::restoreOverrideCursor();

    Timer::stop_timer(m_timer, COLOR_YELLOW);
    update();

    if (m_verbose > 0)
    {
        n1 = m_scene.count_visible_sites();
        E1 = m_scene.compute_wcvt_energy();
        std::cout << "XNorm: " << norm << std::endl;
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ") " << E0 << " -> " << E1 << std::endl;
        if (n0 != n1) std::cout << "Visible: " << red << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizeWeightsGD_triggered()
{
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene.compute_wcvt_energy();
        n0 = m_scene.count_visible_sites();
    }

    Timer::start_timer(m_timer, COLOR_GREEN, "opt W via GD");
    std::cout << std::endl;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene.optimize_weights_via_gradient_descent(stepW(), true);
    QApplication::restoreOverrideCursor();

    Timer::stop_timer(m_timer, COLOR_GREEN);
    update();

    if (m_verbose > 0)
    {
        n1 = m_scene.count_visible_sites();
        E1 = m_scene.compute_wcvt_energy();
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ")" << std::endl;
        if (n0 != n1) std::cout << "Visible: " << red << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizeWeightsNewton_triggered()
{
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene.compute_wcvt_energy();
        n0 = m_scene.count_visible_sites();
    }

    Timer::start_timer(m_timer, COLOR_GREEN, "opt W via Newton");
    std::cout << std::endl;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    FT norm = m_scene.optimize_weights_via_newton(stepW(), true);
    QApplication::restoreOverrideCursor();

    Timer::stop_timer(m_timer, COLOR_GREEN);
    update();

    if (m_verbose > 0)
    {
        n1 = m_scene.count_visible_sites();
        E1 = m_scene.compute_wcvt_energy();
        std::cout << "WGrad: " << norm << std::endl;
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ") " << E0 << " -> " << E1 << std::endl;
        if (n0 != n1) std::cout << "Visible: " << red << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizeWeightsGDUntil_triggered()
{
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene.compute_wcvt_energy();
        n0 = m_scene.count_visible_sites();
    }

    Timer::start_timer(m_timer, COLOR_GREEN, "opt W via GD until");
    std::cout << std::endl;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    FT threshold = m_scene.compute_weight_threshold(epsilon());    
    unsigned iters = m_scene.optimize_weights_via_gradient_descent_until_converge(stepW(),
                                                                                  threshold,
                                                                                  frequency(),
                                                                                  max_iters());
    QApplication::restoreOverrideCursor();

    Timer::stop_timer(m_timer, COLOR_GREEN);
    update();

    if (m_verbose > 0)
    {
        n1 = m_scene.count_visible_sites();
        E1 = m_scene.compute_wcvt_energy();
        std::cout << "Iters: " << iters << " (" << max_iters() << ")" << std::endl;
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ") " << E0 << " -> " << E1 << std::endl;
        if (n0 != n1) std::cout << "Visible: " << red << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizeWeightsNewtonUntil_triggered()
{
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene.compute_wcvt_energy();
        n0 = m_scene.count_visible_sites();
    }

    Timer::start_timer(m_timer, COLOR_GREEN, "opt W via Newton until");
    std::cout << std::endl;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    FT threshold = m_scene.compute_weight_threshold(epsilon());    
    unsigned iters = m_scene.optimize_weights_via_newton_until_converge(stepW(),
                                                                        threshold,
                                                                        frequency(),
                                                                        max_iters());
    QApplication::restoreOverrideCursor();

    Timer::stop_timer(m_timer, COLOR_GREEN);
    update();

    if (m_verbose > 0)
    {
        n1 = m_scene.count_visible_sites();
        E1 = m_scene.compute_wcvt_energy();
        std::cout << "Iters: " << iters << " (" << max_iters() << ")" << std::endl;
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ") " << E0 << " -> " << E1 << std::endl;
        if (n0 != n1) std::cout << "Visible: " << red << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionFullOptimization_triggered()
{
    Timer::start_timer(m_timer, COLOR_RED, "full opt"); 
    std::cout << std::endl;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    unsigned iters = m_scene.optimize_all(stepW(), stepX(), 500,
                                          epsilon(), max_iters(), 
                                          std::cout);
    QApplication::restoreOverrideCursor();
    
    std::cout << "Iters: " << iters << " (" << max_iters() << ")" << std::endl;
    Timer::stop_timer(m_timer, COLOR_RED);
    update();
}

void MainWindow::on_actionBreak_Regularity_triggered()
{
    bool ok;
    double radius = QInputDialog::getDouble(this, tr("Radius"), tr("radius:"), 
                                            1.0, 0.0, 1.0, 5, &ok);
    if (!ok) return;
    
    int max_teleport = QInputDialog::getInteger(this, tr("Teleport"), tr("max teleport:"), 
                                                100, 0, 10000, 1, &ok);
    if (!ok) return;

    Timer::start_timer(m_timer, COLOR_BLUE, "Breaking Regularity");
    std::cout << std::endl;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene.detect_and_break_regularity(radius, max_teleport);
    QApplication::restoreOverrideCursor();

    Timer::stop_timer(m_timer, COLOR_BLUE);
    update();
}

//////////
// DATA //
//////////

void MainWindow::on_actionGenerateRandomPoints_triggered()
{
    bool ok;
    int nb = QInputDialog::getInteger(this, tr("Nb Sites"), tr("n:"), 1024, 3, 1000000, 1, &ok);
    if (!ok) return;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene.generate_random_sites(unsigned(nb));
    QApplication::restoreOverrideCursor();
    update();
}

void MainWindow::on_actionGenerateGrid_triggered()
{
    bool ok;
    int nx = QInputDialog::getInteger(this, tr("Nx"), tr("nx:"), 32, 3, 10000, 1, &ok);
    if (!ok) return;
    
    int ny = QInputDialog::getInteger(this, tr("Ny"), tr("ny:"), 32, 3, 10000, 1, &ok);
    if (!ok) return;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene.generate_regular_grid(unsigned(nx), unsigned(ny));
    QApplication::restoreOverrideCursor();
    update();
}

void MainWindow::on_actionGenerateHextille_triggered()
{
    bool ok;
    int nb = QInputDialog::getInteger(this, tr("Nb"), tr("nb:"), 32, 1, 10000, 1, &ok);
    if (!ok) return;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene.generate_hextille_grid(unsigned(nb));
    QApplication::restoreOverrideCursor();
    update();
}

//////////
// VIEW //
//////////

void MainWindow::on_actionViewDomain_toggled()
{
    viewer->toggle_view_domain();
    update();
}

void MainWindow::on_actionViewPoints_toggled()
{
    viewer->toggle_view_points();
    update();
}

void MainWindow::on_actionViewVertices_toggled()
{
    viewer->toggle_view_vertices();
    update();
}

void MainWindow::on_actionViewEdges_toggled()
{
    viewer->toggle_view_edges();
    update();
}

void MainWindow::on_actionViewWeights_toggled()
{
    viewer->toggle_view_weights();
    update();
}

void MainWindow::on_actionViewDual_toggled()
{
    viewer->toggle_view_dual();
    update();
}

void MainWindow::on_actionViewCapacity_toggled()
{
    viewer->toggle_view_capacity();
    update();
}

void MainWindow::on_actionViewVariance_toggled()
{
    viewer->toggle_view_variance();
    update();
}

void MainWindow::on_actionViewBarycenter_toggled()
{
    viewer->toggle_view_barycenter();
    update();
}

void MainWindow::on_actionViewWeightHistogram_toggled()
{
    viewer->toggle_view_weight_histogram();
    update();
}

void MainWindow::on_actionViewCapacityHistogram_toggled()
{
    viewer->toggle_view_capacity_histogram();
    update();
}

void MainWindow::on_actionViewRegularity_toggled()
{
    viewer->toggle_view_regularity();
    update();
}

void MainWindow::on_actionViewRegularSites_toggled()
{
    viewer->toggle_view_regular_sites();
    update();
}
