// STL
#include <ctime>
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
#include "timer.h"
#include "window.h"
#include "dialog.h"
#include "scene.h"
#include "voronoi_creation.h"
#include "optimal_transport.h"
#include "interpolation.h"
#include "config.h"

int nbpoints; // nb of centroids

MainWindow::MainWindow() : QMainWindow(), Ui_MainWindow(), maxNumRecentFiles(15), recentFileActs(15)
{
	setupUi(this);    
    m_scene = new Scene;
	viewer->set_scene(m_scene);

    source_scene = new Scene;
    compute_scene = new Scene;

    viewer_2->set_scene(source_scene);

    voronoicreator= new VoronoiCreator(m_scene);
    
    m_verbose = 1;
    m_stepX = 0.0;
    m_stepW = 0.0;
    m_epsilon = EPSILON;
    m_frequency = 0;
    m_max_iters = 500;
    m_site_amount = 5000; //camille test interpolation
    //m_site_amount = 1000;
    m_level_max = 3;
    
	// accepts drop events
	setAcceptDrops(true);
	addRecentFiles(menuFile);
    connect(this, SIGNAL(openRecentFile(QString, bool)),
            this, SLOT(open(QString, bool)));


    //open(QString("/home/p/rect2958.png"), false);
    //open(QString("/home/p/Pictures/einstein.png"), false);
    //m_scene->load_image("/home/p/Pictures/einstein.png", false);
    //open(QString("/home/p/Pictures/einstein_2000.dat"), false);

    //open(QString("/home/p/Pictures/white_small.png"), true);
    //open(QString("/home/p/Pictures/white_2000.dat"), true);

}

MainWindow::~MainWindow()
{
    if (m_scene) delete(m_scene);
    if (source_scene) delete(source_scene);
    if (compute_scene) delete(compute_scene);
 
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
            open(filename, false);
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
        emit openRecentFile(action->data().toString(), false);
}

void MainWindow::updateRecentFileActions()
{
	QSettings settings;
	QStringList files = settings.value("recentFileList").toStringList();
    
	int numRecentFiles = qMin(files.size(), int(maxNumberOfRecentFiles()));
	for (int i = 0; i < numRecentFiles; ++i) 
    {
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
    
	for (unsigned int i = 0; i < maxNumberOfRecentFiles(); ++i) 
    {
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

void MainWindow::open(const QString& filename, const bool open_source)
{
    std::cerr << "open ...";
	QApplication::setOverrideCursor(Qt::WaitCursor);
    if (open_source)
    {
        if (is_image(filename)) source_scene->load_image(filename);
        else                    source_scene->load_points(filename);
    }else
    {
        if (is_image(filename)) {
            m_scene->load_image(filename);
            compute_scene->load_image(filename);
        }
        else                    {
            m_scene->load_points(filename);
            compute_scene->load_points(filename);
        }
    }
    QApplication::restoreOverrideCursor();
    std::cerr << "done" << std::endl;

    addToRecentFiles(filename);
	update();
}

void MainWindow::save(const QString& filename, const bool target) const
{
    std::cerr << "save ...";
    QApplication::setOverrideCursor(Qt::WaitCursor);
    if(target)
        m_scene->save_points(filename);
    else
        source_scene->save_points(filename);
    QApplication::restoreOverrideCursor();
    std::cerr << "done" << std::endl;
}

void MainWindow::update()
{
	viewer->repaint();
    viewer_2->repaint();
}

bool MainWindow::is_image(const QString& filename) const
{
    if (filename.contains(".dat", Qt::CaseInsensitive))
        return false;
    return true;
}

void MainWindow::on_actionClear_triggered()
{
    std::cerr << "Clear scene...";
	m_scene->clear();
    std::cerr << "done" << std::endl;
	update();
}

void MainWindow::on_actionLoadWeights_triggered()
{
    QString fileName =
    QFileDialog::getOpenFileName(this, tr("Open Weights"), ".");
    if(fileName.isEmpty()) return;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    std::vector<FT> weights = m_scene->load_weights(fileName);
    m_scene->update_weights(weights, false);
    m_scene->update_triangulation();
    QApplication::restoreOverrideCursor();
    update();
}

void MainWindow::on_actionLoadSingularities_triggered()
{
    QString fileName =
            QFileDialog::getOpenFileName(this, tr("Open Singularities"), ".");

    if(fileName.isEmpty()) return;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    std::vector<PointSingularity> ps;
    std::vector<CurveSingularity> cs;
    m_scene->load_singularities(fileName, ps, cs);
    m_scene->update_singularities(ps, cs);
    QApplication::restoreOverrideCursor();
    update();
}

void MainWindow::on_actionOpenImage_triggered()
{
	QString fileName = 
    QFileDialog::getOpenFileName(this, tr("Open image"), ".");
	if (fileName.isEmpty()) return;
    open(fileName, false);
}

void MainWindow::on_actionLoadSourceImage_triggered()
{
    QString fileName =
    QFileDialog::getOpenFileName(this, tr("Open image"), ".");
    if (fileName.isEmpty()) return;
    open(fileName, true);
}


void MainWindow::on_actionOpenPoints_triggered()
{
	QString fileName = 
    QFileDialog::getOpenFileName(this, tr("Open pointset"), ".dat");
	if (fileName.isEmpty()) return;
    open(fileName, false);
}

void MainWindow::on_actionOpenSourceDAT_triggered()
{
    QString fileName =
            QFileDialog::getOpenFileName(this, tr("Open pointset"), ".dat");
    if (fileName.isEmpty()) return;
    open(fileName, true);
}

void MainWindow::on_actionSavePoints_triggered()
{
    QString filename = 
    QFileDialog::getSaveFileName(this, tr("Save pointset"), ".dat");
	if (filename.isEmpty()) return;
    save(filename);
}

void MainWindow::on_actionSaveSourceDAT_triggered()
{
    QString filename =
    QFileDialog::getSaveFileName(this, tr("Save pointset"), ".dat");
    if (filename.isEmpty()) return;
    save(filename, false);
}

void MainWindow::on_actionSaveEPS_triggered()
{
    QString filename = 
    QFileDialog::getSaveFileName(this, tr("Save pointset"), ".eps");
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

void MainWindow::on_actionResetWeights_triggered()
{
    if (!m_scene->is_valid()) return;
    Timer::start_timer(m_timer, COLOR_GREEN, "reset W");
	
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene->reset_weights();
    source_scene->reset_weights();
	QApplication::restoreOverrideCursor();
    
    Timer::stop_timer(m_timer, COLOR_GREEN);
	update();    
}

void MainWindow::on_actionOptimizePointsLloyd_triggered()
{
    if (!m_scene->is_valid()) return;
    
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene->compute_wcvt_energy();
        n0 = m_scene->count_visible_sites();        
    }
	
    Timer::start_timer(m_timer, COLOR_YELLOW, "opt X via Lloyd");
    std::cout << std::endl;
 
    QApplication::setOverrideCursor(Qt::WaitCursor);
    FT norm = m_scene->optimize_positions_via_lloyd(true);
	QApplication::restoreOverrideCursor();
    
    Timer::stop_timer(m_timer, COLOR_YELLOW);
	update();        

    if (m_verbose > 0)
    {
        n1 = m_scene->count_visible_sites();
        E1 = m_scene->compute_wcvt_energy();
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ")" << std::endl;
        std::cout << "Norm: " << norm << std::endl;
        if (n0 != n1) std::cout << red << "Visible: " << n0 << " -> " << n1 << white << std::endl;
    }    
}

void MainWindow::on_actionOptimizePointsGD_triggered()
{
    if (!m_scene->is_valid()) return;

    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene->compute_wcvt_energy();
        n0 = m_scene->count_visible_sites();        
    }
    
    Timer::start_timer(m_timer, COLOR_YELLOW, "opt X via GD");
    std::cout << std::endl;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    FT norm = m_scene->optimize_positions_via_gradient_ascent(stepX(), true);
	QApplication::restoreOverrideCursor();
    
    Timer::stop_timer(m_timer, COLOR_YELLOW);
	update();        
    
    if (m_verbose > 0)
    {
        n1 = m_scene->count_visible_sites();
        E1 = m_scene->compute_wcvt_energy();
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ")" << std::endl;
        std::cout << "Norm: " << norm << std::endl;
        if (n0 != n1) std::cout << red << "Visible: " << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizeWeightsGD_triggered()
{
    if (!m_scene->is_valid()) return;

    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene->compute_wcvt_energy();
        n0 = m_scene->count_visible_sites();        
    }
    
    Timer::start_timer(m_timer, COLOR_GREEN, "opt W via GD");
    std::cout << std::endl;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene->optimize_weights_via_gradient_descent(stepW(), true);
	QApplication::restoreOverrideCursor();
    
    Timer::stop_timer(m_timer, COLOR_GREEN);
	update();    
    
    if (m_verbose > 0)
    {
        n1 = m_scene->count_visible_sites();
        E1 = m_scene->compute_wcvt_energy();
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ")" << std::endl;
        if (n0 != n1) std::cout << red << "Visible: " << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizeWeightsNewton_triggered()
{
    if (!m_scene->is_valid()) return;
    
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene->compute_wcvt_energy();
        n0 = m_scene->count_visible_sites();        
    }
    
    Timer::start_timer(m_timer, COLOR_GREEN, "opt W via Newton");
    std::cout << std::endl;

	QApplication::setOverrideCursor(Qt::WaitCursor);
    FT norm = m_scene->optimize_weights_via_newton(stepW(), true);
    QApplication::restoreOverrideCursor();
    
    Timer::stop_timer(m_timer, COLOR_GREEN);
	update();        
    
    if (m_verbose > 0)
    {
        n1 = m_scene->count_visible_sites();
        E1 = m_scene->compute_wcvt_energy();
        std::cout << "Norm: " << norm << std::endl;
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ")" << std::endl;
        if (n0 != n1) std::cout << red << "Visible: " << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizeWeightsGDUntil_triggered()
{
    if (!m_scene->is_valid()) return;
    
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene->compute_wcvt_energy();
        n0 = m_scene->count_visible_sites();        
    }
    
    Timer::start_timer(m_timer, COLOR_GREEN, "opt W via GD until");
    std::cout << std::endl;

	QApplication::setOverrideCursor(Qt::WaitCursor);
    FT threshold = m_scene->compute_weight_threshold(epsilon());    
    unsigned iters = m_scene->optimize_weights_via_gradient_descent_until_converge(stepW(), 
                                                                                  threshold, 
                                                                                  frequency(),
                                                                                  max_iters());
    QApplication::restoreOverrideCursor();
    
    Timer::stop_timer(m_timer, COLOR_GREEN);
	update();
    
    if (m_verbose > 0)
    {
        n1 = m_scene->count_visible_sites();
        E1 = m_scene->compute_wcvt_energy();
        std::cout << "Iters: " << iters << " (" << max_iters() << ")" << std::endl;
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ")" << std::endl;
        if (n0 != n1) std::cout << red << "Visible: " << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionOptimizeWeightsNewtonUntil_triggered()
{
    if (!m_scene->is_valid()) return;
    
    FT E0, E1;
    unsigned n0, n1;
    if (m_verbose > 0)
    {
        E0 = m_scene->compute_wcvt_energy();
        n0 = m_scene->count_visible_sites();        
    }
    
    Timer::start_timer(m_timer, COLOR_GREEN, "opt W via Newton until");
    std::cout << std::endl;

	QApplication::setOverrideCursor(Qt::WaitCursor);
    FT threshold = m_scene->compute_weight_threshold(epsilon());
    unsigned iters = m_scene->optimize_weights_via_newton_until_converge(stepW(), 
                                                                         threshold, 
                                                                         frequency(),
                                                                         max_iters());
    QApplication::restoreOverrideCursor();
    
    Timer::stop_timer(m_timer, COLOR_GREEN);
	update();
    
    if (m_verbose > 0)
    {
        n1 = m_scene->count_visible_sites();
        E1 = m_scene->compute_wcvt_energy();
        std::cout << "Iters: " << iters << " (" << max_iters() << ")" << std::endl;
        std::cout << "DE: " << E1-E0 << " (" << (E0 > E1) << ")" << std::endl;
        if (n0 != n1) std::cout << red << "Visible: " << n0 << " -> " << n1 << white << std::endl;
    }
}

void MainWindow::on_actionFullOptimization_triggered()
{
    if (!m_scene->is_valid()) return;
    Timer::start_timer(m_timer, COLOR_RED, "full opt");
    std::cout << std::endl;
	
    QApplication::setOverrideCursor(Qt::WaitCursor);
    unsigned iters = m_scene->optimize_all(stepW(), stepX(), 500,
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
    
    int max_teleport = QInputDialog::getInt(this, tr("Teleport"), tr("max teleport:"),
                                                100, 0, 10000, 1, &ok);

    if (!ok) return;
    
    Timer::start_timer(m_timer, COLOR_BLUE, "Breaking Regularity");
    std::cout << std::endl;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_scene->detect_and_break_regularity(radius, max_teleport);
    QApplication::restoreOverrideCursor();
    
    Timer::stop_timer(m_timer, COLOR_BLUE);
    update();
}

//////////
// DATA //
//////////

void MainWindow::on_actionToggleInvert_toggled()
{
    m_scene->toggle_invert();
    update();
}

void MainWindow::on_actionGenerateVariablePoints_triggered()
{
	bool ok;
        int nb = QInputDialog::getInt(this, tr("Nb Sites"), tr("n:"), 1024, 3, 1500000, 1, &ok);
    if (!ok) return;
    
	QApplication::setOverrideCursor(Qt::WaitCursor);
    Timer::start_timer(m_timer, COLOR_WHITE, "sample domain adapted to image");
    m_scene->generate_random_sites_based_on_image(unsigned(nb));
    Timer::stop_timer(m_timer, COLOR_WHITE);
    QApplication::restoreOverrideCursor();
    update();

    std::cout << "Insert " << m_scene->count_visible_sites() << " sites adapted to image" << std::endl;
}

//////////
// VIEW //
//////////

void MainWindow::on_actionViewImageGrid_toggled()
{
    viewer->toggle_view_image_grid();
    viewer_2->toggle_view_image_grid();
    update();
}

void MainWindow::on_actionViewImage_toggled()
{
    viewer->toggle_view_image();
    viewer_2->toggle_view_image();
    update();
}

void MainWindow::on_actionViewDomain_toggled()
{
    viewer->toggle_view_domain();
    viewer_2->toggle_view_domain();
    update();
}

void MainWindow::on_actionViewPoints_toggled()
{
    viewer->toggle_view_points();
    viewer_2->toggle_view_points();
    update();
}

void MainWindow::on_actionViewVertices_toggled()
{
    viewer->toggle_view_vertices();
    viewer_2->toggle_view_vertices();
    update();
}

void MainWindow::on_actionViewEdges_toggled()
{
    viewer->toggle_view_edges();
    viewer_2->toggle_view_edges();
    update();
}

void MainWindow::on_actionViewFaces_toggled()
{
    viewer->toggle_view_faces();
    viewer_2->toggle_view_faces();
    update();
}

void MainWindow::on_actionViewWeights_toggled()
{
    viewer->toggle_view_weights();
    viewer_2->toggle_view_weights();
    update();
}

void MainWindow::on_actionViewDual_toggled()
{
    viewer->toggle_view_dual();
    viewer_2->toggle_view_dual();
    update();
}

void MainWindow::on_actionViewBoundedDual_toggled()
{
    viewer->toggle_view_bounded_dual();
    viewer_2->toggle_view_bounded_dual();
    update();
}

void MainWindow::on_actionViewPixels_toggled()
{
    viewer->toggle_view_pixels();
    viewer_2->toggle_view_pixels();
    update();
}

void MainWindow::on_actionViewCapacity_toggled()
{
    viewer->toggle_view_capacity();
    viewer_2->toggle_view_capacity();
    update();
}

void MainWindow::on_actionViewRegularity_toggled()
{
    viewer->toggle_view_regularity();
    viewer_2->toggle_view_regularity();
    update();
}

void MainWindow::on_actionViewRegularSites_toggled()
{
    viewer->toggle_view_regular_sites();
    viewer_2->toggle_view_regular_sites();
    update();
}

void MainWindow::on_actionViewVariance_toggled()
{
    viewer->toggle_view_variance();
    viewer_2->toggle_view_variance();
    update();
}

void MainWindow::on_actionViewBarycenter_toggled()
{
    viewer->toggle_view_barycenter();
    viewer_2->toggle_view_barycenter();
    update();
}

void MainWindow::on_actionDrawMovement_toggled()
{
    viewer->toggle_view_movement();
    viewer_2->toggle_view_movement();
    update();
}

void MainWindow::on_actionViewWeightHistogram_toggled()
{
    viewer->toggle_view_weight_histogram();
    viewer_2->toggle_view_weight_histogram();
    update();
}

void MainWindow::on_actionViewCapacityHistogram_toggled()
{
    viewer->toggle_view_capacity_histogram();
    viewer_2->toggle_view_capacity_histogram();
    update();
}

void MainWindow::on_actionViewGradient_toggled()
{
    viewer->toggle_view_gradient();
    viewer_2->toggle_view_gradient();
    update();
}

void MainWindow::on_actionViewPopups_toggled()
{
    viewer->toggle_view_newly_visible();
    viewer_2->toggle_view_newly_visible();
    update();
}

void MainWindow::on_actionSetParameters_triggered()
{
    Dialog dlg;
    dlg.set_all_ranges();
    
    dlg.set_verbose(m_verbose);
    dlg.set_stepX(m_stepX);
    dlg.set_stepW(m_stepW);
    dlg.set_epsilon(m_epsilon);
    dlg.set_frequency(m_frequency);
    dlg.set_max_iters(m_max_iters);
    dlg.set_point_size(viewer->point_size());
    dlg.set_vertex_size(viewer->vertex_size());
    dlg.set_line_thickness(viewer->line_thickness());
    dlg.set_hist_range(viewer->histogram_range());
    dlg.set_hist_nbins(viewer->histogram_nbins());
    dlg.set_tau(m_scene->get_tau());
    dlg.set_site_amount(m_site_amount);
    dlg.set_level_max(m_level_max);
    
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
        m_scene->set_tau(dlg.get_tau());
        m_level_max = dlg.get_level_max();
        m_site_amount = dlg.get_site_amount();
        update();
    }
}

void MainWindow::on_actionToggleTimer_toggled()
{
    m_scene->toggle_timer();
}

void MainWindow::on_actionToggleFixedConnectivity_toggled()
{
    m_scene->toggle_connectivity();
}

void MainWindow::on_actionCountSitesPerBin_triggered()
{
    bool ok;
    int nb = QInputDialog::getInt(this, tr("NbBins"), tr("bins:"), 6, 1, 100, 1, &ok);
    if (!ok) return;
    m_scene->count_sites_per_bin(nb);
}

void MainWindow::on_actionVoronoiCreation_triggered(){
    bool ok;
    nbpoints = QInputDialog::getInt(this, tr("NBpoint"), tr("Number of voronoi Centroids:"), 20, 20, 200000, 100, &ok);
    if (!ok) return;


    /*int nbiter = QInputDialog::getInt(this, tr("NBllyod"), tr("Number of Llyod simplification:"), 5, 1, 40, 1, &ok);
    if (!ok) return;
    voronoicreator->init_points(nbpoints,m_scene);
    voronoicreator->init_points(nbpoints,compute_scene);
    for (uint i=0; i<nbiter; i++){
        std::cout << "(" << (i+1) << "/" << nbiter << "): ";
        voronoicreator->apply_lloyd_optimization(m_scene);
        voronoicreator->apply_lloyd_optimization(compute_scene);
    }*/

    if(voronoicreator->generate_voronoi(m_scene, nbpoints, epsilon(), viewer))
        this->on_actionSavePoints_triggered();
}

void MainWindow::on_actionComputeInterpolation_triggered(){
    QString image = QFileDialog::getOpenFileName(this, tr("Open image"), ".");
    if (image.isEmpty()) return;
    open(image, false);
    open(image, true);
    QString dats = QFileDialog::getOpenFileName(this, tr("Open pointset"), ".dat");
    if (dats.isEmpty()) return;
    source_scene->load_points(dats);
    m_scene->load_points(dats);
    compute_scene->load_points(dats);
    on_actionLoadWeights_triggered();



     QString fileName = QFileDialog::getOpenFileName(this, tr("Load light origin points"), ".");
    if (fileName.isEmpty()) return;
    else m_scene->load_interpolation_points(fileName, m_scene->getLightPointsSource());



    std::cout << "onActionComputeInterpolation" << std::endl;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Interpolation inter = Interpolation(source_scene, m_scene, compute_scene, m_site_amount,this);
    inter.runInterpolation();
    QApplication::restoreOverrideCursor();
    //viewer_2->toggle_view_Xr();
    //update();
}

void MainWindow::on_actionCalculateOptimalTransport_triggered()
{
    std::cout << "onActionComputeOptimalTransport" << std::endl;

    QMessageBox msgBox;
    msgBox.setText("Running Optimal Transport Calculation.");
    msgBox.setInformativeText("Do you want to run the gradient descent approach?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
    msgBox.setDefaultButton(QMessageBox::Yes);
    int ret = msgBox.exec();

    bool gradient_descent;

    switch (ret) {
      case QMessageBox::Yes:
        gradient_descent = true;
        break;
      case QMessageBox::No:
        gradient_descent = false;
        break;
      case QMessageBox::Cancel:
        return;
    }

    Timer::start_timer(m_timer, COLOR_BLUE, "OTM");

    QApplication::setOverrideCursor(Qt::WaitCursor);
    OptimalTransport ot = OptimalTransport(m_scene, source_scene, this, viewer_2, m_level_max, m_site_amount);
    ot.runOptimalTransport(gradient_descent);

    Timer::stop_timer(m_timer, COLOR_BLUE);

    QString filename =
    QFileDialog::getSaveFileName(this, tr("Save weights"), ".weight");
    if (!filename.isEmpty())
        source_scene->save_weights(filename);
    QApplication::restoreOverrideCursor();
    update();

}



