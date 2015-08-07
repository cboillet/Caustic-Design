#include "global.h"
// Qt
#include <QtGui>
#include <QDialog>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QClipboard>

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "SurfaceModel.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), Ui_MainWindow()//,render(new Renderer(5, this, "Target Surface"))
{
    setupUi(this);
    //setModel();
    //viewer = new Renderer(5, this, "Target Surface");
    optimizer = new TargetOptimization();
}

MainWindow::~MainWindow()
{
    //delete ui;
    delete optimizer;
}

void MainWindow::setModel()
{
    QString modelToLoad = QFileDialog::getOpenFileName(this, tr("Open 3D model to load"), ".");
    if(modelToLoad.isEmpty()) return;
    viewer->model.loadModel(modelToLoad.toStdString());
    update();
}

void MainWindow::on_actionSaveModel_triggered()
{
    std::cout << "save model" << std::endl;
    QString saveModelName = QFileDialog::getSaveFileName(this, tr("Save 3D model"));
    if(saveModelName.isEmpty()) return;
    viewer->model.exportModel(saveModelName.toUtf8().constData());
}

void MainWindow::on_actionSave_Vertices_triggered()
{
    std::cout << "save vertices" << std::endl;
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save vertices"));
    if(fileName.isEmpty()) return;
    viewer->model.meshes[0].exportVertices(fileName);
}

void MainWindow::on_actionLoadModel_triggered()
{
    setModel();
}

void MainWindow::on_actionGenerateTriangles_triggered()
{
    viewer->repaint();
    //update();
}

void MainWindow::on_actionRunTargetOptimization_triggered()
{
    optimizer->runCeresTest();
}

void MainWindow::on_actionExit_triggered()
{
    QApplication::quit();
}
