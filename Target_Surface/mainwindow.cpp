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
}

MainWindow::~MainWindow()
{
    //delete ui;
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

void MainWindow::on_actionLoadModel_triggered()
{
    setModel();
}

//void MainWindow::on_actionCreate_Object_triggered(){
//    bool ok;
//    int height = QInputDialog::getInt(this, tr("height surface (in mm)"), tr("h:"), 1000, 0, 1500000, 1, &ok);
//    int width = QInputDialog::getInt(this, tr("width surface"), tr("w:"), 1000, 0, 1500000, 1, &ok);
//    int depth = QInputDialog::getInt(this, tr("depth surface"), tr("d:"), 1000, 0, 1500000, 1, &ok);
//    //QString filename = QFileDialog::getSaveFileName(this, tr("model"), ".obj");
//    //SurfaceModel(width,height,depth, filename::toStdString());
//}
