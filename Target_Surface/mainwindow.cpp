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
    QMainWindow(parent),
    ui(new Ui::MainWindow),render(new Renderer(5, this, "Target Surface"))
{
    ui->setupUi(this);
    //setModel();

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setModel()
{
    //render->setModel();

    QString modelToLoad = QFileDialog::getOpenFileName(this, tr("Open 3D model to load"), ".");
    if(modelToLoad.isEmpty()) return;
    //render->model.loadModel(modelToLoad.toStdString());
    update();
    //render->model.printAllVertices();
}

//void MainWindow::on_actionCreate_Object_triggered(){
//    bool ok;
//    int height = QInputDialog::getInt(this, tr("height surface (in mm)"), tr("h:"), 1000, 0, 1500000, 1, &ok);
//    int width = QInputDialog::getInt(this, tr("width surface"), tr("w:"), 1000, 0, 1500000, 1, &ok);
//    int depth = QInputDialog::getInt(this, tr("depth surface"), tr("d:"), 1000, 0, 1500000, 1, &ok);
//    //QString filename = QFileDialog::getSaveFileName(this, tr("model"), ".obj");
//    //SurfaceModel(width,height,depth, filename::toStdString());
//}
