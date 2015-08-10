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

    connect(meshHeightLineEdit, SIGNAL(returnPressed()), this, SLOT(newMeshHeight()));
    connect(focalLengthLineEdit, SIGNAL(returnPressed()), this, SLOT(newFocalLength()));

    //setModel();
    //viewer = new Renderer(5, this, "Target Surface");
    optimizer = new TargetOptimization(viewer->model);
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
    viewer->model.meshes[0].exportVertices(fileName, viewer->model.surfaceSize);
}

void MainWindow::on_actionLoadModel_triggered()
{
    setModel();
}

void MainWindow::on_actionGenerateTriangles_triggered()
{
    viewer->model.modifyMesh();
    viewer->update();
}

void MainWindow::on_actionRunTargetOptimization_triggered()
{
    optimizer->runOptimization(viewer->model);
}

void MainWindow::on_actionLoadLightRayReceiverPosition_triggered()
{
    QString filename = QFileDialog::getOpenFileName(this, tr("Open .dat File"), ".");
    if(filename.isEmpty()) return;
    viewer->model.loadReceiverLightPoints(filename);
    viewer->sceneUpdate();
}

void MainWindow::on_actionAxis_toggled()
{
    viewer->toggleDrawAxis();
    viewer->update();
}

void MainWindow::on_actionNormals_toggled()
{
    viewer->toggleDrawNormals();
    viewer->update();
}

void MainWindow::newMeshHeight()
{
    std::string txt = meshHeightLineEdit->text().toStdString();
    bool ok;

    float newScale = meshHeightLineEdit->text().toFloat(&ok);
    if(!ok)
    {
        std::cerr << "illegal value" << std::endl;
        return;
    }

    meshWidthLineEdit->setText(meshHeightLineEdit->text());

    newScale *= 0.5; // we go from -0.5 to 0.5
    std::cout << "New Mesh Height: " << newScale << std::endl;

    viewer->model.rescaleMeshes(newScale);
    viewer->model.surfaceSize = newScale;

    viewer->update();
}

void MainWindow::newFocalLength()
{
    std::string txt = focalLengthLineEdit->text().toStdString();
    bool ok;

    float newLength = focalLengthLineEdit->text().toFloat(&ok);
    if(!ok)
    {
        std::cerr << "illegal value" << std::endl;
        return;
    }

    viewer->model.setFocalLength(newLength);
    viewer->sceneUpdate();

}

void MainWindow::on_actionExit_triggered()
{
    QApplication::quit();
}
