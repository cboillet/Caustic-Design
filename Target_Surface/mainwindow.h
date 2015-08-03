#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "global.h"
#include <QMainWindow>
#include "SurfaceModel.h"
#include "rendering.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT


public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void setModel();

private:
    Ui::MainWindow *ui;
    Renderer* render;

protected slots:
    //void on_actionCreate_Object_triggered();
};

#endif // MAINWINDOW_H
