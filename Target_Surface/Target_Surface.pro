#-------------------------------------------------
#
# Project created by QtCreator 2015-06-29T16:42:29
#
#-------------------------------------------------

QT       += core gui
QT       += opengl
LIBS     += -L/usr/local/bar/libs -lglut -lGL -lGLU -lGLEW -L/home/camille/projects/tools/SOIL/lib

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Target_Surface
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    SurfaceMesh.cpp \
    SurfaceModel.cpp

HEADERS  += mainwindow.h \
    SurfaceMesh.h \
    glm/glm.hpp \
    SurfaceModel.h \
    SOIL.h

FORMS    += mainwindow.ui
