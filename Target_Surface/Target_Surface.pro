#-------------------------------------------------
#
# Project created by QtCreator 2015-06-29T16:42:29
#
#-------------------------------------------------

QT       += core gui
QT       += opengl
LIBS     += -L/home/camille/projects/tools/libsoil/unix64 -lSOIL \
            #-L/usr/local/lib -lglfw3 \
            -L/usr/lib/x86_64-linux-gnu -lGLEW \
            -L/usr/local/lib -lglfw3 -lpthread \
            -L/usr/lib/x86_64-linux-gnu -lglut -lGL -lGLU



#-L/usr/local/bar/libs -lGLEW  -lglut -lGL -lGLU

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Target_Surface
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    SurfaceMesh.cpp \
    SurfaceModel.cpp \
    rendering.cpp \
    utils.cpp

HEADERS  += mainwindow.h \
    SurfaceMesh.h \
    glm/glm.hpp \
    SurfaceModel.h \
    SOIL.h \
    global.h \
    glew/freeglut.h \
    glew/freeglut_ext.h \
    glew/freeglut_std.h \
    glew/gl.h \
    glew/glew.h \
    glew/glext.h \
    glew/gl_mangle.h \
    glew/glu.h \
    glew/glu_mangle.h \
    glew/glut.h \
    glew/glx.h \
    glew/glxew.h \
    glew/glxext.h \
    glew/glxint.h \
    glew/glx_mangle.h \
    glew/glxmd.h \
    glew/glxproto.h \
    glew/glxtokens.h \
    glew/wglew.h \
    glew/internal/dri_interface.h \
    glew/internal/glcore.h \
    rendering.h \
    GLFW/glfw3.h \
    GLFW/glfw3native.h \
    utils.h

FORMS    += mainwindow.ui

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/release/ -lassimp
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../usr/local/lib/debug/ -lassimp
else:unix: LIBS += -L$$PWD/../../../../../usr/local/lib/ -lassimp

INCLUDEPATH += $$PWD/../../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../../usr/local/include

QMAKE_CXXFLAGS += -std=c++11
