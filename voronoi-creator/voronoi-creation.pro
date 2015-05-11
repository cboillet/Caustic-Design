#-------------------------------------------------
#
# Project created by QtCreator 2015-05-11T11:38:48
#
#-------------------------------------------------

QT       += core
QT       += gui

TARGET = voronoi-creation
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    sites.cpp \
    matrix/sparse_array.cpp \
    matrix/sparse_matrix.cpp \
    init.cpp \
    assign.cpp \
    voronoicreator.cpp \
    io.cpp \
    optimizer.cpp \
    energy.cpp

INCLUDEPATH +=   /usr/include/
LIBS        += -L/usr/include/
LIBS        += -lCGAL
LIBS        += -lboost_thread
LIBS        += -lgmp
LIBS        += -lmpfr # not really needed for me, but added since gmp had to be added too
# suitesparse libs
LIBS        += -lm
LIBS        += -lamd
LIBS        += -lcamd
LIBS        += -lcolamd
LIBS        += -lccolamd
LIBS        += -lcholmod
LIBS        += -lspqr
LIBS        += -ltbb
LIBS        += -lmetis
LIBS        += -lblas
LIBS        += -llapack
QMAKE_CXXFLAGS += -frounding-math -O3

HEADERS += \
    convex_polygon.h \
    domain.h \
    enriched_segment.h \
    grid.h \
    pixel.h \
    primitives.h \
    rt2.h \
    tst.h \
    types.h \
    util.h \
    scene.h \
    matrix/sparse_array.h \
    matrix/sparse_matrix.h \
    matrix/suite_sparse_qr.h \
    pgm.h \
    random.h \
    timer.h \
    console_color.h \
    voronoicreator.h \
    line_search.h \
    pw_line_search.h

