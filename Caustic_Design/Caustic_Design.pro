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
# L-BFGS
LIBS        += -llbfgs
QMAKE_CXXFLAGS += -frounding-math -O3

HEADERS += \
    matrix/sparse_array.h \
    matrix/sparse_matrix.h \
    matrix/suite_sparse_qr.h \
    console_color.h \
    convex_polygon.h \
    dialog.h \
    domain.h \
    enriched_segment.h \
    glviewer.h \
    grid.h \
    line_search.h \
    pgm.h \
    pixel.h \
    primitives.h \
    pw_line_search.h \
    ramp.h \
    random.h \
    rt2.h \
    scene.h \
    timer.h \
    types.h \
    util.h \
    window.h \
    interpolation.h \
	optimal_transport.h

SOURCES += \
    matrix/sparse_array.cpp \
    matrix/sparse_matrix.cpp \
    assign.cpp \
    energy.cpp \
    glviewer.cpp \
    histogram.cpp \
    init.cpp \
    io.cpp \
    main.cpp \
    optimizer.cpp \
    regularity.cpp \
    render.cpp \
    sites.cpp \
    window.cpp \
    interpolation.cpp \
    optimal_transport.cpp


FORMS += \
    caustic.ui \
    dialog.ui
QT += opengl
CONFIG += c++11
