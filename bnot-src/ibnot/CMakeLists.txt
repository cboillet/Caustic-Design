# CMake Project file for image BN-OT
# Fernando de Goes (fdegoes@caltech.edu)
# Copyright @ Aug 2012
project (ibnot)

cmake_minimum_required(VERSION 2.4.5)
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
if( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
set( CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} -framework Accelerate" )

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(CGAL COMPONENTS Qt4)
find_package(OpenGL)
find_package(Qt4)

set( SUITESPARSE_LIBRARIES 
     -lm -lamd -lcamd -lcolamd -lccolamd -lcholmod -lspqr -ltbb -lmetis -lblas -llapack)

set( HDRS
     matrix/sparse_array.h matrix/sparse_matrix.h matrix/suite_sparse_qr.h
     console_color.h convex_polygon.h dialog.h domain.h enriched_segment.h 
     glviewer.h grid.h line_search.h pgm.h pixel.h primitives.h pw_line_search.h
     ramp.h random.h rt2.h scene.h timer.h types.h util.h window.h )

set( SRCS 
     matrix/sparse_array.cpp matrix/sparse_matrix.cpp
     assign.cpp energy.cpp glviewer.cpp histogram.cpp init.cpp 
     io.cpp main.cpp optimizer.cpp regularity.cpp render.cpp sites.cpp window.cpp )

set( MOCS moc_dialog.cxx moc_glviewer.cxx moc_window.cxx )

if( CGAL_FOUND AND QT4_FOUND AND OPENGL_FOUND AND CGAL_Qt4_FOUND )
    include(${QT_USE_FILE})
    include(${CGAL_USE_FILE})

    # UI files (Qt Designer files)
    qt4_wrap_ui( DT_UI_FILES dialog.ui ibnot.ui )

    # use the Qt MOC preprocessor on classes that derives from QObject
    qt4_generate_moc( dialog.h moc_dialog.cxx )
    qt4_generate_moc( window.h moc_window.cxx )
    qt4_generate_moc( glviewer.h moc_glviewer.cxx )

    # The executable itself.
    add_executable( ibnot ${SRCS} ${MOCS} ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

    # Link with Qt libraries
    target_link_libraries( ibnot ${QT_LIBRARIES} )
    
    # Link with Glut and OpenGL
    target_link_libraries( ibnot ${GLUT_LIBRARY} ${OPENGL_LIBRARY} )

    # Link with CGAL
    target_link_libraries( ibnot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
    
    # Link with SuiteSparse
    target_link_libraries( ibnot ${SUITESPARSE_LIBRARIES} )

    include_directories(.)
else()
    set( OT_MISSING_DEPS "")

    if(NOT OPENGL_FOUND)
        set(OT_MISSING_DEPS "OpenGL, ${OT_MISSING_DEPS}")
    endif()

    if(NOT QT4_FOUND)
        set(OT_MISSING_DEPS "Qt4, ${OT_MISSING_DEPS}")
    endif()

    if(NOT CGAL_FOUND)
        set(OT_MISSING_DEPS "the CGAL library, ${OT_MISSING_DEPS}")
    endif()

    if(NOT CGAL_Qt4_FOUND)
        set(OT_MISSING_DEPS "the CGAL Qt4 library, ${OT_MISSING_DEPS}")
    endif()

    message(STATUS "NOTICE: This project requires ${OT_MISSING_DEPS} and will not be compiled.")
endif()
