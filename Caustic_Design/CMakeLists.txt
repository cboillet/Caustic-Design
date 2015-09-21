cmake_minimum_required( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Release )
#add_definitions ( -Wall )
find_package ( Qt4 REQUIRED )

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

################################
# lbfgs
################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")
FIND_PACKAGE(LBFGS)

################################
# suitesparse
################################
FIND_PACKAGE(SUITESPARSE REQUIRED)

###############################
# tinyxml
###############################
FIND_PACKAGE(TinyXML REQUIRED)


set(QT_USE_QTOPENGL TRUE)
set ( CGAL_FLAGS "-lCGAL -lboost_thread -lgmp -lmpfr -lm -lamd -lcamd -lcolamd -lccolamd -lcholmod -lspqr -ltbb -lmetis -lblas -llapack" )
set (CMAKE_CXX_FLAGS "-frounding-math -O3 -std=c++11 ${CGAL_FLAGS}")
include ( ${QT_USE_FILE} )
include_directories (
	${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR} ${QT_QTOPENGL_INCLUDE_DIR}
	${OPENGL_INCLUDE_DIR}
	${LBFGS_INCLUDE_DIRS}
	${SUITESPARSE_INCLUDE_DIRS}
	${TinyXML_INCLUDE_DIRS}
	/usr/include/)

set ( Caustic_Design_HDRS
	matrix/sparse_array.h
	matrix/sparse_matrix.h
	matrix/suite_sparse_qr.h
	console_color.h
	convex_polygon.h
	domain.h
	enriched_segment.h
	grid.h
	line_search.h
	pgm.h
	pixel.h
	primitives.h
	pw_line_search.h
	ramp.h
	random.h
	rt2.h
	scene.h
	timer.h
	types.h
	util.h
	interpolation.h
	optimal_transport.h
	voronoi_creation.h
	gradientdescent.h
	)

set ( Caustic_Design_SRCS
	matrix/sparse_array.cpp
	matrix/sparse_matrix.cpp
	assign.cpp
	energy.cpp
	glviewer.cpp
	histogram.cpp
	init.cpp
	io.cpp
	main.cpp
	optimizer.cpp
	regularity.cpp
	render.cpp
	sites.cpp
	window.cpp
	interpolation.cpp
	optimal_transport.cpp
	voronoi_creation.cpp
	gradientdescent.cpp
	)

set ( Caustic_Design_UIS
	caustic.ui
	dialog.ui
	)
QT4_WRAP_UI(UIS ${Caustic_Design_UIS})

set ( Caustic_Design_MOCS
	dialog.h
	glviewer.h
	window.h
	)
QT4_WRAP_CPP(MOCS ${Caustic_Design_MOCS})

add_executable ( Caustic_Design ${Caustic_Design_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )

target_link_libraries ( Caustic_Design  ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${OPENGL_LIBRARIES} ${GLU_LIBRARY} ${CGAL_FLAGS} ${LBFGS_LIBRARIES} ${SUITESPARSE_LIBRARIES} ${TinyXML_LIBRARIES})
