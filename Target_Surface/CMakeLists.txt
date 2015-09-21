cmake_minimum_required( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Release )
#add_definitions ( -Wall )
find_package ( Qt4 REQUIRED )


##########################
# C++11
#########################
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

##########################
# OpenGL
#########################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

####################################
# Ceres 
####################################
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

set(ASSIMP_FLAGS "-lassimp")

set(QT_USE_QTOPENGL TRUE)
set (CMAKE_CXX_FLAGS "-g -frounding-math -O3 ${ASSIMP_FLAGS}")
include ( ${QT_USE_FILE} )
include_directories (
	${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR} ${QT_QTOPENGL_INCLUDE_DIR}
	${OPENGL_INCLUDE_DIR}
	/usr/include/
	/usr/local/lib
)

set ( Target_Surface_Headers
    mainwindow.h
    SurfaceMesh.h
    glm/glm.hpp
    SurfaceModel.h
    global.h
    costFunctor.h
    rendering.h
    utils.h
    targetoptimization.h
	)

set ( Target_Surface_Src
    main.cpp
    mainwindow.cpp
    SurfaceMesh.cpp
    SurfaceModel.cpp
    rendering.cpp
    utils.cpp
    targetoptimization.cpp
	)

set ( Target_Surface_UIS
	mainwindow.ui
	)
QT4_WRAP_UI(UIS ${Target_Surface_UIS})

set ( Target_Surface_MOCS
	rendering.h
	mainwindow.h
	)
QT4_WRAP_CPP(MOCS ${Target_Surface_MOCS})

add_executable ( Target_Surface ${Target_Surface_Src} ${UIS} ${RSCS} ${TRS} ${MOCS} )
target_link_libraries ( Target_Surface  ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${OPENGL_LIBRARIES} ${GLU_LIBRARY} ${CERES_LIBRARIES} ${ASSIMP_FLAGS})
