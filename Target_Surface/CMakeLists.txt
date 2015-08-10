cmake_minimum_required( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Release )
#add_definitions ( -Wall )
find_package ( Qt4 REQUIRED )



##########################
# OpenGL, GLUT
#########################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
link_directories(${GLEW_LIBRARIES})

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
	${GLUT_INCLUDE_DIR}
	${GLEW_INDLUDE_DIRS}
	/usr/include/
	/usr/local/lib
)

set ( Target_Surface_Headers
    mainwindow.h
    SurfaceMesh.h
    glm/glm.hpp
    SurfaceModel.h
    SOIL.h
    global.h
    glew/freeglut.h
    glew/freeglut_ext.h
    glew/freeglut_std.h
    glew/gl.h
    glew/glew.h
    glew/glext.h
    glew/gl_mangle.h
    glew/glu.h
    glew/glu_mangle.h
    glew/glut.h
    glew/glx.h
    glew/glxew.h
    glew/glxext.h
    glew/glxint.h
    glew/glx_mangle.h
    glew/glxmd.h
    glew/glxproto.h
    glew/glxtokens.h
    glew/wglew.h
    glew/internal/dri_interface.h
    glew/internal/glcore.h
    rendering.h
    GLFW/glfw3.h
    GLFW/glfw3native.h
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
target_link_libraries ( Target_Surface  ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${GLU_LIBRARY} ${CERES_LIBRARIES} ${GLEW_LIBRARIES} ${ASSIMP_FLAGS})
