cmake_minimum_required(VERSION 2.8.4)
project(cg1_ex5)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")


########################################################
#Options go here
########################################################
option(WITH_OPENMP "use openMP" OFF)
option(WITH_WARNINGS "display warnings" OFF)
option(WITH_OPENIMAGEIO "use oiio" OFF)

#INCLUDE_DIRECTORIES($(ex5_SOURCE_DIR)/glm)
#add_subdirectory(glm)

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

#########################################################
# FIND GLEW
#########################################################
find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
link_directories(${GLEW_LIBRARY_DIRS})
add_definitions(${GLEW_DEFINITIONS})
if(NOT GLEW_FOUND)
    message(ERROR " GLEW not found!")
endif(NOT GLEW_FOUND)


#########################################################
# FIND OpenImageIO
#########################################################
if(WITH_OPENIMAGEIO)
    find_package(OpenImageIO REQUIRED)
    include_directories(${OPENIMAGEIO_INCLUDE_DIRS})
    link_directories(${OPENIMAGEIO_LIBRARIES})
    add_definitions(-DWITH_OPENIMAGEIO)
endif()



#########################################################
# find OpenMP
#########################################################
if(WITH_OPENMP)
    FIND_PACKAGE(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        add_definitions(-DWITH_OPENMP)
    endif()
endif(WITH_OPENMP)


if(WITH_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

include_directories(
    glm/
    .
)

set(SOURCE_FILES
    Assets.cpp
    Assets.hpp
    Context.cpp
    Context.hpp
    GLSLShader.cpp
    GLSLShader.hpp
    Scene.cpp
    Scene.hpp
    TriMesh.cpp
    TriMesh.hpp
    View3D.cpp
    View3D.hpp
    ViewRender.cpp
    ViewRender.hpp
    Window.cpp
    Window.hpp
    Model.cpp
    Model.hpp
    Image.cpp
    Image.hpp
    Ray.hpp
    Ray.cpp)

add_executable(cg1_ex5 ${SOURCE_FILES})

target_link_libraries(cg1_ex5 ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${GLEW_LIBRARIES} ${GLU_LIBRARY} ${OPENIMAGEIO_LIBRARIES})
