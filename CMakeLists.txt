# CmakeLists for ICStoolbox project:
# 0 : macros and settings
# 1 : Commons library
# 2 : Graphics library and executable
# 3 : Other projects (Navier stokes, Velocity extension...)



##############################################################
# 0 : MACROS AND SETTINGS

# 0.1 : General CMake settings
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
set(CMAKE_CXX_FLAGS   "-w -O3")
set(CMAKE_C_FLAGS   "-w -O3")
add_definitions("-w")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

project(libMedit)
file(                  GLOB_RECURSE sources src/*.c)
file(                  GLOB_RECURSE hea src/*.h)
add_library(           Medit3D SHARED ${sources})

IF(APPLE)
	INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
	FIND_LIBRARY(COCOA_LIBRARY Cocoa)
	FIND_LIBRARY(GLUT_LIBRARY GLUT )
	FIND_LIBRARY(OpenGL_LIBRARY OpenGL )
	MARK_AS_ADVANCED (COCOA_LIBRARY
		GLUT_LIBRARY
	OpenGL_LIBRARY)
	SET(GL_LIBS ${COCOA_LIBRARY} ${GLUT_LIBRARY} ${OpenGL_LIBRARY})

elseif(UNIX)
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
	set(GL_LIBS ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )

ENDIF()

target_link_libraries( Medit3D ${GL_LIBS})
INSTALL(               TARGETS Medit3D LIBRARY  DESTINATION "$ENV{HOME}/lib")
INSTALL(               FILES ${hea} DESTINATION "$ENV{HOME}/include")

project(medit)
add_executable(        medit src/main.c)
target_link_libraries( medit Medit3D ${GL_LIBS})
INSTALL(               TARGETS medit RUNTIME DESTINATION "$ENV{HOME}/bin")
