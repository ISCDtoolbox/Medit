
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
set(CMAKE_CXX_FLAGS   "-w -O3")
set(CMAKE_C_FLAGS   "-w -O3")
add_definitions("-w")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "$ENV{HOME}/lib")
set(CMAKE_INSTALL_NAME_DIR "$ENV{HOME}/lib")
#INCLUDE_DIRECTORIES(  "$ENV{HOME}/include")
#LINK_DIRECTORIES(     "$ENV{HOME}/lib")
#SET(CMAKE_SKIP_BUILD_RPATH  TRUE)



macro(DATE RESULT)
    if(WIN32)
        execute_process(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE ${RESULT} OUTPUT_STRIP_TRAILING_WHITESPACE)
        string(REGEX REPLACE "(..)/(..)/(....).*" "\\3-\\2-\\1" ${RESULT} ${${RESULT}})
    elseif(UNIX)
        execute_process(COMMAND "date" "+%Y-%m-%d" OUTPUT_VARIABLE ${RESULT} OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
        message(SEND_ERROR "Unable to detect date")
        set(${RESULT} UNKNOWN)
    endif()
endmacro()
macro(TIME RESULT)
    if(WIN32)
        execute_process(COMMAND "cmd" " /C echo %TIME%" OUTPUT_VARIABLE ${RESULT} OUTPUT_STRIP_TRAILING_WHITESPACE)
        string(REGEX REPLACE "(..:..:..),(..)" "\\1" ${RESULT} ${${RESULT}})
    elseif(UNIX)
        execute_process(COMMAND "date" "+%H:%M:%S" OUTPUT_VARIABLE ${RESULT} OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
        message(SEND_ERROR "Unable to detect time")
        set(${RESULT} UNKNOWN)
    endif()
endmacro()



project(libMedit)
#string(TIMESTAMP time "%Y-%m-%d %H:%M:%S")
DATE(day)
TIME(hour)
file(WRITE src/compil.date "#define COMPIL  \"${day} ${hour}\"")
file(                  GLOB_RECURSE sources src/*.c)
file(                  GLOB_RECURSE hea src/*.h)


IF(APPLE)
	add_library(           Medit3D STATIC ${sources})
	INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
	FIND_LIBRARY(COCOA_LIBRARY Cocoa)
	FIND_LIBRARY(GLUT_LIBRARY GLUT )
	FIND_LIBRARY(OpenGL_LIBRARY OpenGL )
	MARK_AS_ADVANCED (COCOA_LIBRARY
		GLUT_LIBRARY
	OpenGL_LIBRARY)
	SET(GL_LIBS ${COCOA_LIBRARY} ${GLUT_LIBRARY} ${OpenGL_LIBRARY})

	target_link_libraries( Medit3D ${GL_LIBS})

elseif(UNIX)
	add_library(           Medit3D SHARED ${sources})
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

	target_link_libraries( Medit3D ${GL_LIBS})
	INSTALL(               TARGETS Medit3D LIBRARY  DESTINATION "$ENV{HOME}/lib")
	INSTALL(               FILES ${hea} DESTINATION "$ENV{HOME}/include")

ENDIF()

target_link_libraries( Medit3D ${GL_LIBS})


project(medit)
add_executable(        medit src/medit.c)
target_link_libraries( medit Medit3D)# ${GL_LIBS})
INSTALL(               TARGETS medit RUNTIME DESTINATION "$ENV{HOME}/bin")
