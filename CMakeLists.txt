cmake_minimum_required(VERSION 3.5)

project(LagrangeTest)

enable_testing()

enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(dialect "-ffree-form -std=f2008 -fimplicit-none")
	set(bounds " -Wall -fbounds-check -g -C -fcheck=all -fbacktrace")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE " -O3 " CACHE STRING
	" Flags used by the compiler during release builds. ")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

ADD_CUSTOM_TARGET(debug
	COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
	COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
	COMMENT "Switch CMAKE_BUILD_TYPE to debug"
	)

ADD_CUSTOM_TARGET(release
	COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
	COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all 
	COMMENT "Switch CMAKE_BUILD_TYPE to Release"
	)

add_subdirectory(source)


