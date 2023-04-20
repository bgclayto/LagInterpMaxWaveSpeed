include_directories(${CMAKE_BINARY_DIR}/source/)

#file(GLOB_RECURSE myModules ${CMAKE_BINARY_DIR}/../source/*.f90)

add_executable(a.exe
	main.f90
	arbitrary_eos.f90
	)
