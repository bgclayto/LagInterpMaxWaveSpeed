include_directories(${CMAKE_BINARY_DIR}/source/)

add_executable(lambda
	main.f90
	arbitrary_eos.f90
	)
