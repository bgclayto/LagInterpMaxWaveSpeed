include_directories(${CMAKE_BINARY_DIR}/source/)

#file(GLOB_RECURSE myModules ${CMAKE_BINARY_DIR}/../source/*.f90)

add_executable(lambda
	main.f90
	arbitrary_eos.f90
	)

add_test(NAME case1 COMMAND ${CMAKE_SOURCE_DIR}/tests/test1.sh
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_test(NAME case2 COMMAND ${CMAKE_SOURCE_DIR}/tests/test2.sh
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
