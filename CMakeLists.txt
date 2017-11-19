cmake_minimum_required(VERSION 2.8.3)
PROJECT (python_cr3bp)

SET(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")

FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(PythonLibs REQUIRED)
if (${PYTHON_VERSION_MAJOR} EQUAL 3)
    FIND_PACKAGE(Boost COMPONENTS python3)
else()
    FIND_PACKAGE(Boost COMPONENTS python)
endif()

FIND_PACKAGE(PythonLibs 2.7 REQUIRED)
include_directories(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIR})

# Build and link the pylib module
add_library(pylib SHARED boost_cr3bp.cpp)
target_link_libraries(pylib ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

ENABLE_TESTING()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) # Deprecated but so convenient!

PYTHON_ADD_MODULE(boost_cr3bp boost_cr3bp.cpp)
FILE(COPY boost_cr3bp.py DESTINATION .)
ADD_TEST(NAME BOOST_CR3BP COMMAND ${PYTHON_EXECUTABLE} boost_cr3bp.py)

ADD_EXECUTABLE(BoostCr3bp boost_cr3bp.cpp)