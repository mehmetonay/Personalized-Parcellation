#Umut Küçükaslan
#Mehmet Said Onay

cmake_minimum_required(VERSION 2.8.12)
 
project(Vavlab)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#find_package(VTK REQUIRED)
#include(${VTK_USE_FILE})

add_executable(Vavlab main.cxx Parcellation.cpp Parcellation.h)
target_link_libraries(Vavlab  ${ITK_LIBRARIES} Eigen3::Eigen)
