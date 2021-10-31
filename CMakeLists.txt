cmake_minimum_required(VERSION 2.8)
project( 3dgraphics )
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
add_executable(3dgraphics 3dengine.cpp )
target_link_libraries(3dgraphics ${OpenCV_LIBS} )


#cmake_minimum_required(VERSION 2.8)
#project( Record )
#find_package( OpenCV REQUIRED )
#include_directories( ${OpenCV_INCLUDE_DIRS} )
#add_executable( Fourier_data data.cpp )
#target_link_libraries( Fourier_data ${OpenCV_LIBS} )
