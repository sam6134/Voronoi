cmake_minimum_required(VERSION 3.1...3.15)
project(Kernel_23_Examples)

find_package(CGAL REQUIRED)
create_single_source_cgal_program("voronoi.cpp")