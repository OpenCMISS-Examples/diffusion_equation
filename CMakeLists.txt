cmake_minimum_required(VERSION 3.26 FATAL_ERROR)

if(DEFINED OpenCMISS_INSTALL_ROOT)
  if(EXISTS "${OpenCMISS_INSTALL_ROOT}/libOpenCMISSSetup.cmake")
    set(CMAKE_PROJECT_INCLUDE_BEFORE "${OpenCMISS_INSTALL_ROOT}/libOpenCMISSSetup.cmake")
  endif()
endif()  

set(ExampleName "diffusion_equation")

project(${ExampleName} VERSION 1.0 LANGUAGES C Fortran)

find_package(libOpenCMISS 1.0.0 REQUIRED
  COMPONENTS Fortran
)

add_subdirectory(src/fortran)
