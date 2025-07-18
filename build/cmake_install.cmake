# Install script for directory: /home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build2")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/generic_dummy/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/portlib/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/vaxonly/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/lsode_linpack/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/lsode/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/comput/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/smlib/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/fluxav/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/r8bloat/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/mdstransp/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/esc/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/ezcdf/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/nscrunch/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/pspline/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/geqdsk_mds/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/xplasma2/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/old_xplasma/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/i2mex/cmake_install.cmake")
  include("/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/pest3/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/stubenj9/ROCK/PEST3/PEST3_OMFIT_BUILD/pest3code/code/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
