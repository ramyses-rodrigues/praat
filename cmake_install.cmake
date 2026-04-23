# Install script for directory: C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/Praat")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "C:/msys64/ucrt64/bin/objdump.exe")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/clapack/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/gsl/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/glpk/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/lame/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/mp3/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/num/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/flac/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/portaudio/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/espeak/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/vorbis/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/opusfile/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/external/whispercpp/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/kar/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/melder/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/sys/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/dwsys/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/stat/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/fon/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/foned/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/dwtools/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/LPC/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/EEG/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/sensors/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/gram/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/FFNet/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/artsynth/cmake_install.cmake")
  include("C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/main/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "C:/Users/ramyses.rmr/OneDrive/Documentos/Projects/Praat/Praat/praat/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
