cmake_minimum_required(VERSION 3.10)
project(Praat VERSION 1.0 LANGUAGES C CXX)

# 1. Gerenciamento de DEBUG/RELEASE (Substitui DEBUG_ON)
# O CMake já tem suporte nativo: use -DCMAKE_BUILD_TYPE=Debug ou Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(DEBUG_TAG "_Debug")
else()
    set(DEBUG_TAG "_Release")
endif()

# 2. Configurações de Compilação (Substitui COMMONFLAGS, CFLAGS, CXXFLAGS)
add_definitions(
    -D_FILE_OFFSET_BITS=64 
    -DWORDS_BIGENDIAN=0 
    -DPA_LITTLE_ENDIAN
    -D_UNICODE -DUNICODE # Equivalente ao -municode
)

# Flags globais para C e C++
add_compile_options(-municode -O3)

# Padrões de linguagem
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON) # Necessário para o -std=gnu++17 (M_PI)
set(CMAKE_VERBOSE_MAKEFILE OFF)
# set(CMAKE_CXX_FLAGS "-std=c++11")

# 3. Tratamento do Ícone (Windows)
# O comando windres converte o .rc em um objeto que o linker entende
set(RESOURCE_FILES "main/praat_win.rc") # Presumindo que exista um .rc
if(WIN32)
    enable_language(RC)
    set(ICON_RESOURCE "main/praat_win.rc")
endif()

# 4. Adicionar Subdiretórios (Conforme vimos antes)
add_subdirectory(external/clapack)
add_subdirectory(external/gsl)
add_subdirectory(external/glpk)
add_subdirectory(external/lame)
add_subdirectory(external/mp3)
add_subdirectory(external/num)
add_subdirectory(external/flac)
add_subdirectory(external/portaudio)
add_subdirectory(external/espeak)
add_subdirectory(external/vorbis)
add_subdirectory(external/opusfile)
add_subdirectory(external/whispercpp)

add_subdirectory(kar)
add_subdirectory(melder)
add_subdirectory(sys)
add_subdirectory(dwsys)
add_subdirectory(stat)
add_subdirectory(fon)
add_subdirectory(foned)
add_subdirectory(dwtools)
add_subdirectory(LPC)
add_subdirectory(EEG)
add_subdirectory(sensors)
add_subdirectory(gram)
add_subdirectory(FFNet)
add_subdirectory(artsynth)
add_subdirectory(main)

# 5. Definição do Executável
set(EXECUTABLE_NAME "Praat_Mod${DEBUG_TAG}")

add_executable(${EXECUTABLE_NAME} main/main_Praat.cpp ${ICON_RESOURCE})

set_target_properties( ${EXECUTABLE_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# 6. Linkagem de Bibliotecas do Sistema (Substitui LIBS)
target_link_libraries(${EXECUTABLE_NAME}
    PRIVATE
    # Suas bibliotecas compiladas
    fon foned dwtools LPC EEG sensors gram FFNet artsynth
    stat dwsys sys melder kar
    # Bibliotecas externas
    clapack gsl glpk lame mp3 num flac portaudio espeak vorbis opusfile whisper
    # Bibliotecas do Windows (LIBS)
    # winmm wsock32 comctl32 ole32 gdi32 gdiplus comdlg32 winspool winpthread
    # -static-libgcc -static-libstdc++ -static -mwindows

    -lwinmm -lwsock32 -lcomctl32 -lole32 -lgdi32 -lgdiplus -lcomdlg32 -lwinspool
      -static -lc++ -lc++abi -mwindows
)

# set(dir "../melder")
# FILE(GLOB MELDER_HEADERS ${dir}/*.h)
# foreach(file ${MELDER_HEADERS})
#     get_filename_component(filebase ${file} NAME)
#     message(STATUS "Found header: ${filebase}")
# endforeach()

# find_file(FILE dir/melder.h)
# message("FILE points to ${dir}. Status: ${FILE}")

# # 4. Adicionar Subdiretórios (Conforme vimos antes)
# target_include_directories(
#     external/clapack
#     external/gsl
#     external/glpk
#     external/lame   
#     external/mp3
#     external/num
#     external/flac
#     external/portaudio
#     external/espeak
#     external/vorbis
#     external/opusfile
#     external/whispercpp
#     kar
#     melder
#     sys
#     dwsys
#     stat
#     fon
#     foned
#     dwtools
#     LPC
#     EEG
#     sensors
#     gram
#     FFNet
#     artsynth
#     main
# )
