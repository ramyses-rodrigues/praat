cmake_minimum_required(VERSION 3.10)
project(NUM_Library C)

# 1. Configurações de compilação (equivalente ao makefile.defs)
set(CMAKE_C_STANDARD 99)
add_definitions(${CPPFLAGS}) # Caso existam definições extras no makefile.defs

# 2. Lista exaustiva de arquivos fonte
FILE(GLOB NUM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# 3. Criação da biblioteca estática (libnum.a)
add_library(libnum STATIC ${NUM_SOURCES})

# 4. Inclusão de diretórios (Equivalente ao comportamento do Make)
target_include_directories(libnum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

message(STATUS "Saindo do diretório ${CMAKE_CURRENT_SOURCE_DIR}")
