cmake_minimum_required(VERSION 3.10)
project(SYS_Library C CXX)

# 1. Configurações de compilação (equivalente ao makefile.defs)
set(CMAKE_C_STANDARD 99)
add_definitions(${CPPFLAGS}) # Caso existam definições extras no makefile.defs

# 2. Lista exaustiva de arquivos fonte
FILE(GLOB SYS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# 3. Criação da biblioteca estática (libsys.a)
add_library(libsys STATIC ${SYS_SOURCES})

# # 4. Inclusão de diretórios (Equivalente ao comportamento do Make)
# target_include_directories(libsys PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# 3. Incluir diretórios de busca (Equivalente ao CPPFLAGS -I)
target_include_directories(libsys PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}        # Pasta sys
    ${CMAKE_CURRENT_SOURCE_DIR}/../kar # Pasta kar
    ${CMAKE_CURRENT_SOURCE_DIR}/../melder
    ${CMAKE_CURRENT_SOURCE_DIR}/../dwsys
    ${CMAKE_CURRENT_SOURCE_DIR}/..        # Para permitir include "fon/Function.h"
)

message(STATUS "Saindo do diretório ${CMAKE_CURRENT_SOURCE_DIR}")
