project(CLAPACK_Library C CXX)

# 1. Configurações de compilação (equivalente ao makefile.defs)
set(CMAKE_C_STANDARD 99)
add_definitions(${CPPFLAGS}) # Caso existam definições extras no makefile.defs

# 2. Lista exaustiva de arquivos fonte
FILE(GLOB CLAPACK_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# 3. Criação da biblioteca estática (libclapack.a)
add_library(libclapack STATIC ${CLAPACK_SOURCES})

# 4. Inclusão de diretórios (Equivalente ao comportamento do Make)
target_include_directories(libclapack PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ../../melder *.h ../../melder/*.h)

message(STATUS "Saindo do diretório ${CMAKE_CURRENT_SOURCE_DIR}")
