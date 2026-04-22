message(STATUS "Entrando no diretório ${CMAKE_CURRENT_SOURCE_DIR}")

cmake_minimum_required(VERSION 3.10)
project(ClapackLibrary C)

# Define o diretório de inclusão (equivalente ao CPPFLAGS)
include_directories(../../melder)

# Define os arquivos de origem (OBJECTS)
set(SOURCES
    blas.c
    lapack.c
    lapack_dg.c
    lapack_dlaq.c
    lapack_dlar.c
    lapack_ds.c
    lapack_dt.c
)

# Cria a biblioteca estática libclapack.a
add_library(libclapack STATIC ${SOURCES})

# Garante que a biblioteca dependa dos headers (equivalente a $(OBJECTS): *.h)
# O CMake já monitora mudanças nos headers incluídos nos fontes automaticamente,
# mas você pode listar os headers explicitamente se desejar.




# # Exemplo para o diretório 'clapack'. 
# #Repita este padrão para os outros subdiretórios (gsl, glpk, lame, mp3, num, flac, portaudio, espeak, vorbis, opusfile, whispercpp).   

# project(libclapack)
# # message("CMAKE_PROJECT_NAME = ${CMAKE_PROJECT_NAME}")
# # message("PROJECT_NAME = ${PROJECT_NAME}")

# # Coleta todos os arquivos .cpp e .c do diretório
# file(GLOB CLAPACK_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
# foreach(file ${CLAPACK_SOURCES})
#     get_filename_component(filename ${file} NAME)
#     message(STATUS "Arquivo encontrado: ${filename}")
# endforeach()

# set(LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
# set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# # message(STATUS "Diretório de saída da lib: ${LIBRARY_OUTPUT_DIRECTORY}")

# # Cria a biblioteca estática (ex: libclapack.a)
# message(STATUS "Criando biblioteca ${PROJECT_NAME}.a")
# # message(STATUS "Fontes da biblioteca ${CLAPACK_SOURCES}")
# add_library(${PROJECT_NAME} STATIC ${CLAPACK_SOURCES})

# # Garante que as flags globais sejam aplicadas
# target_compile_options(${PROJECT_NAME} PRIVATE ${DEBUGFLAGS} ${COMMONFLAGS})

# target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${PROJECT_NAME})

message(STATUS "Saindo do diretório ${CMAKE_CURRENT_SOURCE_DIR}")
