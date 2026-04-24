project(SYS_Library CXX)

# 1. Configurações de compilação (equivalente ao makefile.defs)
set(CMAKE_C_STANDARD 99)
add_definitions(${CPPFLAGS}) # Caso existam definições extras no makefile.defs

# 2. Lista exaustiva de arquivos fonte
FILE(GLOB SYS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# 3. Criação da biblioteca estática (libsys.a)
add_library(sys STATIC ${SYS_SOURCES})

# # 4. Inclusão de diretórios (Equivalente ao comportamento do Make)
# target_include_directories(sys PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# 3. Incluir diretórios de busca (Equivalente ao CPPFLAGS -I)
target_include_directories(sys PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/
    ../kar 
    ../melder
    ../dwsys
)
message(STATUS "Saindo do diretório ${CMAKE_CURRENT_SOURCE_DIR}")
