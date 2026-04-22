cmake_minimum_required(VERSION 3.10)
project(Praat C CXX)

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
set(CMAKE_CXX_EXTENSIONS ON) # Necessário para o -std=gnu++17 (M_PI)

# 3. Tratamento do Ícone (Windows)
# O comando windres converte o .rc em um objeto que o linker entende
set(RESOURCE_FILES "main/praat_win.rc") # Presumindo que exista um .rc
if(WIN32)
    enable_language(RC)
    set(ICON_RESOURCE "main/praat_win.rc")
endif()

# 4. Adicionar Subdiretórios (Conforme vimos antes)
add_subdirectory(external/clapack)
# ... adicione os outros aqui ...
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

# 6. Linkagem de Bibliotecas do Sistema (Substitui LIBS)
target_link_libraries(${EXECUTABLE_NAME}
    PRIVATE
    # Suas bibliotecas compiladas
    fon foned dwtools LPC EEG sensors gram FFNet artsynth
    stat dwsys sys melder kar
    # Bibliotecas externas
    clapack gsl glpk lame mp3 num flac portaudio espeak vorbis opusfile whisper
    # Bibliotecas do Windows (LIBS)
    winmm wsock32 comctl32 ole32 gdi32 gdiplus comdlg32 winspool winpthread
    -static-libgcc -static-libstdc++ -static -mwindows
)


# cmake_minimum_required(VERSION 3.10)
# project(Praat)

# # 1. Configurações Globais (Substitui o makefile.defs)
# set(CMAKE_CXX_STANDARD 17)
# add_compile_options(-O3) # Ou as flags que estavam no seu .defs

# # 2. Adicionar todos os subdiretórios das bibliotecas externas
# # Cada pasta deve ter seu próprio CMakeLists.txt com add_library()
# add_subdirectory(external/clapack)


# # 3. Adicionar os subdiretórios do núcleo (all-self)
# add_subdirectory(kar)
# add_subdirectory(melder)
# add_subdirectory(sys)
# add_subdirectory(dwsys)
# add_subdirectory(stat)
# add_subdirectory(fon)
# add_subdirectory(foned)
# add_subdirectory(dwtools)
# add_subdirectory(LPC)
# add_subdirectory(EEG)
# add_subdirectory(sensors)
# add_subdirectory(gram)
# add_subdirectory(FFNet)
# add_subdirectory(artsynth)
# add_subdirectory(main)

# # 4. Definir o Executável Principal
# # O arquivo main_Praat.o (ou .cpp) deve estar em main/
# add_executable(praat main/main_Praat.cpp)

# # 5. Linkar todas as bibliotecas (Substitui o bloco $(LINK))
# target_link_libraries(praat
#     PRIVATE
#     fon foned dwtools LPC EEG sensors gram FFNet artsynth
#     stat dwsys sys melder kar
#     espeak portaudio flac lame mp3 glpk clapack gsl num vorbis opusfile whisper
#     ${CMAKE_DL_LIBS} # Substitui flags de sistema como -ldl se necessário
# )



# cmake_minimum_required(VERSION 3.15)
# set(CMAKE_VERBOSE_MAKEFILE ON)
# if(CMAKE_GENERATOR MATCHES "Ninja")
#     if(DEFINED CMAKE_GENERATOR_PLATFORM)
#         unset(CMAKE_GENERATOR_PLATFORM CACHE)
#     endif()
# endif()

# # Setup output directories.
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# project(Praat_Mod VERSION 1.0 LANGUAGES C CXX)

# # --- Configurações de Compilação (makefile.defs) ---

# # Opção para alternar entre Debug e Release[span_0](end_span)
# option(DEBUG_ON "Ativar modo de depuração" OFF)

# set(CMAKE_C_STANDARD 99)
# set(CMAKE_CXX_STANDARD 17)

# # Definições comuns (COMMONFLAGS)
# add_definitions(-municode -D_FILE_OFFSET_BITS=64 -DWORDS_BIGENDIAN=0 -DPA_LITTLE_ENDIAN)

# if(DEBUG_ON)
#     set(DEBUG_TAG "_Debug")
#     set(CMAKE_BUILD_TYPE Debug)
#     add_compile_options(-g)
# else()
#     set(DEBUG_TAG "_Release")
#     set(CMAKE_BUILD_TYPE Release)
#     add_compile_options(-O3) 
# endif()

# # Nome do executável conforme definido no Makefile
# set(EXECUTABLE_NAME "Praat_Mod${DEBUG_TAG}")

# # --- Bibliotecas Externas e Internas (Makefile) ---

# # Nota: No CMake, o ideal é que cada subdiretório tenha seu próprio CMakeLists.txt.
# # Aqui definimos a ordem de inclusão baseada nos comandos 'all-external' e 'all-self'.

# set(SUBDIRS_EXTERNAL
#     external/clapack external/gsl external/glpk external/lame external/mp3
#     external/num external/flac external/portaudio external/espeak 
#     external/vorbis external/opusfile external/whispercpp
# )

# set(SUBDIRS_SELF
#     kar melder sys dwsys stat fon foned dwtools 
#     LPC EEG sensors gram FFNet artsynth main
# )

# foreach(dir ${SUBDIRS_EXTERNAL} ${SUBDIRS_SELF})
# 	message(STATUS "Adicionando subdiretório ${dir}")
# 	add_subdirectory(${dir})
# endforeach()

# # --- Definição do Executável e Links ---

# # O arquivo principal e o ícone
# add_executable(${EXECUTABLE_NAME} main/main_Praat.o main/praat_win.o)

# # Linkagem das bibliotecas estáticas geradas nos subdiretórios
# target_link_libraries(${EXECUTABLE_NAME} PRIVATE
#     fon artsynth FFNet gram EEG LPC dwtools sensors foned stat
#     dwsys sys melder kar espeak portaudio flac lame mp3
#     glpk clapack gsl num vorbis opusfile whisper
#     # Bibliotecas do sistema Windows (LIBS)
#     winmm wsock32 comctl32 ole32 gdi32 gdiplus comdlg32 winspool
#     winpthread
# )

# # Flags de linkagem estática e interface Windowss
# set_target_properties(${EXECUTABLE_NAME} PROPERTIES 
#     LINK_FLAGS "-static-libgcc -static-libstdc++ -mwindows -static"
# )

# file(GLOB FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${dir}/*.cpp ${dir}/*.c)	
	
	# foreach(file ${FILES})
	# 	message(STATUS "Arquivo a ser acrescido: ${file}")
	# 	list(LENGTH FILES MY_LIST_LEN)
	# endforeach()

	# # message(STATUS "Variavel ${DIR}_SCRS possui ${MY_LIST_LEN} arquivos c/cpp")
	
	# get_filename_component(path "${dir}" NAME)
	# set(LIBa "lib${path}")
	# file(REMOVE LIBa)
	# message(STATUS "Lib gerada: ${LIBa}")
	
	# add_library(${LIBa} STATIC ${DIR_SCRS})


# Adiciona o executável e os arquivos fonte
# add_executable(${Exec} 
#         fon/libfon.a 
#         artsynth/libartsynth.a 
#         FFNet/libFFNet.a 
# 		gram/libgram.a 
#         EEG/libEEG.a 
# 		LPC/libLPC.a 
#         dwtools/libdwtools.a 
#         sensors/libsensors.a 
# 		foned/libfoned.a        
#         stat/libstat.a 
# 		dwsys/libdwsys.a 
#         sys/libsys.a 
#         melder/libmelder.a 
#         kar/libkar.a 
# 		external/espeak/libespeak.a 
# 		external/portaudio/libportaudio.a 
# 		external/flac/libflac.a 
#         external/lame/liblame.a 
#         external/mp3/libmp3.a 
# 		external/glpk/libglpk.a 
# 		external/clapack/libclapack.a 
# 		external/gsl/libgsl.a 
# 		external/num/libnum.a 
# 		external/vorbis/libvorbis.a 
# 		external/opusfile/libopusfile.a 
# 		external/whispercpp/libwhisper.a
#         )





# File: Makefile

# Makefile for Praat
# Paul Boersma & David Weenink, 7 July 2025

# System-dependent definitions of CC, LIBS, ICON and MAIN_ICON should be in
# makefile.defs, which has to be copied and renamed
# from a suitable makefile.defs.XXX file in the makefiles directory,
# Perhaps that file requires some editing.
# include makefile.defs

# .PHONY: all clean install

# # Makes the Praat executable in the source directory.
# all: all-external all-self
# 	$(LINK) -o $(EXECUTABLE) main/main_Praat.o $(MAIN_ICON) fon/libfon.a \
# 		artsynth/libartsynth.a FFNet/libFFNet.a \
# 		gram/libgram.a EEG/libEEG.a \
# 		LPC/libLPC.a dwtools/libdwtools.a sensors/libsensors.a \
# 		foned/libfoned.a fon/libfon.a stat/libstat.a \
# 		dwsys/libdwsys.a sys/libsys.a melder/libmelder.a kar/libkar.a \
# 		external/espeak/libespeak.a \
# 		external/portaudio/libportaudio.a \
# 		external/flac/libflac.a external/lame/liblame.a external/mp3/libmp3.a \
# 		external/glpk/libglpk.a \
# 		external/clapack/libclapack.a \
# 		external/gsl/libgsl.a \
# 		external/num/libnum.a \
# 		external/vorbis/libvorbis.a \
# 		external/opusfile/libopusfile.a \
# 		external/whispercpp/libwhisper.a \
# 		$(LIBS)

# all-external:
# 	$(MAKE) -C external/clapack
# 	$(MAKE) -C external/gsl
# 	$(MAKE) -C external/glpk
# 	$(MAKE) -C external/lame
# 	$(MAKE) -C external/mp3
# 	$(MAKE) -C external/num
# 	$(MAKE) -C external/flac
# 	$(MAKE) -C external/portaudio
# 	$(MAKE) -C external/espeak
# 	$(MAKE) -C external/vorbis
# 	$(MAKE) -C external/opusfile
# 	$(MAKE) -C external/whispercpp

# all-self:
# 	$(MAKE) -C kar
# 	$(MAKE) -C melder
# 	$(MAKE) -C sys
# 	$(MAKE) -C dwsys
# 	$(MAKE) -C stat
# 	$(MAKE) -C fon
# 	$(MAKE) -C foned
# 	$(MAKE) -C dwtools
# 	$(MAKE) -C LPC
# 	$(MAKE) -C EEG
# 	$(MAKE) -C sensors
# 	$(MAKE) -C gram
# 	$(MAKE) -C FFNet
# 	$(MAKE) -C artsynth
# 	$(MAKE) -C main main_Praat.o $(ICON)

# clean: clean-external clean-self
# 	$(RM) praat

# clean-external:
# 	$(MAKE) -C external/clapack clean
# 	$(MAKE) -C external/gsl clean
# 	$(MAKE) -C external/glpk clean
# 	$(MAKE) -C external/lame clean
# 	$(MAKE) -C external/mp3 clean
# 	$(MAKE) -C external/num clean
# 	$(MAKE) -C external/flac clean
# 	$(MAKE) -C external/portaudio clean
# 	$(MAKE) -C external/espeak clean
# 	$(MAKE) -C external/vorbis clean
# 	$(MAKE) -C external/opusfile clean
# 	$(MAKE) -C external/whispercpp clean

# clean-self:
# 	$(MAKE) -C kar clean
# 	$(MAKE) -C melder clean
# 	$(MAKE) -C sys clean
# 	$(MAKE) -C dwsys clean
# 	$(MAKE) -C stat clean
# 	$(MAKE) -C fon clean
# 	$(MAKE) -C foned clean
# 	$(MAKE) -C dwtools clean
# 	$(MAKE) -C LPC clean
# 	$(MAKE) -C EEG clean
# 	$(MAKE) -C sensors clean
# 	$(MAKE) -C gram clean
# 	$(MAKE) -C FFNet clean
# 	$(MAKE) -C artsynth clean
# 	$(MAKE) -C main clean

# install:
# 	$(INSTALL)