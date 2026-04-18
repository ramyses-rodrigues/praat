# Versão mínima do CMake
cmake_minimum_required(VERSION 3.10)
set(MODE "Release" CACHE STRING "Build mode: Debug or Release")
set(Exec "Praat_Mod")
SET(LIBS -lwinmm -lwsock32 -lcomctl32 -lole32 -lgdi32 -lgdiplus -lcomdlg32 -lwinspool 
	-static-libgcc -static-libstdc++ -mwindows -static -lwinpthread)

# Nome do projeto e linguagem
project(${Exec} LANGUAGES C CXX)

# Define o padrão C++ (opcional)
set(CMAKE_CXX_STANDARD 17)

include_directories( 
        fon 
        artsynth 
        FFNet 
		gram 
        EEG 
		LPC 
        dwtools 
        sensors 
		foned        
        stat 
		dwsys 
        sys 
        melder 
        kar 
		external/espeak 
		external/portaudio 
		external/flac
        external/lame
        external/mp3 
		external/glpk 
		external/clapack 
		external/gsl 
		external/num 
		external/vorbis 
		external/opusfile 
		external/whispercpp
        )

set(LIST_DIR_PROJECT 
		fon 
        artsynth 
        FFNet 
		gram 
        EEG 
		LPC 
        dwtools 
        sensors 
		foned        
        stat 
		dwsys 
        sys 
        melder 
        kar 
		external/espeak 
		external/portaudio 
		external/flac
        external/lame
        external/mp3 
		external/glpk 
		external/clapack 
		external/gsl 
		external/num 
		external/vorbis 
		external/opusfile 
		external/whispercpp
		)

foreach(FILE IN LISTS LIST_DIR_PROJECT)
    file(GLOB ${FILE}_SCRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${FILE}/*.cpp ${FILE}/*.c)
	list(LENGTH  ${FILE}_SCRS MY_LIST_LEN)
	message(STATUS "Variavel ${FILE}_SCRS possui ${MY_LIST_LEN} arquivos c/cpp")
endforeach()

# file(GLOB artsynth_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} artsynth/*.cpp)
# list(LENGTH artsynth_SRCS MY_LIST_LEN)
# message(STATUS "Variavel artsynthFiles: ${artsynth_SRCS} - ${MY_LIST_LEN} elementos")

# file(GLOB dwsys_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} dwsys/*.cpp)
# list(LENGTH dwsys_SRCS MY_LIST_LEN)
# message(STATUS "Variavel dwsysFiles: ${dwsys_SRCS} - ${MY_LIST_LEN} elementos")


# Setup output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)


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