# Versão mínima do CMake
cmake_minimum_required(VERSION 3.10)

# set(MODE "RELEASE" CACHE STRING "Build mode: Debug or Release")
set(MODE "DEBUG" CACHE STRING "Build mode: Debug or Release")

set(EXEC_NAME "Praat_Mod_{MODE}")

# Nome do projeto e linguagem
project(${EXEC_NAME} LANGUAGES C CXX)

SET(LIBS -lwinmm -lwsock32 -lcomctl32 -lole32 -lgdi32 -lgdiplus -lcomdlg32 -lwinspool 
	-static-libgcc -static-libstdc++ -mwindows -static -lwinpthread)

add_definitions(-DUNICODE -D_UNICODE)
# add_definitions(-D_CRT_SECURE_NO_WARNINGS)
# add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
# add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
# add_definitions(-D_USE_MATH_DEFINES)
# add_definitions(-DWIN32_LEAN_AND_MEAN)

# set(CMAKE_GENERATOR_PLATFORM x64 CACHE INTERNAL "Force 64-bit")

# if(MSVC)
#     target_compile_options(${EXEC_NAME} PRIVATE _SBCS)
# endif()

# Define o padrão C++ (opcional)
set(CMAKE_CXX_STANDARD 17)
# target_compile_options(${EXEC_NAME} PRIVATE -std=gnu++17 CMAKE_CL_64 Wshadow)
# set(COMMON_FLAGS -municode -D_FILE_OFFSET_BITS=64 -DWORDS_BIGENDIAN=0 -DPA_LITTLE_ENDIAN -O3)
# set(CMAKE_CXX_COMPILE_FEATURESFLAGS "${CMAKE_CXX_FLAGS} std=gnu++17")
# message(STATUS "Variavel: ${CMAKE_CXX_EXTENSIONS} ${CMAKE_CXX_EXTENSIONS}")

# Setup output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

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


foreach(DIR IN LISTS LIST_DIR_PROJECT)
    
	set(DIR_SCRS ${DIR}_SCRS)
	file(GLOB DIR_SCRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*.cpp ${DIR}/*.c)
	# message(STATUS "Variavel ${DIR}_SCRS: ${DIR_SCRS}")
	
	list(LENGTH  DIR_SCRS MY_LIST_LEN)
	# message(STATUS "Variavel ${DIR}_SCRS possui ${MY_LIST_LEN} arquivos c/cpp")
	
	get_filename_component(path "${DIR}" NAME)
	set(LIBa "lib${path}")
	file(REMOVE LIBa)
	message(STATUS "Lib gerada: ${LIBa}")
	
	add_library(${LIBa} STATIC ${DIR_SCRS})
endforeach()




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