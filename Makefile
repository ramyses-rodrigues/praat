# File: Makefile

# Makefile for Praat
# Paul Boersma, 17 May 2026

##########################
#
# To make Praat from the command line on Windows or Linux, type
#
# make
#   will make Praat with default settings
#
# There are several settings that you can either add on the command line
# or define as an environment variable.
#
##########################
#
# PRAAT_ARCH
#
# determines the processor architecture for which Praat is built.
# The default is what GCC and Clang will use if an architecture is not specified, which is:
#   - armv8-a (simplest ARM64) for the MSYS2-CLANGARM64 shell (Windows) and the OrbStack ARM64 shell (Linux)
#   - x86-64 (i.e. x64-v1) for the MSYS2-CLANG64 shell (Windows) and the OrbStack x64 shell (Linux)
#   - pentium4 (i.e. basic i686) for the MSYS2-MINGW32 shell (Windows)
#   - armv6+fp for some Raspberry Pi computers (you can check this using: gcc -Q --help=target)
# Note that x64-v1 is the original 64-bit AMD/Intel chip;
# for post-2013 computers you'd opt for x64-v3 instead.
#
# make PRAAT_ARCH=arm64
#   will make Praat with the "armv8-a" setting for the architecture,
#   producing an ARM64 edition of Praat.
#
# make PRAAT_ARCH=x64v1
#   will make Praat with the "x86-64" setting for the architecture,
#   producing an x64-v1 edition of Praat.
#
# make PRAAT_ARCH=x64v3
#   will make Praat with the "x86-64-v3" setting for the architecture,
#   i.e. with support for AVX2, FMA, F16C and more,
#   producing an x64-v3 edition of Praat.
#
# make PRAAT_ARCH=i686
#   will make Praat with the "i686" setting for the architecture,
#   producing a 32-bits Windows edition of Praat.
#
# make PRAAT_ARCH=native
#   will make Praat for the architecture of the shell,
#   which will often produce a version of Praat that is fastest on your computer.
#
# Tips and tricks
#
#   To make PRAAT_ARCH=native the default, you can add an environment variable, e.g. you can put
#     export PRAAT_ARCH=native
#   in your ~/.bashrc file (or so). Once you have done that,
#     make
#   will from then on behave the same as if you typed
#     make PRAAT_ARCH=native
#   Note that an explicit command-line-specified setting, such as
#     make PRAAT_ARCH=x64v1
#   will override the environment variable.
#
##########################
#
# PRAAT_CHROME=1
#
# Add this in case you want to compile for Linux on a Chromebook (gives extra window titles).
#
##########################

ifeq ($(OS),Windows_NT)
  OPERATING_SYSTEM := Windows
else
  OPERATING_SYSTEM := $(shell uname -s)
  ifeq ($(OPERATING_SYSTEM),Linux)
    OPTIONAL_RASPBERRY_PI_MODEL := $(shell cat /proc/device-tree/model 2>/dev/null)
    ifeq ($(findstring Raspberry Pi,$(OPTIONAL_RASPBERRY_PI_MODEL)),Raspberry Pi)
      WE_HAVE_RASPBERRY_PI := 1
    endif
  endif
endif
$(info OPERATING SYSTEM: "$(OPERATING_SYSTEM)")

ifeq ($(PRAAT_ARCH),arm64)
  ARCH_COMPILER_FLAGS = -march=armv8-a
else ifeq ($(PRAAT_ARCH),x64v1)
  ARCH_COMPILER_FLAGS = -march=x86-64
else ifeq ($(PRAAT_ARCH),x64v3)
  ARCH_COMPILER_FLAGS = -march=x86-64-v3
else ifeq ($(PRAAT_ARCH),i686)
  ARCH_COMPILER_FLAGS = -march=i686
else ifeq ($(PRAAT_ARCH),s390x)
  ARCH_COMPILER_FLAGS = -march=arch11 -mtune=arch13
else ifeq ($(PRAAT_ARCH),native)
  ARCH_COMPILER_FLAGS = -march=native
else
  $(info No architecture given)
  ARCH_COMPILER_FLAGS =
endif
$(info ARCHITECTURE: given as "$(PRAAT_ARCH)", resolved as "$(ARCH_COMPILER_FLAGS)")

ifeq ($(OPERATING_SYSTEM),Windows)
ifeq (1,1)
  #
  # Determine the compiler: either GCC or Clang.
  #
  # First try: explicit setting of the compiler via argument or environment variable PRAAT_COMPILER.
  ifeq ($(PRAAT_COMPILER),gcc)
    WE_HAVE_GCC := 1
  else ifeq ($(PRAAT_COMPILER),clang)
    WE_HAVE_GCC := 0
  # Second try: default setting, depending on the name of the shell.
  else ifeq ($(findstring MINGW,$(MSYSTEM)),MINGW)
    WE_HAVE_GCC := 1
  else ifeq ($(findstring CLANG,$(MSYSTEM)),CLANG)
    WE_HAVE_GCC := 0
  else
    $(error Unknown compiler)
  endif

  #
  # Compiler-specific variables:
  # - CC is a standard variable for the command that will be used by implicit `make` rules to compile C files
  # - CXX is a standard variable for the command that will be used by implicit `make` rules to compile C++ files
  # - LINKER_COMMAND is our command to invoke the linker; it should match the compiler command
  # - NON_PRAAT_LIBRARIES is our listing of the Windows libraries and the GCC or Clang libraries
  #
  ifeq ($(WE_HAVE_GCC),1)
    CC := gcc
    CXX := g++
    LINKER_COMMAND := g++
    NON_PRAAT_LIBRARIES = -lwinmm -lwsock32 -lcomctl32 -lole32 -lgdi32 -lgdiplus -lcomdlg32 -lwinspool \
      -static-libgcc -static-libstdc++ -mwindows -static -lwinpthread
      # Note: winpthread is linked statically, because it may not be available on all users' computers;
      # for the same reason, GCC's C and C++ libraries are also statically linked.
  else
    CC := clang
    CXX := clang
    LINKER_COMMAND := clang
    NON_PRAAT_LIBRARIES = -lwinmm -lwsock32 -lcomctl32 -lole32 -lgdi32 -lgdiplus -lcomdlg32 -lwinspool \
      -static -lc++ -lc++abi -mwindows
      # Note: Clang's C++ libraries are not available on all users' computers, so they are linked statically.
  endif

  #
  # Compiler flags:
  # - CFLAGS is a standard variable for the settings that will be used by implicit `make` rules to compile C files
  # - CXXFLAGS is a standard variable for the settings that will be used by implicit `make` rules to compile C++ files
  # - SHARED_COMPILER_FLAGS is our variable for the settings that will be shared between C and C++ compilation
  #
  SHARED_COMPILER_FLAGS := -municode -D_FILE_OFFSET_BITS=64 \
    $(ARCH_COMPILER_FLAGS) -O3
    # Note: the leaves many settings implicit, such as -mwin32 and -mtune=generic;
    # also -m64 (for x64-v1, x64-v3 and arm64) and -march=xxx (if not set in ARCH_COMPILER_FLAGS).
  CFLAGS = -std=gnu99 $(SHARED_COMPILER_FLAGS)
  CXXFLAGS = -std=gnu++17 $(SHARED_COMPILER_FLAGS) -Wshadow
    # Note: gnu++17 instead of c++17 is necessary to define M_PI in external code.

  EXECUTABLE_FILE = Praat.exe

  AR = ar
  RANLIB = ranlib
  WINDRES = windres
  ICON = praat_win.o
  MAIN_ICON = main/praat_win.o
else
  # Use clang on Windows under MSYS2

  # This is the makefile for building all three Windows edition of Praat:
  # for ARM64 (aarch64), x86_64 ("Intel64") and i686 ("Intel32") processors.
  # We therefore assume that you have installed three MSYS2 shells: CLANGARM64, CLANG64 and CLANG32.
  #
  # For reasons of speed, all object files are compiled with a compiler for the native physical processor,
  # i.e. with /clangarm64/bin/clang or /clang64/bin/clang, depending on whether your hardware is ARM64 or x86_64.
  # This means we have to do cross-compilation, with explicit specification
  # of the target architecture (aarch64-w64-windows-gnu, x86_64-w64-windows-gnu, or i686-w64-windows-gnu)
  # as well as of sysroot (/clangarm64, /clang64, or /clang32), which is the folder where
  # the include, lib and bin directories reside.
  #
  # Linking, on the other hand, takes place with a target-native compiler,
  # so that the file `libclang_rt.builtins.a` can be found.
  # On an ARM64 processor, this will work correctly, because such a system can emulate Intel64 and Intel32 binaries.
  # On an x86_64 processor, we don't precisely know how an ARM64 target can be linked, because only Intel32 binaries
  # can be emulated. Please tell us when you figure out how to link an ARM64 target on an x86_64 processor;
  # perhaps it doesn't exist, in which case we must advise you to work on an ARM64 computer if
  # you want to produce all three Windows editions of Praat.
  #
  # For the ARM64 target, the tools ar, ranlib and windres can be in irregular places and have irregular names.

  #
  # This makefile relies on the environment variables MSYSTEM (the target you want to build for)
  # and PROCESSOR_IDENTIFIER (a string describing the actual physical hardware processor on your computer).
  # If either variable is not set on your computer, hard-code it, as explained in the next two paragraphs.
  #

  # If MSYSTEM is not set on your computer, then that is probably a mistake in the set-up of your MSYS2 shell;
  # in case of emergency, however, you can hard-code MSYSTEM by uncommenting one of the following three lines:
  #MSYSTEM := CLANGARM64
  #MSYSTEM := CLANG64
  #MSYSTEM := CLANG32
  $(info Making: the MSYSTEM is $(MSYSTEM))

  # If PROCESSOR_IDENTIFIER is not set on your computer, then hard-code it by uncommenting one of the following two lines:
  #PROCESSOR_IDENTIFIER := ARM64
  #PROCESSOR_IDENTIFIER := x86_64
  $(info Making: the processor identifier is: $(PROCESSOR_IDENTIFIER))

  #
  # What the fastest compiler is, depends on the host computer, i.e. on PROCESSOR_IDENTIFIER.
  #
  ifeq ($(findstring ARM,$(PROCESSOR_IDENTIFIER)),ARM)
    PHYSICAL_PROCESSOR := ARM64
  else ifeq ($(findstring arm,$(PROCESSOR_IDENTIFIER)),arm)
    PHYSICAL_PROCESSOR := ARM64
  else
    PHYSICAL_PROCESSOR := x86_64
  endif
  $(info Making: the physical processor is $(PHYSICAL_PROCESSOR))
  ifeq ($(PHYSICAL_PROCESSOR),ARM64)
    FASTEST_COMPILER := /clangarm64/bin/clang
  else
    FASTEST_COMPILER := /clang64/bin/clang
  endif
  $(info Making: the fastest compiler is $(FASTEST_COMPILER))

  #
  # What the compilation target is, depends on the target computer architecture, i.e. on MSYSTEM.
  #
  ifeq ($(MSYSTEM),CLANGARM64)
    TARGET = aarch64-w64-windows-gnu
    SYSROOT = /clangarm64
  else ifeq ($(MSYSTEM),CLANG64)
    TARGET = x86_64-w64-windows-gnu
    SYSROOT = /clang64
  else ifeq ($(MSYSTEM),CLANG32)
    #
    # Building Praat for 32-bit Windows with clang requires having the folder C:/msys2/clang32.
    # Around 2025 MSYS stopped supporting 32-bit Windows,
    # so we work with version 18, from 2024, which we happened to have installed in time.
    # Fortunately, version 21 of clang on ARM64 and Intel64 does still support cross-compilation for the 32-bit target.
    #
    TARGET = i686-w64-windows-gnu
    SYSROOT = /clang32
  endif
  $(info Making: the target is $(TARGET))
  $(info Making: sysroot is $(SYSROOT))
  CC := $(FASTEST_COMPILER) --target=$(TARGET) --sysroot=$(SYSROOT)
  CXX := $(FASTEST_COMPILER) --target=$(TARGET) --sysroot=$(SYSROOT)
  $(info Making: the compiler set-up is $(CXX))

  #
  # What the best tools are (ar, ranlib, windres) depends crucially on the target computer architecture, because
  # the tools are specific to the target computer architecture. If we have a choice,
  # we prefer cross-tooling, because that will be slightly faster, so we then also depend on the host computer.
  #
  ifeq ($(MSYSTEM),CLANGARM64)
    ifeq ($(PHYSICAL_PROCESSOR),ARM64)
      # If you're on a physical ARM64 host, you'll like to install ar/ranlib/windres for ARM64 this way:
      #    pacman -S mingw-w64-clang-aarch64-llvm-tools
      # The tools will be ARM64 executables that create ARM64 object code;
      # ranlib will be at /clangarm64/bin/llvm-ranlib.exe.
      TOOLS_PREFIX = /clangarm64/bin/llvm-
    else
      # If you're on a physical x86_64 host, you'll like to install ar/ranlib/windres for ARM64 this way:
      #    pacman -S mingw-w64-cross-mingwarm64-binutils
      # or even
      #    pacman -S mingw-w64-cross-toolchain
      # The tools will be x64_64 executables that create ARM64 object code;
      # ranlib will be at /opt/bin/aarch64-w64-mingw32-ranlib.exe.
      TOOLS_PREFIX = /opt/bin/aarch64-w64-mingw32-
    endif
  else ifeq ($(MSYSTEM),CLANG64)
    # ar/ranlib/windres are in /usr/bin, even if you are on an ARM64 host.
    # ar/ranlib/windres are also in /clang64/bin/, so either of the following may be correct.
    TOOLS_PREFIX =
    #TOOLS_PREFIX = /clang64/bin/
  else ifeq ($(MSYSTEM),CLANG32)
    TOOLS_PREFIX = /clang32/bin/
  endif
  $(info Making: ar/ranlib/windres are in $(TOOLS_PREFIX))

  # The linker should usually be native, so as to be able to find the native runtime library, which is at:
  #    /clangarm64/lib/clang/21/lib/windows/libclang_rt.builtins-aarch64.a
  #    /clang64/lib/clang/21/lib/windows/libclang_rt.builtins-x86_64.a
  #    /clang32/lib/clang/18/lib/windows/libclang_rt.builtins-i386.a

  ifeq ($(MSYSTEM),CLANGARM64)
    # Cross-linking on an x86_64 host apparently tries to find the runtime library at
    #    /clang<64|arm64>/lib/clang/21/lib/aarch64-windows-gnu/libclang_rt.builtins.a
    # but we don't know how to install one there, so we use native ARM64 clang for linking,
    # making it (for the moment) impossible to build Praat for ARM64 on an x86_64 computer.
    LINKER_COMMAND := /clangarm64/bin/clang.exe --target=$(TARGET) --sysroot=$(SYSROOT)
    # Once we do know how to cross-link, we'll be able to do something like
    #LINKER_COMMAND = $(FASTEST_COMPILER) --target=$(TARGET) --sysroot=$(SYSROOT)
  else ifeq ($(MSYSTEM),CLANG64)
    # Cross-linking on an ARM64 host apparently tries to find the runtime library at
    #    /clang<64|arm64>/lib/clang/21/lib/x86_64-windows-gnu/libclang_rt.builtins.a
    # but we don't know how to install one there, so we use native x86_64 clang for linking;
    # fortunately, we can still build Praat for x86_64 on an ARM64 computer,
    # because ARM64 computers can run x86_64-native code.
    LINKER_COMMAND := /clang64/bin/clang.exe --target=$(TARGET) --sysroot=$(SYSROOT)
    # Once we do know how to cross-link, we'll be able to do something like
    #LINKER_COMMAND = $(FASTEST_COMPILER) --target=$(TARGET) --sysroot=$(SYSROOT)
    # We tried without avail:
    #    pacman -S mingw-w64-cross-clang
  else ifeq ($(MSYSTEM),CLANG32)
    # Cross-linking on an ARM64 or x86_64 host apparently tries to find the runtime library at
    #    /clang<32|64|arm64>/lib/clang/21/lib/i686-windows-gnu/libclang_rt.builtins.a
    # but we don't know how to install one there, so we use native i686 clang for linking;
    # fortunately, we can still build Praat for i686 on an ARM64 or x86_64 computer,
    # because ARM64 and x86_64 computers can run i686-native code.
    LINKER_COMMAND := /clang32/bin/clang.exe --target=$(TARGET) --sysroot=$(SYSROOT)
    # Once we do know how to cross-link, we'll be able to do something like
    #LINKER_COMMAND = $(FASTEST_COMPILER) --target=$(TARGET) --sysroot=$(SYSROOT)
  else
    $(warning Unknown MSYS environment: $(MSYSTEM))
  endif
  $(info Making: the linker set-up is $(LINKER_COMMAND))

  #
  # Now starts the minimally needed part of the makefile.
  # If you just want to do native compilation for an x86_64 target on an x86_64 host,
  # you can erase all of the code above.
  #

  CC ?= clang
  CXX ?= clang
  TOOLS_PREFIX ?=
  LINKER_COMMAND ?= $(CXX)

  SHARED_COMPILER_FLAGS := $(ARCH_COMPILER_FLAGS) -municode -D_FILE_OFFSET_BITS=64 \
    -O3
  # Probably implicit: -m64 -mwin32 -march=... -mtune=generic

  CFLAGS := -std=gnu99 $(SHARED_COMPILER_FLAGS)

  # gnu++17 instead of c++17 is necessary to define M_PI in external code
  CXXFLAGS := -std=gnu++17 $(SHARED_COMPILER_FLAGS) -Wshadow

  EXECUTABLE_FILE = Praat.exe

  NON_PRAAT_LIBRARIES = -lwinmm -lwsock32 -lcomctl32 -lole32 -lgdi32 -lgdiplus -lcomdlg32 -lwinspool -static-libgcc -lc++ -lc++abi -mwindows -static

  AR := $(TOOLS_PREFIX)ar
  RANLIB := $(TOOLS_PREFIX)ranlib
  WINDRES := $(TOOLS_PREFIX)windres
  ICON := praat_win.o
  MAIN_ICON := main/praat_win.o
endif
else ifeq ($(OPERATING_SYSTEM),Linux)

  PKG_CONFIG ?= pkg-config

  ifeq ($(PRAAT_GRAPHICS),barren)
    GRAPHICS_COMPILER_FLAGS := -DNO_GRAPHICS
    GRAPHICS_LINKER_FLAGS :=
    EXECUTABLE_FILE = praat_barren
  else ifeq ($(PRAAT_GRAPHICS),nogui)
    GRAPHICS_COMPILER_FLAGS := -DNO_GUI `$(PKG_CONFIG) --cflags pangocairo`
    GRAPHICS_LINKER_FLAGS := `$(PKG_CONFIG) --libs pangocairo`
    EXECUTABLE_FILE = praat_nogui
  else
    GRAPHICS_COMPILER_FLAGS := `$(PKG_CONFIG) --cflags gtk+-3.0`
    GRAPHICS_LINKER_FLAGS := `$(PKG_CONFIG) --libs gtk+-3.0`
    EXECUTABLE_FILE = praat
  endif

  ifeq ($(PRAAT_AUDIO),none)
    AUDIO_COMPILER_FLAGS :=
    AUDIO_LINKER_FLAGS :=
  else ifeq ($(PRAAT_AUDIO),jack)
    AUDIO_COMPILER_FLAGS := -DJACK
    AUDIO_LINKER_FLAGS := `$(PKG_CONFIG) --libs jack`
  else ifeq ($(PRAAT_AUDIO),alsa)
    AUDIO_COMPILER_FLAGS := -DALSA
    AUDIO_LINKER_FLAGS := -lasound
  else ifeq ($(WE_HAVE_RASPBERRY_PI),1)
    AUDIO_COMPILER_FLAGS := -DJACK
    AUDIO_LINKER_FLAGS := -ljack
  else
    AUDIO_COMPILER_FLAGS := -DALSA -DJACK -DHAVE_PULSEAUDIO
    AUDIO_LINKER_FLAGS := `$(PKG_CONFIG) --libs jack` -lpulse -lasound
  endif

  ifeq ($(WE_HAVE_RASPBERRY_PI),1)
    SYSTEM_COMPILER_FLAGS = -DUNIX -Dlinux -Draspberrypi
  else
    SYSTEM_COMPILER_FLAGS = -DUNIX -Dlinux
  endif

  SHARED_COMPILER_FLAGS := $(SYSTEM_COMPILER_FLAGS) -D_FILE_OFFSET_BITS=64 \
    $(ARCH_COMPILER_FLAGS) \
    $(GRAPHICS_COMPILER_FLAGS) $(AUDIO_COMPILER_FLAGS) \
    -Wreturn-type -Wunused -Wunused-parameter -Wuninitialized -O3 -g1 -pthread
  ifeq ($(PRAAT_FOR_CHROME),1)
    SHARED_COMPILER_FLAGS += -Dchrome
  endif

  CFLAGS := -std=gnu99 $(SHARED_COMPILER_FLAGS) -Werror=missing-prototypes -Werror=implicit

  CXXFLAGS := -std=c++17 $(SHARED_COMPILER_FLAGS) -Wshadow

  ifeq ($(PRAAT_COMPILER),clang)
    CC := clang
    CXX := clang
    LINKER_COMMAND := clang
    CXXFLAGS += -stdlib=libc++
    NON_PRAAT_LIBRARIES := $(GRAPHICS_LINKER_FLAGS) -no-pie -lc++ -lm $(AUDIO_LINKER_FLAGS) -lpthread
  else
    CC := gcc
    CXX := g++
    LINKER_COMMAND := g++
    ifeq ($(WE_HAVE_RASPBERRY_PI),1)
      # The static library /usr/lib/gcc/arm-linux-gnueabihf/8/libstdc++fs.a is needed with GCC 8,
      # and can be removed when using GCC 9
      NON_PRAAT_LIBRARIES := $(GRAPHICS_LINKER_FLAGS) -no-pie -lm $(AUDIO_LINKER_FLAGS) \
        -static-libgcc -static-libstdc++ -lpthread -latomic -ldl -lstdc++fs
    else ifeq ($(PRAAT_FOR_CHROME),1)
      NON_PRAAT_LIBRARIES := $(GRAPHICS_LINKER_FLAGS) -no-pie -lm $(AUDIO_LINKER_FLAGS) \
        -static-libgcc -static-libstdc++ -lpthread -L /usr/lib/x86_64-linux-gnu
    else ifeq ($(PRAAT_GRAPHICS),barrenXXX)
      NON_PRAAT_LIBRARIES := $(GRAPHICS_LINKER_FLAGS) -no-pie -lm $(AUDIO_LINKER_FLAGS) \
        -static-libgcc -static-libstdc++ -lpthread -L /usr/lib/x86_64-linux-gnu
    else ifeq ($(PRAAT_GRAPHICS),noguiXXX)
      NON_PRAAT_LIBRARIES := $(GRAPHICS_LINKER_FLAGS) -no-pie -lm $(AUDIO_LINKER_FLAGS) \
        -static-libgcc -static-libstdc++ -lpthread
    else
      NON_PRAAT_LIBRARIES := $(GRAPHICS_LINKER_FLAGS) -no-pie -lm $(AUDIO_LINKER_FLAGS) -lpthread
    endif
  endif

  AR = ar
  RANLIB = ls
  ICON =
  MAIN_ICON =

  INSTALL = install -p praat /usr/local/bin
endif

export

.PHONY: all clean install

$(info CC      = $(CC))
$(info CXX     = $(CXX))
$(info LINKER_COMMAND    = $(LINKER_COMMAND))
$(info LD      = $(LD))
$(info origin LINKER_COMMAND = $(origin LINKER_COMMAND))
$(info flavor LINKER_COMMAND = $(flavor LINKER_COMMAND))

# Makes the Praat executable in the source directory.
all: all-external all-self
	$(LINKER_COMMAND) -o $(EXECUTABLE_FILE) main/main_Praat.o $(MAIN_ICON) fon/libfon.a \
		artsynth/libartsynth.a FFNet/libFFNet.a \
		gram/libgram.a EEG/libEEG.a \
		LPC/libLPC.a dwtools/libdwtools.a sensors/libsensors.a \
		foned/libfoned.a fon/libfon.a stat/libstat.a \
		dwsys/libdwsys.a sys/libsys.a melder/libmelder.a kar/libkar.a \
		external/espeak/libespeak.a \
		external/portaudio/libportaudio.a \
		external/flac/libflac.a external/lame/liblame.a external/mp3/libmp3.a \
		external/glpk/libglpk.a \
		external/clapack/libclapack.a \
		external/gsl/libgsl.a \
		external/num/libnum.a \
		external/vorbis/libvorbis.a \
		external/opusfile/libopusfile.a \
		external/whispercpp/libwhisper.a \
		$(NON_PRAAT_LIBRARIES)

all-external:
	$(MAKE) -C external/clapack
	$(MAKE) -C external/gsl
	$(MAKE) -C external/glpk
	$(MAKE) -C external/lame
	$(MAKE) -C external/mp3
	$(MAKE) -C external/num
	$(MAKE) -C external/flac
	$(MAKE) -C external/portaudio
	$(MAKE) -C external/espeak
	$(MAKE) -C external/vorbis
	$(MAKE) -C external/opusfile
	$(MAKE) -C external/whispercpp

all-self:
	$(MAKE) -C kar
	$(MAKE) -C melder
	$(MAKE) -C sys
	$(MAKE) -C dwsys
	$(MAKE) -C stat
	$(MAKE) -C fon
	$(MAKE) -C foned
	$(MAKE) -C dwtools
	$(MAKE) -C LPC
	$(MAKE) -C EEG
	$(MAKE) -C sensors
	$(MAKE) -C gram
	$(MAKE) -C FFNet
	$(MAKE) -C artsynth
	$(MAKE) -C main main_Praat.o $(ICON)

clean: clean-external clean-self
	$(RM) praat

clean-external:
	$(MAKE) -C external/clapack clean
	$(MAKE) -C external/gsl clean
	$(MAKE) -C external/glpk clean
	$(MAKE) -C external/lame clean
	$(MAKE) -C external/mp3 clean
	$(MAKE) -C external/num clean
	$(MAKE) -C external/flac clean
	$(MAKE) -C external/portaudio clean
	$(MAKE) -C external/espeak clean
	$(MAKE) -C external/vorbis clean
	$(MAKE) -C external/opusfile clean
	$(MAKE) -C external/whispercpp clean

clean-self:
	$(MAKE) -C kar clean
	$(MAKE) -C melder clean
	$(MAKE) -C sys clean
	$(MAKE) -C dwsys clean
	$(MAKE) -C stat clean
	$(MAKE) -C fon clean
	$(MAKE) -C foned clean
	$(MAKE) -C dwtools clean
	$(MAKE) -C LPC clean
	$(MAKE) -C EEG clean
	$(MAKE) -C sensors clean
	$(MAKE) -C gram clean
	$(MAKE) -C FFNet clean
	$(MAKE) -C artsynth clean
	$(MAKE) -C main clean

install:
	$(INSTALL)
