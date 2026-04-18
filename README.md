# Praat: doing phonetics by computer

Welcome to Praat! Praat is a speech analysis tool used for doing phonetics by computer.
Praat can analyse, synthesize, and manipulate speech, and create high-quality pictures for your publications.
Praat was created by Paul Boersma and David Weenink of the Institute of Phonetics Sciences of the University of Amsterdam.

Some of Praat’s most prominent features are:

#### Speech analysis

Praat allows you to analyze different aspects of speech including pitch, formant, intensity, and voice quality.
You have access to spectrograms (a visual representation of sound changing over time)
and cochleagrams (a specific type of spectrogram more closely resembling how the inner ear receives sound).

#### Speech synthesis

Praat allows you to generate speech from a pitch curve and filters that you create (acoustic synthesis),
or from muscle activities (articulatory synthesis).

#### Speech manipulation

Praat gives you the ability to modify existing speech utterances. You can alter pitch, intensity, and duration of speech.

#### Speech labelling

Praat allows you to custom-label your samples using the IPA (International Phonetics Alphabet),
and annotate your sound segments based on the particular variables you are seeking to analyze.
Multi-language text-to-speech facilities allow you to segment the sound into words and phonemes.

#### Grammar models

With Praat, you can try out Optimality-Theoretic and Harmonic-Grammar learning,
as well as several kinds of neural-network models.

#### Statistical analysis

Praat allows you to perform several statistical techniques, among which
multidimensional scaling, principal component analysis, and discriminant analysis.

For more information, consult the extensive manual in Praat (under Help),
and the website [praat.org](https://praat.org), which has Praat tutorials in several languages.

## 1. Binary executables

While the [Praat website](https://praat.org) contains the latest executable for all platforms that we support
(or used to support), the [releases on GitHub](https://github.com/praat/praat.github.io/releases) contain many older executables as well.

The meaning of the names of binary files available on GitHub is as follows (editions that currently receive updates are in bold):

### 1.1. Windows binaries
- **`praatXXXX_win-arm64.zip`: zipped executable for ARM64 Windows (11 and higher)**
- **`praatXXXX_win-intel64.zip`: zipped executable for Intel64/AMD64 Windows (7 and higher)**
- **`praatXXXX_win-intel32.zip`: zipped executable for Intel32 Windows (7 and higher)**
- `praatXXXX_win64.zip`: zipped executable for Intel64/AMD64 Windows (XP and higher, or 7 and higher)
- `praatXXXX_win32.zip`: zipped executable for Intel32 Windows (XP and higher, or 7 and higher)
- `praatconXXXX_win64.zip`: zipped executable for Intel64/AMD64 Windows, console edition
- `praatconXXXX_win32.zip`: zipped executable for Intel32 Windows, console edition
- `praatconXXXX_win32sit.exe`: self-extracting StuffIt archive with executable for Intel32 Windows, console edition
- `praatXXXX_win98.zip`: zipped executable for Windows 98
- `praatXXXX_win98sit.exe`: self-extracting StuffIt archive with executable for Windows 98

### 1.2. Mac binaries
- **`praatXXXX_mac.dmg`: disk image with universal executable for (64-bit) Intel and Apple Silicon Macs (Cocoa)**
- **`praatXXXX_xcodeproj.zip`: zipped Xcode project file for the universal (64-bit) edition (Cocoa)**
- `praatXXXX_mac64.dmg`: disk image with executable for 64-bit Intel Macs (Cocoa)
- `praatXXXX_xcodeproj64.zip`: zipped Xcode project file for the 64-bit edition (Cocoa)
- `praatXXXX_mac32.dmg`: disk image with executable for 32-bit Intel Macs (Carbon)
- `praatXXXX_xcodeproj32.zip`: zipped Xcode project file for the 32-bit edition (Carbon)
- `praatXXXX_macU.dmg`: disk image with universal executable for (32-bit) PPC and Intel Macs (Carbon)
- `praatXXXX_macU.sit`: StuffIt archive with universal executable for (32-bit) PPC and Intel Macs (Carbon)
- `praatXXXX_macU.zip`: zipped universal executable for (32-bit) PPC and Intel Macs (Carbon)
- `praatXXXX_macX.zip`: zipped executable for MacOS X (PPC)
- `praatXXXX_mac9.sit`: StuffIt archive with executable for MacOS 9
- `praatXXXX_mac9.zip`: zipped executable for MacOS 9
- `praatXXXX_mac7.sit`: StuffIt archive with executable for MacOS 7

### 1.3. Linux binaries
- **`praatXXXX_linux-s390x-barren.tar.gz`: gzipped tarred executable for s390x Linux, without GUI, sound and graphics**
- **`praatXXXX_linux-s390x.tar.gz`: gzipped tarred executable for s390x Linux (GTK 3)**
- **`praatXXXX_linux-arm64-barren.tar.gz`: gzipped tarred executable for ARM64 Linux (Ubuntu, Debian...), without GUI, sound and graphics**
- **`praatXXXX_linux-arm64.tar.gz`: gzipped tarred executable for ARM64 Linux (Ubuntu, Debian...) (GTK 3)**
- **`praatXXXX_linux-intel64-barren.tar.gz`: gzipped tarred executable for Intel64/AMD64 Linux (Ubuntu, Debian...), without GUI, sound and graphics**
- **`praatXXXX_linux-intel64.tar.gz`: gzipped tarred executable for Intel64/AMD64 Linux (Ubuntu, Debian...) (GTK 3)**
- `praatXXXX_linux-arm64-nogui.tar.gz`: gzipped tarred executable for ARM64 Linux, without GUI and sound but with graphics (Cairo and Pango)
- `praatXXXX_linux-intel64-nogui.tar.gz`: gzipped tarred executable for Intel64/AMD64 Linux, without GUI and sound but with graphics (Cairo and Pango)
- `praatXXXX_linux64barren.tar.gz`: gzipped tarred executable for Intel64/AMD64 Linux, without GUI, sound and graphics
- `praatXXXX_linux64nogui.tar.gz`: gzipped tarred executable for Intel64/AMD64 Linux, without GUI and sound but with graphics (Cairo and Pango)
- `praatXXXX_linux64.tar.gz`: gzipped tarred executable for Intel64/AMD64 Linux (GTK 2 or 3)
- `praatXXXX_linux32.tar.gz`: gzipped tarred executable for Intel32 Linux (GTK 2)
- `praatXXXX_linux_motif64.tar.gz`: gzipped tarred executable for Intel64/AMD64 Linux (Motif)
- `praatXXXX_linux_motif32.tar.gz`: gzipped tarred executable for Intel32 Linux (Motif)

### 1.4. Chromebook binaries
- **`praatXXXX_chrome-arm64.tar.gz`: gzipped tarred executable for Linux on ARM64 Chromebooks (GTK 3)**
- **`praatXXXX_chrome-intel64.tar.gz`: gzipped tarred executable for Intel64/AMD64 Linux on Intel64/AMD64 Chromebooks (GTK 3)**
- `praatXXXX_chrome64.tar.gz`: gzipped tarred executable for 64-bit Linux on Intel64/AMD64 Chromebooks (GTK 2 or 3)

### 1.5. Raspberry Pi binaries
- **`praatXXXX_rpi-armv7.tar.gz`: gzipped tarred executable for (32-bit) ARMv7 Linux on the Raspberry Pi 4B (GTK 3)**
- `praatXXXX_rpi_armv7.tar.gz`: gzipped tarred executable for (32-bit) ARMv7 Linux on the Raspberry Pi 4B (GTK 2 or 3)

### 1.6. Other Unix binaries (all obsolete)
- `praatXXXX_solaris.tar.gz`: gzipped tarred executable for Sun Solaris
- `praatXXXX_sgi.tar.gz`: gzipped tarred executable for Silicon Graphics Iris
- `praatXXXX_hpux.tar.gz`: gzipped tarred executable for HP-UX (Hewlett-Packard Unix)

## 2. Compiling the source code

You need the Praat source code only in the following cases:

1. you want to extend Praat’s functionality by adding C or C++ code to it; or
2. you want to understand or reuse Praat’s source code; or
3. you want to compile Praat for a computer for which we do not provide binary executables,
e.g. Linux for some non-Intel computers, FreeBSD, HP-UX, SGI, or SPARC Solaris.

Before trying to dive into Praat’s source code, you should be familiar with the working of the Praat program
and with writing Praat scripts. The Praat program can be downloaded from
https://praat.org or https://www.fon.hum.uva.nl/praat.

### 2.1. License

Most of the source code of Praat is distributed on GitHub under the General Public License,
[version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html) or later,
or [version 3](https://praat.org/manual/General_Public_License__version_3.html) or later.
However, as Praat includes software written by others,
the whole of Praat is distributed under the General Public License,
[version 3](https://praat.org/manual/General_Public_License__version_3.html) or later.
See [Acknowledgments](https://praat.org/manual/Acknowledgments.html) for details on the licenses
of software libraries by others that are included in Praat.
Of course, any improvements in the Praat source code are welcomed by the authors.

### 2.2. Downloading the archive

To download the latest source code of Praat from GitHub,
click on the *zip* or *tar.gz* archive at the latest release,
or fork ("clone") the praat/praat repository at any later change.

### 2.3. Steps to take if you want to extend Praat

First make sure that the source code can be compiled as is.
Then add your own buttons by editing `main/main_Praat.cpp` or `fon/praat_Fon.cpp`.
Consult the manual page on [Programming](https://praat.org/manual/Programming_with_Praat.html).

### 2.4. The programming language

Most of the source code is written in C++, but some parts are written in C.
The code requires that your compiler supports C99 and C++17.

## 3. Developing Praat for one platform

Developing Praat means two things: *building* the Praat executable from the Praat
source code, and *testing* the correctness of the Praat executable.

**Building** is largely automated:

- On the Mac, you use an existing Xcode project (these files are included in the releases).
- On other platforms (Windows and Unixes), you use existing makefiles
  (these files are included in the source tree).

**Testing** on a platform can be done by starting up Praat on that platform,
and then go through to types of tests. Basic GUI functionality is tested as follows:

1. record a sound (New -> `Record mono Sound...`)
2. open the sound (`View & Edit`)
3. select a part of the sound (drag the mouse across the waveform)
4. play the sound (click on the rectangle below or above the selection)

The integrity of Praats’s algorithms (e.g. signal processing)
and of the Praat scripting language is tested as follows:

1. open the script `test/runAlltests.praat` (Praat -> `Open Praat script...`)
2. run the script (Run -> `Run`)
3. after 2 to 10 minutes, the Info window should contain a big “OK” graph
4. go through steps 1 through 3 for `dwtest/runAllTests.praat`
5. if you feel adventurous, try some tests in the folder `test/manually`

### 3.1. Developing Praat for Windows

On Windows, Praat is **built** through the makefiles provided in Praat’s source tree.

One could use Cygwin or MSYS2. As we like to provide not only an Intel64/AMD64 and Intel32 edition,
but an ARM64 edition as well, and Cygwin has no toolchains for ARM64, we work with MSYS2 instead.

After installing MSYS2, we see that a `mingw64` toolchain (for Praat’s Intel64/AMD64 edition)
and a `mingw32` toolchain (for Praat’s Intel32 edition) are already available.
Make sure you have installed at least `make`, `gcc`, `g++` and `pkg-config` to make those work.
To also install a `clangarm64` toolchain (for Praat’s ARM64 edition),
run `clangarm64.exe` to get a `clangarm64` shell. In that shell, run `pacman -Suy` to update and
`pacman -S mingw-w64-clang-aarch64-clang` to install the build tools package.
In the same way you can create a `clang64` toolchain and a `clang32` toolchain
(`pacman -S mingw-w64-clang-x86_64-clang` and `pacman -S mingw-w64-i686-clang`),
which are good alternatives to `mingw64` and `mingw32`;
in 2025, MSYS2 ended support for Intel32, so you may have to use `mingw32` for that,
unless you have an old Intel32 toolchain lying around.

Ramyses: para instalar o toolchain mingw64:
1. Instalar a extensão C/C++ (Microsoft);
1-a. instalar MSYS2 do site (https://www.msys2.org/#installation)
2. instalar pacote de compilação UCRT64, com o comando "pacman -S mingw-w64-ucrt-x86_64-gcc" (antes era "$ pacman -S mingw-w64-x86_64-gcc")
3. instalar make: "pacman -S mingw-w64-ucrt-x86_64-make" (antes era "pacman -S make")
4. instalar gdb: "pacman -S mingw-w64-ucrt-x86_64-gdb" (antes era "$ pacman -S mingw-w64-x86_64-gdb")
5. atualizar Msys2: "pacman -Suy"
6. configurar a extensão Makefile, inserindo o caminho do arquivo Makefile na pasta raiz do Praat;
7. configurar a extensão Makefile, inserindo o caminho do executável MAKE : Apontar para o arquivo "C:\msys64\ucrt64\bin\mingw32-make.exe"
8. Configurar o Path do windows com os caminhos: "\msys2\ucrt64\bin" (onde estão os executáveis do compilador - g++.exe, gdb.exe e gcc.exe) e "\msys2\usr\bin" (onde estão os arquivos do coreutils (touch.exe) do msys2).
9. O arquivo makefile.defs já está editado para usar com o msys2 (File: makefile.defs.msys-mingw64)

Usando a biblioteca UCRT64, os arquivos baixados pelo Pacman serão armazenados na pasta  ou "\msys2\mingw64\bin".

Ramyses: adicionalmente, para usar o CMAKE:
1. instalar a extensão CMake Tools (Microsoft) no vscode;
2. Instalar o pacote cmake e ninja, pelo comando "pacman -S --needed mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja";
3. 
4. 

Move the Praat sources folders somewhere in your `/home/yourname` tree,
perhaps even in three places, e.g. as `/home/yourname/praats-arm64`,
`/home/yourname/praats-intel64` and `/home/yourname/praats-intel32`;
the folders `fon` and `sys` should be visible within each of these folders.

If you now want to build Praat’s ARM64 edition, start the shell `clangarm64` and type

    cd ~/praats-arm64
    cp makefiles/makefile.defs.msys-clang ./makefile.defs
    make -j12

If you want to build Praat’s Intel64/AMD64 edition, start the shell `clang64` and type

    cd ~/praats-intel64
    cp makefiles/makefile.defs.msys-clang ./makefile.defs
    make -j12

or start the shell `mingw64` and type

    cd ~/praats-intel64
    cp makefiles/makefile.defs.msys-mingw64 ./makefile.defs
    make -j12

If you want to build Praat’s Intel32 edition, start the shell `clang32` and type

    cd ~/praats-intel32
    cp makefiles/makefile.defs.msys-clang ./makefile.defs
    make -j12

or start the shell `mingw32` and type

    cd ~/praats-intel32
    cp makefiles/makefile.defs.msys-mingw32 ./makefile.defs
    make -j12

(With Cygwin, you would install the Devel package `mingw64-x86_64-gcc-g++`
for Praat’s Intel64/AMD64 edition and `mingw64-i686-gcc-g++` for Praat’s Intel32 edition,
plus perhaps `make` and `pkg-config` if you dont’t have those yet.)

**Code-signing.** From version 6.4.25 on, we have signed the three Praat executables
with an “open-source code-signing certificate” (by Certum)
under the name “Paulus Boersma” (the Dutch-legal name of one of the authors).
This is designed to make it easier for Praat to pass the SmartScreen checks
on Windows 11. Early testing shows that the signature is seen by SmartScreen,
but that SmartScreen can still block Praat, so that users still have to
click “Run Anyway” (or “Unblock” under Properties).
It seems that the card reader for this certificate cannot be used yet for code-signing on ARM64 Windows
or on MacOS Sonoma/Sequoia, so for the moment we have to fall back on an obsolete Intel64/AMD64 Windows 10 machine
and are looking for a more robust solution.

**Testing** on multiple platform versions can be done with virtual machines
for Windows 7 (64-bit), Windows 8.1 (64-bit), 64-bit Windows 10 (1507, 1803, 22H2) and Windows 11,
for instance on an Intel64 Mac with Parallels Desktop.
On an ARM64 Mac with Parallels Desktop, you can test only on Windows 11.

### 3.2. Compiling for Macintosh

To **build** Praat on the Mac, extract the *praatXXXX_xcodeproj.zip* file
from [Praat’s latest release](https://github.com/praat/praat.github.io/releases)
into the folder that contains `sys`, `fon`, `dwtools` and so on (e.g. `~/Dropbox/Praats/src`).
Then open the project `praat.xcodeproj` in Xcode 16.3 (or later),
and edit the Intermediate and Product build paths to something that suits you
(Xcode -> Settings... -> Locations -> Derived Data -> Advanced... -> Custom -> Absolute,
then type something after Products, e.g. `~/Dropbox/Praats/bin/macos`,
as well as something after Intermediates, e.g. `~/builds/mac_intermediates`, then click Done).
After this preliminary work, choose Build or Run for the target `praat_mac`.
You can compile with the 14.2 SDK, which will work as far back as macOS 10.11 El Capitan,
which is our deployment target, and will look good even on macOS 14 Sonoma.

If you get an error message like “Code Signing Identity xxx does not match any valid, non-expired,
code-signing certificate in your keychain”, then select the target `praat_mac`, go to Info → Build,
and switch “Code Signing Identity” to “Don’t Code Sign”,
or sign with your own certificate if you have one as a registered Apple developer.

If you get lots of errors saying “Expected unqualified-id” or “Unknown type name NSString”,
then you may have to switch the Type of some .cpp file from “C++ Source” to “Objective-C++ Source”
(under “Identity and Type” in the righthand sidebar).

If you want to build Praat as a library instead of as an executable,
try the target `praat_mac_a` (static) or `praat_mac_so` (dynamic).

**Notarization.** If you want others to be able to use your Mac app,
you will probably have to not only *sign* the executable, but also *notarize* it. To this end,
do Xcode (version 16) -> Product -> Archive -> Distribute App -> Developer ID -> Upload ->
Automatically manage signing -> Upload -> ...wait... (“Package Approved”) ...wait...
(“Ready to distribute”) -> Export Notarized App). If your Praat.app was built into
`~/Dropbox/Praats/bin/macos/Configuration64`, then you can save the notarized
`Praat.app` in `~/Dropbox/Praats/bin/macos`, then drag it in the Finder to
`~/Dropbox/Praats/bin/macos/Configuration64`, overwriting the non-notarized
Praat.app that was already there. If on the way you receive an error
“App Store Connect Operation Error -- You must first sign the relevant contracts online”,
or “Couldn’t communicate with a helper application“,
you will have to log in to `developer.apple.com` and do Review Agreement -> Agree;
you may then also have to either wait for a couple of minutes,
and or go (or log in) to App Store Connect, then Agreements, Tax, and Banking
-> Paid Apps -> View or Terms (even if you have no paid apps).

**Testing** on multiple Intel64/AMD64 platform versions can be done on older Intel64 Macs,
using virtual machines with Parallels Desktop. For instance, a 2013 Macbook Pro can handle
OS X 10.11 El Capitan, 10.12 Sierra, 10.13 High Sierra, macOS 10.14 Mojave, 10.15 Catalina,
and macOS 11 Big Sur, while a 2018 Macbook Pro can handle macOS 10.14 Mojave, 10.15 Catalina,
macOS 11 Big Sur, macOS 12 Monterey, and macOS 13 Ventura (and macOS 14 Sonoma natively).
Testing on multiple ARM64 platform versions can be done on an older ARM64 Mac,
using virtual machines with Parallels Desktop. For instance, a 2020 Mac Mini could handle
macOS 11 Big Sur, macOS 12 Monterey, and macOS 13 Ventura (and macOS 14 Sonoma natively),
while a 2023 Macbook Pro can do macOS 14 Sonoma or macOS 15 Sequoia natively.

### 3.3. Compiling on Linux and other Unixes

To set up the system libraries required for **building** with the Clang or GCC compiler,
install the necessary build tools as well as some graphics and sound packages:

    sudo apt install make rsync pkg-config
    # either:
        sudo apt install clang libc++-dev libc++abi-dev
    # or:
        sudo apt install gcc g++
    sudo apt install libgtk-3-dev
    sudo apt install libasound2-dev
    sudo apt install libpulse-dev
    sudo apt install libjack-dev

To set up your source tree for Linux, go to Praat's sources directory (where the folders `fon` and `sys` are)
and type one of the four following commands:

    # on Ubuntu command line (Intel64/AMD64 or ARM64 processor)
    # either:
        cp makefiles/makefile.defs.linux.pulse-clang ./makefile.defs
    # or:
        cp makefiles/makefile.defs.linux.pulse-gcc ./makefile.defs

    # on Ubuntu command line (s390x processor)
    cp makefiles/makefile.defs.linux.s390x.pulse ./makefile.defs

    # on Chromebook command line
    cp makefiles/makefile.defs.chrome64 ./makefile.defs

    # on Raspberry Pi command line
    cp makefiles/makefile.defs.linux.rpi ./makefile.defs
    
    # on FreeBSD command line
    cp makefiles/makefile.defs.freebsd.alsa ./makefile.defs

To build the Praat executable, type `make -j15` or so.
If your Unix isn’t Linux, you may have to edit the library names in the makefile
(you may need pthread, gtk-3, gdk-3, atk-1.0, pangoft2-1.0, gdk_pixbuf-2.0, m, pangocairo-1.0,
cairo-gobject, cairo, gio-2.0, pango-1.0, freetype, fontconfig, gobject-2.0, gmodule-2.0, 
gthread-2.0, rt, glib-2.0, asound, jack).

When compiling Praat on an external supercomputer or so, you will not have sound.
If you do have `libgtk-3-dev` (and its dependencies), do

    # on Ubuntu command line (Intel64/AMD64 or ARM64 processor)
    cp makefiles/makefile.defs.linux.silent ./makefile.defs

Then type `make -j12` or so to build the program. If your Unix isn’t Linux,
you may have to edit the library names in the makefile (you may need pthread, gtk-3, gdk-3, atk-1.0,
pangoft2-1.0, gdk_pixbuf-2.0, m, pangocairo-1.0, cairo-gobject, cairo, gio-2.0, pango-1.0, 
freetype, fontconfig, gobject-2.0, gmodule-2.0, gthread-2.0, rt, glib-2.0).

When compiling Praat for use as a server for commands from your web pages,
you may not need sound, a GUI, amd graphics. In that case, do

    # on Ubuntu command line (Intel64/AMD64 or ARM64 processor)
    # either:
        cp makefiles/makefile.defs.linux.barren-clang ./makefile.defs
    # or:
        cp makefiles/makefile.defs.linux.barren-gcc ./makefile.defs

    # on Ubuntu command line (s390x processor)
    cp makefiles/makefile.defs.linux.s390x.barren ./makefile.defs

which creates the executable `praat_barren`. Then type `make` or `make -j15` to build the program.
If your Unix isn’t Linux, you may have to edit the library names in the makefile.

The above works exactly the same for Intel64/AMD64 and ARM64 processors, with the same makefiles.

**Testing** on multiple platform versions can be done with virtual machines
for e.g. Ubuntu 20.04, Ubuntu 22.04, Fedora 35, Fedora 37, Mint 20.2,
Debian GNU Linux 10.10, CentOS 8.4, and CentOS Stream 9, 
for instance on an Intel64 Mac with Parallels Desktop.
On an ARM64 Mac, we test with virtual machines for Ubuntu 22.04, Fedora 38,
and Debian GNU Linux 12 ARM64.
See [HOW_TO_BUILD_ONE.md](HOW_TO_BUILD_ONE.md).

## 4. Developing Praat on all platforms simultaneously

See [HOW_TO_BUILD_ALL.md](HOW_TO_BUILD_ALL.md).
