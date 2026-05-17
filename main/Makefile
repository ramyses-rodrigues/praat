# Makefile of the library "main"
# Paul Boersma, 16 May 2026

CPPFLAGS = -I ../kar -I ../melder -I ../sys -I ../fon

.PHONY: clean

main_Praat.o: ../sys/praat.h ./main_Praat.h ../sys/Graphics.h

clean:
	$(RM) main_*.o
	$(RM) praat_win.o

# The following is only for Windows.

praat_win.o: praat_win.rc
	$(WINDRES) -o praat_win.o praat_win.rc
