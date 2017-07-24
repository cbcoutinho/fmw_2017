FFLAGS = -g -O0
FF := gfortran
RM := rm -rf

all: clean main
	./main.exe

main: fd1d_heat_explicit.o mod_types.o
	$(FF) $< -o $@

fd1d_heat_explicit.o: fd1d_heat_explicit.f90 mod_types.o
	$(FF) -c $< -o $@

mod_types.o: mod_types.f90
	$(FF) -c $< -o $@

clean:
	$(RM) *.mod *.o
