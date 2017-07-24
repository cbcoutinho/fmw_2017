FFLAGS = -g -O0
FF := gfortran
RM := rm -rf

all: clean main2
	./main2.exe

main2: fd1d_heat_explicit.o mod_types.o mod_helper.o
	$(FF) $^ -o $@

fd1d_heat_explicit.o: fd1d_heat_explicit.f90 mod_types.o mod_helper.o
	$(FF) -c $< -o $@

mod_types.o: mod_types.f90
	$(FF) -c $< -o $@

mod_helper.o: mod_helper.f90 mod_types.o
	$(FF) -c $< -o $@

clean:
	$(RM) *.mod *.o *.exe *.txt
