FFLAGS = -g -O0
FF := gfortran
RM := rm -rf

objects = mod_io.o mod_types.o mod_rhs.o mod_cfl.o mod_solver.o

all: clean main2
	./main2.exe

plot: all
	python main.py

main2: fd1d_heat_explicit.o $(objects)
	$(FF) $^ -o $@

fd1d_heat_explicit.o: fd1d_heat_explicit.f90 $(objects)
	$(FF) -c $< -o $@

mod_solver.o: mod_solver.f90 mod_types.o mod_rhs.o
	$(FF) -c $< -o $@

mod_types.o: mod_types.f90
	$(FF) -c $< -o $@

mod_io.o: mod_io.f90 mod_types.o mod_rhs.o
	$(FF) -c $< -o $@

mod_cfl.o: mod_cfl.f90 mod_types.o
	$(FF) -c $< -o $@

mod_rhs.o: mod_rhs.f90 mod_types.o
	$(FF) -c $< -o $@

clean:
	$(RM) *.mod *.o *.obj *.exe *.txt
