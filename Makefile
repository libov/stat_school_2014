# Vladyslav Libov, vladyslav.libov@desy.de

FLAGS=`root-config --glibs --cflags` -g

FLAGS_OBJ=`root-config --cflags` -Wall -g

OBJ := tmp
SRC := src
INC := inc

vpath %.o   $(OBJ)
vpath %.cxx $(SRC)
vpath %.C   $(SRC)
vpath %.h   $(INC)

all: mass_peak_fit generate_data

mass_peak_fit.o: mass_peak_fit.cxx
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

generate_data.o: generate_data.cxx
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

TMassPeakFit.o: TMassPeakFit.cxx TMassPeakFit.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

CustomFit.o: CustomFit.cxx TMassPeakFit.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

Fit.o: Fit.cxx TMassPeakFit.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

minimization_functions.o: minimization_functions.cxx TMassPeakFit.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

fit_functions.o: fit_functions.cxx fit_functions.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

mass_peak_fit: mass_peak_fit.o TMassPeakFit.o CustomFit.o fit_functions.o Fit.o minimization_functions.o
	gcc -o mass_peak_fit $(OBJ)/mass_peak_fit.o $(OBJ)/TMassPeakFit.o $(OBJ)/CustomFit.o $(OBJ)/Fit.o $(OBJ)/fit_functions.o $(OBJ)/minimization_functions.o $(FLAGS) -lMinuit

generate_data: generate_data.o
	gcc -o generate_data $(OBJ)/generate_data.o $(FLAGS)

.PHONY: clean
clean:
	-rm -rf $(OBJ)/*.o
	-rm -rf mass_peak_fit
	-rm -rf generate_data
