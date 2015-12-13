CC=icpc
CFLAGS=-std=c++11 -O3 -Wall -Wextra -Werror 
LDFLAGS=-lfftw3
SOURCES=serial.cpp
SERIAL_EXECUTABLE=serial
LIB=-L/opt/local/lib -Lfftw/lib
INC=-I/opt/local/include -Ifftw/include
OBJECTS = marshaller.o
OMP_SOURCES=ppm_omp.cpp
OMP_EXECUTABLE=ppm_omp
OMP_CFLAGS=-std=c++11 -O3 -openmp
OMP_LDFLAGS=-lfftw3_omp -lfftw3 -lm
OMP_LIB=-Lfftw/lib
OMP_INC=-Ifftw/include

all: $(SOURCES) $(OBJECTS) $(SERIAL_EXECUTABLE) $(OMP_EXECUTABLE)

$(SERIAL_EXECUTABLE): $(SOURCES) $(OBJECTS)
	$(CC) $(SOURCES) $(OBJECTS) $(LIB) $(INC) $(LDFLAGS) $(CFLAGS) -o $@

$(OMP_EXECUTABLE): $(OMP_SOURCES) $(OBJECTS)
	$(CC) $(OMP_SOURCES) $(OBJECTS) $(OMP_LIB) $(OMP_INC) $(OMP_LDFLAGS) $(OMP_CFLAGS) -o $@

%.o: %.cpp
	$(CC) -c $(LIB) $(INC) $(LDFLAGS) $(CFLAGS) $<

clean:
	rm -f *.csv
	rm -f $(OBJECTS)
	rm -f $(SERIAL_EXECUTABLE)
	rm -f $(OMP_EXECUTABLE)
	rm -f *.o[0-9][0-9][0-9][0-9][0-9]
