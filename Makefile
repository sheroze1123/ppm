CC=icpc

SERIAL_EXECUTABLE=serial
SOURCES=serial.cpp
CFLAGS=-std=c++11 -O3 -Wall -Wextra -Werror
LDFLAGS=-lfftw3
LIB=-L/opt/local/lib -Lfftw/lib
INC=-I/opt/local/include -Ifftw/include

OMP_EXECUTABLE=ppm_omp
OMP_SOURCES=ppm_omp.cpp
OMP_CFLAGS=$(CFLAGS) -openmp
OMP_LDFLAGS=-lfftw3_omp -lfftw3 -lm
OMP_LIB=-Lfftw/lib
OMP_INC=-Ifftw/include

OBJECTS = marshaller.o

all: $(SOURCES) $(OBJECTS) $(SERIAL_EXECUTABLE) $(OMP_EXECUTABLE)

## Building
$(SERIAL_EXECUTABLE): $(SOURCES) $(OBJECTS)
	$(CC) $(SOURCES) $(OBJECTS) $(LIB) $(INC) $(LDFLAGS) $(CFLAGS) -o $@

$(OMP_EXECUTABLE): $(OMP_SOURCES) $(OBJECTS)
	$(CC) $(OMP_SOURCES) $(OBJECTS) $(OMP_LIB) $(OMP_INC) $(OMP_LDFLAGS) $(OMP_CFLAGS) -o $@

%.o: %.cpp
	$(CC) -c $< $(LIB) $(INC) $(LDFLAGS) $(CFLAGS)

## Miscellaneous
# watch the status of your qsubbed jobs in reverse chronological order
.PHONY: watch
watch:
	watch -d -n 1 'qstat | tac'

# http://blog.jgc.org/2015/04/the-one-line-you-should-add-to-every.html
print-%:
	@echo $*=$($*)

## Cleaning
clean:
	rm -f *.csv
	rm -f $(OBJECTS)
	rm -f $(SERIAL_EXECUTABLE)
	rm -f $(OMP_EXECUTABLE)
	rm -f *.o[0-9][0-9][0-9][0-9][0-9]
