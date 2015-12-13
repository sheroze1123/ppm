CC=icpc

SERIAL_EXECUTABLE=serial
SOURCES=serial.cpp
CFLAGS=-std=c++11 -O3 -Wall -Wextra -Werror
LDFLAGS=-lfftw3
LIB=-L/opt/local/lib -Lfftw/lib
INC=-I/opt/local/include -Ifftw/include
FLAGS=$(LIB) $(INC) $(LDFLAGS) $(CFLAGS)

OMP_EXECUTABLE=ppm_omp
OMP_SOURCES=ppm_omp.cpp
OMP_CFLAGS=$(CFLAGS) -openmp
OMP_LDFLAGS=-lfftw3_omp -lfftw3 -lm
OMP_LIB=-Lfftw/lib
OMP_INC=-Ifftw/include
OMP_FLAGS=$(OMP_LIB) $(OMP_INC) $(OMP_LDFLAGS) $(OMP_CFLAGS)

OBJECTS = marshaller.o

# If you invoke make with no arguments (i.e. `make`), then MARSHAL defaults to
# false and MARSHAL_FLAG defaults to nothing. If you invoke make like this:
# `make MARSHAL=true`, then MARSHAL is overwritten and MARSHAL_FLAG is set to
# `-DMARSHAL`. When the `-DMARSHAL` flag is passed into the compiler, the code
# performs marshalling; otherwise it doesn't.
MARSHAL = false
ifeq ($(MARSHAL), true)
MARSHAL_FLAG = -DMARSHAL
else
MARSHAL_FLAG =
endif

all: $(SOURCES) $(OBJECTS) $(SERIAL_EXECUTABLE) $(OMP_EXECUTABLE)

## Building
$(SERIAL_EXECUTABLE): $(SOURCES) $(OBJECTS)
	$(CC) $(SOURCES) $(OBJECTS) $(FLAGS) $(MARSHAL_FLAG) -o $@

$(OMP_EXECUTABLE): $(OMP_SOURCES) $(OBJECTS)
	$(CC) $(OMP_SOURCES) $(OBJECTS) $(OMP_FLAGS) $(MARSHAL_FLAG) -o $@

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
