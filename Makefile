CC=icpc
CFLAGS=-std=c++11 -O3
LDFLAGS=-lfftw3 -lm
SOURCES=serial.cpp
EXECUTABLE=serial
LIB=-L/opt/local/lib -Lfftw/lib
INC=-I/opt/local/include -Ifftw/include

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CC) $(SOURCES) $(LIB) $(INC) $(LDFLAGS) $(CFLAGS) -o $@

clean:
	rm $(EXECUTABLE)

