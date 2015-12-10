CC=icpc
CFLAGS=-std=c++11 -O3 -lfftw3
LDFLAGS=
SOURCES=serial.cpp
EXECUTABLE=serial
LIB=-L/opt/local/lib -Lfftw/lib
INC=-I/opt/local/include -Ifftw/include

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CC) $(LIB) $(INC) $(LDFLAGS) $(CFLAGS) $(SOURCES) -o $@

clean:
	rm $(EXECUTABLE)

