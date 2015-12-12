CC=icpc
CFLAGS=-std=c++11 -O3 -Wall -Wextra -Werror
LDFLAGS=-lfftw3
SOURCES=serial.cpp
EXECUTABLE=serial
LIB=-L/opt/local/lib -Lfftw/lib
INC=-I/opt/local/include -Ifftw/include
OBJECTS = marshaller.o server.o

all: $(SOURCES) $(OBJECTS) $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES) $(OBJECTS)
	$(CC) $(SOURCES) $(OBJECTS) $(LIB) $(INC) $(LDFLAGS) $(CFLAGS) -o $@

%.o: %.cpp
	$(CC) -c $(LIB) $(INC) $(LDFLAGS) $(CFLAGS) $<

clean:
	rm -f *.csv
	rm -f $(OBJECTS)
	rm -f $(EXECUTABLE)
