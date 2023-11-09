CC = nvcc
CFLAGS = -g -Wall -
LFLAGS = 
BIN = 
OBJ = 

default: $(BIN)

run: default
	./bin/$(BIN)

$(BIN):

clean:
	rm -rf *.o $(BIN)