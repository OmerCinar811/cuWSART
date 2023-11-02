CC = nvcc
CFLAGS = -g
LFLAGS = 
BIN = 
OBJ = 

default: $(BIN)

run: default
	./bin/$(BIN)

$(BIN):

clean:
	rm -rf *.o $(BIN)