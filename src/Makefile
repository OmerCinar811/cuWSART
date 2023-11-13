CC = nvcc
CFLAGS = -O3 -g -I/usr/include/cuda/include
LFLAGS = -lcudart -ltiff -L/usr/local/cuda/lib64

default: wsart

test:
	./wsart ../Shepp_logan.tif out.tif

wsart_main.o: wsart_main.cu wsart_kernel.cu wsart_support.h
	$(CC) -c -o $@ wsart_main.cu $(CFLAGS)

wsart_support.o: wsart_support.cu wsart_support.h
	$(CC) -c -o $@ wsart_support.cu $(NVCC_FLAGS)

wsart:	wsart_main.o wsart_support.o
	$(CC) wsart_main.o wsart_support.o -o wsart $(LFLAGS)

clean:
	rm -rf *.o wsart