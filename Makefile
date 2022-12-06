CC = clang++
CFLAGS = -std=c++2b -O3
LFLAGS =
PROG_NAME = a.out

SOURCE_DIR = src
OBJECT_DIR = obj

source  := $(wildcard $(SOURCE_DIR)/*.cpp)
objects  = $(wildcard $(OBJECT_DIR)/*.o)

default: build_all

build_all:
	make compile
	make link

link: $(objects)
	$(CC) $(CFLAGS) $(LFLAGS) $(objects) -o $(PROG_NAME)

compile: $(source)
	$(CC) $(CFLAGS) -c $(source)
	mv *.o $(OBJECT_DIR)

clean: