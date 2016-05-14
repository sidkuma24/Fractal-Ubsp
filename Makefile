SRCS    = mars_enc.c image_io.c index_func.c coding_func.c miscell.c \
          split_func.c nn_search.c mars_dec.c 

OBJ_ENC = mars_enc.o image_io.o index_func.o coding_func.o miscell.o \
          split_func.o nn_search.o 

OBJ_DEC = mars_dec.o miscell.o image_io.o

CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 

all:  encmars decmars

encmars: $(OBJ_ENC)
	$(CC) -o fcom  $(OBJ_ENC) $(LIBS)

decmars: $(OBJ_DEC)
	$(CC) -o fdcom  $(OBJ_DEC) $(LIBS)


mars_enc.o    :def.h globals.h nn_search.h prot.h triangle.h
index_func.o  :def.h globals.h nn_search.h prot.h
coding_func.o :def.h globals.h nn_search.h prot.h
miscell.o     :def.h globals.h nn_search.h prot.h
split_func.o  :def.h globals.h nn_search.h prot.h
image_io.o    :def.h globals.h nn_search.h prot.h
mars_dec.o    :def.h globals.h prot.h triangle.h
nn_search.o   :nn_search.h 

clean:
	rm -f *.o


