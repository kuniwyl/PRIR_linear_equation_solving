CC = gcc
CFLAGS = -Wall -Wextra

SRCS = main.c gauss_seidel.c gaussian_elimination.c
OBJS = $(SRCS:.c=.o)
HEADERS = gauss_seidel.h gaussian_elimination.h

TARGET = program

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
