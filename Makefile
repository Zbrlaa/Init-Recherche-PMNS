
CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c11
SRC = main.c
TARGET = main

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRC) param.h
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)