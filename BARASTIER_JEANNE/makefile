CC = gcc
TARGET = pmns
SRC = main.c
HEADERS = param.h

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CC) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)

.PHONY: all clean
