CC      = gcc
CFLAGS  = -std=c99 -O3 -Wall -Wextra -pedantic -D_POSIX_C_SOURCE=200809L
TARGET  = motif_finder
SRC     = motif_finder.c


$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)
