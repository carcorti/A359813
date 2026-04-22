# Makefile for A359813 — campaign code for OEIS A359813.
#
# Build targets:
#   debug    : -O0 -g3 -DDEBUG, assert attive
#   release  : -O3 -march=native -DNDEBUG (default)
#   profile  : -O3 -g -fno-omit-frame-pointer (per perf/gprof)
#   sanitize : AddressSanitizer + UBSan, binario separato A359813_san

CC      := gcc
TARGET  := A359813
SRC     := A359813.c
BASE    := -std=c17 -fopenmp -Wall -Wextra -Wpedantic

CFLAGS_DEBUG   := $(BASE) -O0 -g3 -DDEBUG
CFLAGS_RELEASE := $(BASE) -O3 -march=native -DNDEBUG
CFLAGS_PROFILE := $(BASE) -O3 -g -fno-omit-frame-pointer -march=native -DNDEBUG
CFLAGS_SAN     := $(BASE) -O1 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer
LDFLAGS_SAN    := -fsanitize=address,undefined

all: release

debug:
	$(CC) $(CFLAGS_DEBUG) $(SRC) -o $(TARGET)

release:
	$(CC) $(CFLAGS_RELEASE) $(SRC) -o $(TARGET)

profile:
	$(CC) $(CFLAGS_PROFILE) $(SRC) -o $(TARGET)

sanitize:
	$(CC) $(CFLAGS_SAN) $(SRC) -o $(TARGET)_san $(LDFLAGS_SAN)

clean:
	rm -f $(TARGET) $(TARGET)_san *.o A359813_state.bin A359813_state.bin.tmp

.PHONY: all debug release profile sanitize clean
