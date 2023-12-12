TARGET_EXEC := build/main

SRCS := src/*.c

COMPLIE_FLAGS += -lm -fopenmp

.PHONY: debug
debug:
	$(CC) -Og -fsanitize=address -o $(TARGET_EXEC) $(SRCS) $(CFLAGS) $(LDFLAGS) $(COMPLIE_FLAGS)

.PHONY: release
release:
	$(CC) -O3 -o $(TARGET_EXEC) $(SRCS) $(CFLAGS) $(LDFLAGS) $(COMPLIE_FLAGS)
	strip $(TARGET_EXEC)

.PHONY: clean
clean:
	rm main data.txt
