TARGET_EXEC := build/main

SRCS := src/main.c src/utils.c

LDFLAGS += -lm

.PHONY: debug
debug:
	$(CC) -Og -fsanitize=address -o $(TARGET_EXEC) $(SRCS) $(CFLAGS) $(LDFLAGS)

.PHONY: release
release:
	$(CC) -O3 -o $(TARGET_EXEC) $(SRCS) $(CFLAGS) $(LDFLAGS)
	strip $(TARGET_EXEC)

.PHONY: clean
clean:
	rm main data.txt
