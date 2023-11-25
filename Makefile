TARGET_EXEC := main

SRCS := main.c utils.c

LDFLAGS += -lm

.PHONY: debug
debug:
	$(CC) $(CFLAGS) $(LDFLAGS) -Og -fsanitize=address -o $(TARGET_EXEC) $(SRCS)

.PHONY: release
release:
	$(CC) $(CFLAGS) $(LDFLAGS) -O2 -o $(TARGET_EXEC) $(SRCS)

.PHONY: clean
clean:
	rm main
