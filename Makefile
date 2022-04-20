CC = mpicc
CFLAGS = -Wall -Wextra
LDFLAGS = -lm

TARGET = main
BENCHMARK = main_bench
DBG = main_dbg

SRCS = function.c utils.c $(TARGET).c 
OBJS := $(patsubst %.c,%.o,$(SRCS))
OBJSBENCH := $(patsubst %.c,%_bench.o,$(SRCS))
OBJSDBG := $(patsubst %.c,%_dbg.o,$(SRCS))

all: $(TARGET) $(BENCHMARK) $(DBG)

$(TARGET): $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(BENCHMARK): $(OBJSBENCH)
	$(CC) $^ -o $@ $(LDFLAGS)

$(DBG): $(OBJSDBG)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

%_bench.o: %.c
	$(CC) -DBENCHMARK $(CFLAGS) -c $< -o $@

%_dbg.o: %.c
	$(CC) -DDEBUG $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJS) $(OBJSBENCH) $(OBJSDBG) $(TARGET) $(BENCHMARK) $(DBG) output.txt
