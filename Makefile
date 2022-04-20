CC = mpicc
CFLAGS = -Wall -Wextra
LDFLAGS = -lm

TARGET = main
BENCHMARK = main_bench
DBG = main_dbg

SRCSDIR = src

SRCS = $(SRCSDIR)/function.c $(SRCSDIR)/utils.c $(SRCSDIR)/$(TARGET).c 
OBJS := $(patsubst $(SRCSDIR)/%.c,%.o,$(SRCS))
OBJSBENCH := $(patsubst $(SRCSDIR)/%.c,%_bench.o,$(SRCS))
OBJSDBG := $(patsubst $(SRCSDIR)/%.c,%_dbg.o,$(SRCS))

all: $(TARGET) $(BENCHMARK) $(DBG)

$(TARGET): $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(BENCHMARK): $(OBJSBENCH)
	$(CC) $^ -o $@ $(LDFLAGS)

$(DBG): $(OBJSDBG)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: $(SRCSDIR)/%.c
	$(CC) $(CFLAGS) -c $<

%_bench.o: $(SRCSDIR)/%.c
	$(CC) -DBENCHMARK $(CFLAGS) -c $< -o $@

%_dbg.o: $(SRCSDIR)/%.c
	$(CC) -DDEBUG $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJS) $(OBJSBENCH) $(OBJSDBG) $(TARGET) $(BENCHMARK) $(DBG) output.txt
