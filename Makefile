NX     := 128
NY     := 384
LX     := 1.
LY     := 3.
DOMAIN := -DNX=$(NX) -DNY=$(NY) -DLX=$(LX) -DLY=$(LY)

CC     := cc
CFLAG  := -std=c99 -Wall -Wextra -Werror -O3 $(ARG_CFLAG) $(DOMAIN)
INC    := -Iinclude
LIB    := -lm
SRCDIR := src
OBJDIR := obj
SRCS   := $(shell find $(SRCDIR) -type f -name "*.c")
OBJS   := $(patsubst %.c,$(OBJDIR)/%.o,$(SRCS))
DEPS   := $(patsubst %.c,$(OBJDIR)/%.d,$(SRCS))
OUTDIR := output
TARGET := a.out

help:
	@echo "all     : create \"$(TARGET)\""
	@echo "clean   : remove \"$(TARGET)\" and object files under \"$(OBJDIR)\""
	@echo "output  : create \"$(OUTDIR)\" to store output"
	@echo "datadel : clean-up \"$(OUTDIR)\""
	@echo "help    : show this message"

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAG) -o $@ $^ $(LIB)

$(OBJDIR)/%.o: %.c
	@if [ ! -e $(dir $@) ]; then \
		mkdir -p $(dir $@); \
	fi
	$(CC) $(CFLAG) -MMD $(INC) -c $< -o $@

clean:
	$(RM) -r $(OBJDIR) $(TARGET)

output:
	@if [ ! -e $(OUTDIR)/log ]; then \
		mkdir -p $(OUTDIR)/log; \
	fi
	@if [ ! -e $(OUTDIR)/save ]; then \
		mkdir -p $(OUTDIR)/save; \
	fi

datadel:
	$(RM) -r $(OUTDIR)/log/*
	$(RM) -r $(OUTDIR)/save/*

-include $(DEPS)

.PHONY : all clean output datadel help

