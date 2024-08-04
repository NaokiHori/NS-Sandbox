# Makefile for test use

CC     := cc
CFLAG  := -DDCT_TEST -std=c99 -Wall -Wextra -Werror $(ARG_CFLAG)
INC    := -I../../include
LIB    := -lm
SRCS   := ../memory.c test.c main.c
TARGET := a.out

help:
	@echo "all   : create \"$(TARGET)\""
	@echo "clean : remove \"$(TARGET)\""
	@echo "help  : show this message"

all:
	$(CC) $(CFLAG) $(INC) $(SRCS) -o $(TARGET) $(LIB)

clean:
	$(RM) -r $(TARGET)

.PHONY : all clean help

