# Makefile for test use

CC     := cc
CFLAG  := -DTRANSPOSE_TEST -std=c99 -Wall -Wextra -Werror $(ARG_CFLAG)
INC    := -I../../include
LIB    := -lm
SRCS   := test.c main.c
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

