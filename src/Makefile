SRC=MD.c control.c
TEST=test/diff-output.c
OBJ=$(SRC:.c=.o)
TEST_O=$(TEST:.c=.o)
CC=pgcc -Minfo
ARCH=px
CFLAGS=-fast -fastsse -Mipa=fast -tp=$(ARCH)

.PHONY: all test clean submit

all: MD.out test
test: diff-output.out

diff-output.out: $(TEST)
	$(CC) -o $@ $(TEST) -lm

MD.out: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

output.dat: MD.out input.dat
	./MD.out

clean:
	rm -f MD.out $(OBJ) $(TEST_O) *.oo diff-output.out
	rm -rf oprofile_data
	rm -f output.dat*

$(OBJ) : coord.h Makefile
