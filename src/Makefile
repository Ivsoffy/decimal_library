UNAME:=$(shell uname -s)
CC = gcc
CFLAGS  = -Wall -Wextra -Werror -std=c11 -c
GCOVFLAGS = -lgcov --coverage

SOURCES = $(wildcard s21_*.c)
TESTS_SOURSES = $(wildcard ./tests/*.check)
TESTS = $(TESTS_SOURSES:.check=.test)
OBJECTS = $(SOURCES:.c=.o)
ifeq ($(UNAME),Linux)
	TEST_FLAGS = -pthread -lcheck_pic -pthread -lrt -lm -lsubunit
endif
ifeq ($(UNAME),Darwin)
	TEST_FLAGS = $(shell pkg-config --cflags --libs check)
endif

all: s21_decimal.a

rebuild: clean all

s21_decimal.a: $(OBJECTS) s21_decimal.h 
	ar rc $@ $(OBJECTS)
	ranlib $@

lib_cov.a: $(SOURCES) s21_decimal.h
	$(CC) $(CFLAGS) -c *.c -g $(GCOVFLAGS)
	ar rc $@ *.o
	ranlib $@

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ -g

%.test: %.check s21_decimal.a
	checkmk $< > $@.c
	$(CC) $@.c -L. -lcheck -lm s21_decimal.a $(GCOVFLAGS) -o $@
	rm -f $@.c
	./$@ 

gcov_report:
	gcc -fprofile-arcs -ftest-coverage tests.c s21_decimal.c -o gcov_report $(TEST_FLAGS)
	./gcov_report
	lcov --no-external --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory out
	open out/index.html

test: clean $(TESTS_SOURSES) s21_decimal.a
	$(CC) -o test tests.c s21_decimal.a -lcheck $(TEST_FLAGS)
	./test

clean:
	rm -rf *.a *.o *.test *test.c .test.c *.gcda *.gcno *.info test.c out test *.input *.output */*.exe *gcov_report

main: s21_decimal.a
	gcc main.c s21_decimal.a -g

clang_format:
	clang-format -i -n -style=Google *.c *.h

leak:
ifeq ($(UNAME),Darwin)
	CK_FORK=no leaks -atExit -- ./test
else
	valgrind --leak-check=full -s --track-origins=yes ./test
endif
