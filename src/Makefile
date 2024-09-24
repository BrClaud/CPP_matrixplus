CC = g++
CFLAGS =-Wall -Werror -Wextra  -g
COVFLAGS = -fprofile-arcs  -lcheck -ftest-coverage

all: clean s21_matrix_oop.a test

s21_matrix_oop.a:
	$(CC) $(CFLAGS) -c $(COVFLAGS) s21_matrix_oop.cpp
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

gcov_report: test
	gcov -f *.gcda
	lcov -t "test" -o test.info -c -d .
	genhtml -o report test.info 
	open report/index.html

test: s21_matrix_oop.a
	$(CC) $(CFLAGS) -c test/tests*.cpp
	$(CC) $(CFLAGS) --coverage -o test.out tests*.o -lgtest -lgtest_main -L. s21_matrix_oop.a
	./test.out

clean:
	rm -rf *.o *.a *.gcov *.gcno *.gcda ./report *.info *.out