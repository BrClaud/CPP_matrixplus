#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(Overload, Sum) {
  S21Matrix source1(3, 3), source2(3, 3), expected(3, 3);
  source1(0, 0) = 352.412;
  source1(1, 0) = 1;
  source1(2, 0) = 1;
  source1(0, 1) = 1234.8234;
  source1(1, 1) = 0.00001;
  source1(2, 1) = 352.412;
  source1(0, 2) = 1234.8234;
  source1(1, 2) = 352.412;
  source1(2, 2) = 9;

  source2(0, 0) = 352.412;
  source2(1, 0) = 1;
  source2(2, 0) = 0;
  source2(0, 1) = -1234.8234;
  source2(1, 1) = 0.00001;
  source2(2, 1) = 123456789;
  source2(0, 2) = -0.00001;
  source2(1, 2) = 0.412;
  source2(2, 2) = 9;

  S21Matrix result = source1 + source2;

  expected(0, 0) = 352.412 + 352.412;
  expected(1, 0) = 1 + 1;
  expected(2, 0) = 1;
  expected(0, 1) = -1234.8234 + 1234.8234;
  expected(1, 1) = 0.00001 + 0.00001;
  expected(2, 1) = 352.412 + 123456789;
  expected(0, 2) = -0.00001 + 1234.8234;
  expected(1, 2) = 352.412 + 0.412;
  expected(2, 2) = 9 + 9;
  ASSERT_EQ(result.EqMatrix(expected), 1);
}

TEST(Overload, Sub) {
  S21Matrix source1(3, 3), source2(3, 3), expected(3, 3);
  source1(0, 0) = 352.412;
  source1(1, 0) = 1;
  source1(2, 0) = 1;
  source1(0, 1) = 1234.8234;
  source1(1, 1) = 0.00001;
  source1(2, 1) = 352.412;
  source1(0, 2) = 1234.8234;
  source1(1, 2) = 352.412;
  source1(2, 2) = 9;

  source2(0, 0) = 352.412;
  source2(1, 0) = 1;
  source2(2, 0) = 0;
  source2(0, 1) = -1234.8234;
  source2(1, 1) = 0.00001;
  source2(2, 1) = 123456789;
  source2(0, 2) = -0.00001;
  source2(1, 2) = 0.412;
  source2(2, 2) = 9;
  S21Matrix result = source1 - source2;
  expected(0, 0) = 352.412 - 352.412;
  expected(1, 0) = 0;
  expected(2, 0) = 1;
  expected(0, 1) = 1234.8234 - -1234.8234;
  expected(1, 1) = 0.00001 - 0.00001;
  expected(2, 1) = 352.412 - 123456789;
  expected(0, 2) = 1234.8234 - -0.00001;
  expected(1, 2) = 352.412 - 0.412;
  expected(2, 2) = 0;
  ASSERT_EQ(result.EqMatrix(expected), 1);
}

TEST(Overload, MulNum1) {
  S21Matrix source(3, 3);
  source(0, 0) = 1.3;
  source(1, 0) = 89.13;
  source(2, 0) = 1;
  source(0, 1) = 9.1;
  source(1, 1) = 7;
  source(2, 1) = 1389;
  source(0, 2) = 8.123;
  source(1, 2) = 0;
  source(2, 2) = 0.0004;
  S21Matrix result = source * 1.234567;
  S21Matrix expected(3, 3);
  expected(0, 0) = 1.3 * 1.234567;
  expected(1, 0) = 89.13 * 1.234567;
  expected(2, 0) = 1 * 1.234567;
  expected(0, 1) = 9.1 * 1.234567;
  expected(1, 1) = 7 * 1.234567;
  expected(2, 1) = 1389 * 1.234567;
  expected(0, 2) = 8.123 * 1.234567;
  expected(1, 2) = 0 * 1.234567;
  expected(2, 2) = 0.0004 * 1.234567;
  ASSERT_EQ(result.EqMatrix(expected), 1);
}

TEST(Overload, MulNum2) {
  S21Matrix source(3, 3);
  source(0, 0) = 1.3;
  source(1, 0) = 89.13;
  source(2, 0) = 1;
  source(0, 1) = 9.1;
  source(1, 1) = 7;
  source(2, 1) = 1389;
  source(0, 2) = 8.123;
  source(1, 2) = 0;
  source(2, 2) = 0.0004;
  S21Matrix result = 1.234567 * source;
  S21Matrix expected(3, 3);
  expected(0, 0) = 1.3 * 1.234567;
  expected(1, 0) = 89.13 * 1.234567;
  expected(2, 0) = 1 * 1.234567;
  expected(0, 1) = 9.1 * 1.234567;
  expected(1, 1) = 7 * 1.234567;
  expected(2, 1) = 1389 * 1.234567;
  expected(0, 2) = 8.123 * 1.234567;
  expected(1, 2) = 0 * 1.234567;
  expected(2, 2) = 0.0004 * 1.234567;
  ASSERT_EQ(result.EqMatrix(expected), 1);
}

TEST(Overload, MulMat) {
  S21Matrix source1(3, 3), source2(3, 3), expected(3, 3);
  source1(0, 0) = 3;
  source1(0, 1) = 27;
  source1(0, 2) = 18;
  source1(1, 0) = 0;
  source1(1, 1) = 18;
  source1(1, 2) = 15;
  source1(2, 0) = 24;
  source1(2, 1) = 9;
  source1(2, 2) = 9;

  source2(0, 0) = 3;
  source2(0, 1) = 18;
  source2(0, 2) = 24;
  source2(1, 0) = 21;
  source2(1, 1) = 15;
  source2(1, 2) = 24;
  source2(2, 0) = 0;
  source2(2, 1) = 27;
  source2(2, 2) = 3;

  expected(0, 0) = 576;
  expected(0, 1) = 945;
  expected(0, 2) = 774;
  expected(1, 0) = 378;
  expected(1, 1) = 675;
  expected(1, 2) = 477;
  expected(2, 0) = 261;
  expected(2, 1) = 810;
  expected(2, 2) = 819;

  S21Matrix result = source1 * source2;
  ASSERT_EQ(result.EqMatrix(expected), 1);
}

TEST(Overload, SumCut) {
  S21Matrix source1(3, 3), source2(3, 3), expected(3, 3);
  source1(0, 0) = 352.412;
  source1(1, 0) = 1;
  source1(2, 0) = 1;
  source1(0, 1) = 1234.8234;
  source1(1, 1) = 0.00001;
  source1(2, 1) = 352.412;
  source1(0, 2) = 1234.8234;
  source1(1, 2) = 352.412;
  source1(2, 2) = 9;

  source2(0, 0) = 352.412;
  source2(1, 0) = 1;
  source2(2, 0) = 0;
  source2(0, 1) = -1234.8234;
  source2(1, 1) = 0.00001;
  source2(2, 1) = 123456789;
  source2(0, 2) = -0.00001;
  source2(1, 2) = 0.412;
  source2(2, 2) = 9;
  source1 += source2;
  expected(0, 0) = 352.412 + 352.412;
  expected(1, 0) = 1 + 1;
  expected(2, 0) = 1;
  expected(0, 1) = -1234.8234 + 1234.8234;
  expected(1, 1) = 0.00001 + 0.00001;
  expected(2, 1) = 352.412 + 123456789;
  expected(0, 2) = -0.00001 + 1234.8234;
  expected(1, 2) = 352.412 + 0.412;
  expected(2, 2) = 9 + 9;
  ASSERT_EQ(source1.EqMatrix(expected), 1);
}

TEST(Overload, SubCut) {
  S21Matrix source1(3, 3), source2(3, 3), expected(3, 3);
  source1(0, 0) = 352.412;
  source1(1, 0) = 1;
  source1(2, 0) = 1;
  source1(0, 1) = 1234.8234;
  source1(1, 1) = 0.00001;
  source1(2, 1) = 352.412;
  source1(0, 2) = 1234.8234;
  source1(1, 2) = 352.412;
  source1(2, 2) = 9;

  source2(0, 0) = 352.412;
  source2(1, 0) = 1;
  source2(2, 0) = 0;
  source2(0, 1) = -1234.8234;
  source2(1, 1) = 0.00001;
  source2(2, 1) = 123456789;
  source2(0, 2) = -0.00001;
  source2(1, 2) = 0.412;
  source2(2, 2) = 9;
  source1 -= source2;
  expected(0, 0) = 352.412 - 352.412;
  expected(1, 0) = 0;
  expected(2, 0) = 1;
  expected(0, 1) = 1234.8234 - -1234.8234;
  expected(1, 1) = 0.00001 - 0.00001;
  expected(2, 1) = 352.412 - 123456789;
  expected(0, 2) = 1234.8234 - -0.00001;
  expected(1, 2) = 352.412 - 0.412;
  expected(2, 2) = 0;
  ASSERT_EQ(source1.EqMatrix(expected), 1);
}

TEST(Overload, MulNum1Cut) {
  S21Matrix source(3, 3);
  source(0, 0) = 1.3;
  source(1, 0) = 89.13;
  source(2, 0) = 1;
  source(0, 1) = 9.1;
  source(1, 1) = 7;
  source(2, 1) = 1389;
  source(0, 2) = 8.123;
  source(1, 2) = 0;
  source(2, 2) = 0.0004;
  source *= 1.234567;
  S21Matrix expected(3, 3);
  expected(0, 0) = 1.3 * 1.234567;
  expected(1, 0) = 89.13 * 1.234567;
  expected(2, 0) = 1 * 1.234567;
  expected(0, 1) = 9.1 * 1.234567;
  expected(1, 1) = 7 * 1.234567;
  expected(2, 1) = 1389 * 1.234567;
  expected(0, 2) = 8.123 * 1.234567;
  expected(1, 2) = 0 * 1.234567;
  expected(2, 2) = 0.0004 * 1.234567;
  ASSERT_EQ(source.EqMatrix(expected), 1);
}

TEST(Overload, MulMatCut) {
  S21Matrix source1(3, 3), source2(3, 3), expected(3, 3);
  source1(0, 0) = 3;
  source1(0, 1) = 27;
  source1(0, 2) = 18;
  source1(1, 0) = 0;
  source1(1, 1) = 18;
  source1(1, 2) = 15;
  source1(2, 0) = 24;
  source1(2, 1) = 9;
  source1(2, 2) = 9;

  source2(0, 0) = 3;
  source2(0, 1) = 18;
  source2(0, 2) = 24;
  source2(1, 0) = 21;
  source2(1, 1) = 15;
  source2(1, 2) = 24;
  source2(2, 0) = 0;
  source2(2, 1) = 27;
  source2(2, 2) = 3;

  expected(0, 0) = 576;
  expected(0, 1) = 945;
  expected(0, 2) = 774;
  expected(1, 0) = 378;
  expected(1, 1) = 675;
  expected(1, 2) = 477;
  expected(2, 0) = 261;
  expected(2, 1) = 810;
  expected(2, 2) = 819;

  source1 *= source2;
  ASSERT_EQ(source1.EqMatrix(expected), 1);
}

TEST(Overload, Equity) {
  S21Matrix a(5, 5), b(5, 5);
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      a(i, j) = b(i, j) = i + j;
    }
  }
  ASSERT_EQ(a == b, 1);
  ASSERT_EQ(b == a, 1);
}

TEST(Overload, Constant) {
  int size = 5;
  S21Matrix a(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      a(i, j) = i + j;
    }
  }
  S21Matrix b = a;
  ASSERT_DOUBLE_EQ(a(1, 1), b(1, 1));
}