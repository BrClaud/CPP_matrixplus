#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(Constructor, InvalidParams) {
  EXPECT_THROW(S21Matrix matrix(-1, 0), std::invalid_argument);
}

TEST(SumMatrix, DiffSize) {
  S21Matrix source1(3, 2), source2(2, 3);
  EXPECT_THROW(source1.SumMatrix(source2), std::logic_error);
  EXPECT_THROW(source2.SumMatrix(source1), std::logic_error);
}

TEST(SubMatrix, DiffSize) {
  S21Matrix source1(3, 2), source2(2, 3);
  EXPECT_THROW(source1.SubMatrix(source2), std::logic_error);
  EXPECT_THROW(source2.SubMatrix(source1), std::logic_error);
}

TEST(Overload, GiveRvalue) {
  S21Matrix a(9, 9);
  S21Matrix b;
  b = std::move(a);
  EXPECT_EQ(b.EqMatrix(a), 0);
  EXPECT_EQ(a.GetRows(), 0);
  EXPECT_EQ(a.GetCols(), 0);
  EXPECT_THROW(a.GetValue(0, 0), std::range_error);
}

TEST(MulMatrix, WrongDimension) {
  S21Matrix source1(3, 2), source2(3, 2);
  EXPECT_THROW(source1.MulMatrix(source2), std::invalid_argument);
}

TEST(Determinant, WrongSize) {
  S21Matrix source(3, 9);
  EXPECT_THROW(source.Determinant(), std::invalid_argument);
}

TEST(Complements, SizeError) {
  S21Matrix source(2, 1);
  EXPECT_THROW(source.CalcComplements(), std::invalid_argument);
}

TEST(Complements, EmptyMatrix) {
  S21Matrix source;
  EXPECT_THROW(source.CalcComplements(), std::invalid_argument);
}

TEST(InverseMatrix, Uninitialized) {
  S21Matrix source;
  EXPECT_THROW(source.InverseMatrix(), std::runtime_error);
}

TEST(InverseMatrix, AllZero) {
  S21Matrix source(10, 10);
  EXPECT_THROW(source.InverseMatrix(), std::runtime_error);
}

TEST(InverseMatrix, WrongSize) {
  S21Matrix source(10, 9);
  EXPECT_THROW(source.InverseMatrix(), std::invalid_argument);
}

TEST(Inverse, Uncorrect) {
  S21Matrix source1(3, 2);
  EXPECT_THROW(source1.InverseMatrix(), std::invalid_argument);
}

TEST(Inverse, Uncorrect_1) {
  S21Matrix source1(3, 12);
  EXPECT_THROW(source1.InverseMatrix(), std::invalid_argument);
}

TEST(MulMatrix, SizeError) {
  S21Matrix source1(2, 2), source2(5, 5);
  EXPECT_THROW(source1.MulMatrix(source2), std::invalid_argument);
}

TEST(Index, OutOfRange) {
  S21Matrix source1(3, 3);
  EXPECT_THROW(source1(12, 4), std::invalid_argument);
}

TEST(Seter, LessThanRows) {
  S21Matrix source;
  EXPECT_THROW(source.SetRows(-3), std::invalid_argument);
}

TEST(Seter, LessThanCols) {
  S21Matrix source;
  EXPECT_THROW(source.SetCols(-3), std::invalid_argument);
}

TEST(Constructors, WrongIndex) {
  EXPECT_THROW(S21Matrix source(-1, -1), std::invalid_argument);
}

TEST(Constructors, Move) {
  S21Matrix source(3, 3);
  source(0, 0) = 99;
  source(2, 2) = -99;
  source(1, 1) = 13;
  S21Matrix moved_matrix(std::move(source));
  ASSERT_EQ(source.EqMatrix(moved_matrix), 0);
  ASSERT_EQ(moved_matrix(0, 0), 99);
  ASSERT_EQ(moved_matrix(1, 1), 13);
  ASSERT_EQ(moved_matrix(2, 2), -99);
}

TEST(Overload, DifferentDimensionsSum) {
  S21Matrix source1(3, 2), source2(1, 3);
  EXPECT_THROW(source1 + source2, std::logic_error);
}

TEST(Overload, DifferentDimensionsSub) {
  S21Matrix source1(3, 2), source2(1, 3);
  EXPECT_THROW(source1 - source2, std::logic_error);
}

TEST(MulMatrix, SizeErrorMUL) {
  S21Matrix source1(2, 2), source2(5, 5);
  EXPECT_THROW(source1 * source2, std::invalid_argument);
}