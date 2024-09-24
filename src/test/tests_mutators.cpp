#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(Mutators, Setter1) {
  int size = 9;
  S21Matrix source(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      source(i, j) = i + j;
    }
  }
  size--;
  source.SetRows(size);
  ASSERT_EQ(source.GetRows(), size);
}

TEST(Mutators, Setter2) {
  int size = 9;
  S21Matrix source(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      source(i, j) = i + j;
    }
  }
  size++;
  source.SetRows(size);
  ASSERT_EQ(source.GetRows(), size);
  ASSERT_DOUBLE_EQ(source(size - 1, size - 2), 0);
}

TEST(Mutators, Setter3) {
  int size = 9;
  S21Matrix source(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      source(i, j) = i + j;
    }
  }
  size--;
  source.SetCols(size);
  ASSERT_EQ(source.GetCols(), size);
}

TEST(Mutators, Setter4) {
  int size = 9;
  S21Matrix source(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      source(i, j) = i + j;
    }
  }
  size++;
  source.SetCols(size);
  ASSERT_EQ(source.GetCols(), size);
  ASSERT_DOUBLE_EQ(source(size - 2, size - 1), 0);
}