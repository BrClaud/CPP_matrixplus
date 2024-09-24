#ifndef S21_MATRIX_OOP
#define S21_MATRIX_OOP

#include <stdio.h>

#include <cmath>
#include <iostream>
// #include <math.h>

using namespace std;

#define INCORRECT_MATRIX 1
#define CALCULATION_ERROR 2
#define SUCCESS 1
#define FAILURE 0
#define EPS 1e-6

class S21Matrix {
 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  double **matrix_;  // Pointer to the memory where the matrix is allocated

 public:
  // Default constructor 2X2
  S21Matrix();
  // Параметризированный конструктор с количеством строк и столбцов.
  S21Matrix(int rows, int cols);
  // Конструктор копирования.
  S21Matrix(const S21Matrix &other);
  // Конструктор переноса.
  S21Matrix(S21Matrix &&other);
  // создание матрицы заполненной нулями
  void create_matrix(int rows, int cols);
  // Destructor
  ~S21Matrix();

  // Проверяет матрицы на равенство между собой.
  bool EqMatrix(const S21Matrix &other);
  // Прибавляет вторую матрицы к текущей
  void SumMatrix(const S21Matrix &other);
  // Вычитает из текущей матрицы другую
  void SubMatrix(const S21Matrix &other);
  // Умножает текущую матрицу на число.
  void MulNumber(const double num);
  // Умножает текущую матрицу на вторую.
  void MulMatrix(const S21Matrix &other);
  // Создает новую транспонированную матрицу из текущей и возвращает ее.
  S21Matrix Transpose();
  // Вычисляет матрицу алгебраических дополнений текущей матрицы и возвращает ее
  S21Matrix CalcComplements();
  // Вычисляет и возвращает определитель текущей матрицы.
  double Determinant();
  // 2х2
  double Determinant_default();
  // считает минор
  S21Matrix Calc_minor(int ind_i, int ind_j);
  // Вычисляет и возвращает обратную матрицу.
  S21Matrix InverseMatrix();

  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix operator*(const double num);
  bool operator==(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(double num);
  friend S21Matrix operator*(double num, const S21Matrix &other);
  double &operator()(int rows, int cols);
  S21Matrix &operator=(S21Matrix &other);
  S21Matrix &operator=(S21Matrix &&other);
  // информация о матрице
  // void info();
  // проверка на нормальные значения
  // int check_matrix(const S21Matrix& other);

  // int nan_or_inf(const S21Matrix& other);
  // геттеры
  int GetRows();
  int GetCols();
  double GetValue(int i, int j);

  // сеттеры, трогать очень аккуратно)
  void SetValue(double val, int i, int j);
  void SetRows(int new_rows);
  void SetCols(int new_cols);

  // копирование матрицы
  void CopyMatrix(const S21Matrix &other);
};
#endif