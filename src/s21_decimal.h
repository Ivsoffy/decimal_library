#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BITS_IN_BIG 192
#define BITS_IN_INT 32
#define POSITIVE_MASK 2147483648
#define MAX_VALUE 7.922816251426434e+28
#define MIN_VALUE 1e-28

typedef struct {
  unsigned int bits[4];
} s21_decimal;

typedef struct {
  unsigned int bits[7];
} s21_big_decimal;

typedef enum ArithmeticStatus {
  OK,       // 0
  MAX_VAL,  // 1
  MIN_VAL,  // 2
  DIV_ZERO  // 3
} ArithmeticStatus;

// Арифметические операторы
int s21_add(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_sub(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_mul(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_div(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_mod(s21_decimal value_1, s21_decimal value_2, s21_decimal *result);
int s21_floor(s21_decimal value, s21_decimal *result);
// Операторы сравнения
int s21_is_less(s21_decimal, s21_decimal);
int s21_is_less_or_equal(s21_decimal, s21_decimal);
int s21_is_greater(s21_decimal, s21_decimal);
int s21_is_greater_or_equal(s21_decimal, s21_decimal);
int s21_is_equal(s21_decimal, s21_decimal);
int s21_is_not_equal(s21_decimal, s21_decimal);
// Преобразователи
int s21_from_int_to_decimal(int src, s21_decimal *dst);
int s21_from_float_to_decimal(float src, s21_decimal *dst);
int s21_from_decimal_to_int(s21_decimal src, int *dst);
int s21_from_decimal_to_float(s21_decimal src, float *dst);
// Вспомогательные функции
void set_sign(int sign, s21_decimal *bit);
void get_sign(int *sign, s21_decimal bit);
void set_scale(int coef, s21_decimal *bit);
void get_scale(int *coef, s21_decimal bit);
int cmp_scale(s21_decimal bit1, s21_decimal bit2);
int set_bit(int number, s21_decimal *bit);
int get_bit(int number, s21_decimal *bit);
void mul_X(s21_decimal src, s21_decimal *rez);
void copy_decimal(s21_decimal *cp, s21_decimal src);
int compare(s21_decimal number1, s21_decimal number2);
int checkMinusZero(s21_decimal number1, s21_decimal number2);
// Big_decimal функции
int s21_is_greater_or_equal_big(s21_big_decimal value_1,
                                s21_big_decimal value_2);
void scale_norming(s21_decimal value_1, s21_decimal value_2, s21_big_decimal *a,
                   s21_big_decimal *b);
s21_big_decimal shift(s21_big_decimal value_1, int shift);
s21_big_decimal s21_minus_core(s21_big_decimal value_1, s21_big_decimal value_2,
                               s21_big_decimal *result);
int shift_10(int total_scale, s21_big_decimal *result);
int wrt(s21_big_decimal result);
s21_big_decimal s21_add_core(s21_big_decimal value_1, s21_big_decimal value_2,
                             s21_big_decimal *result);
int s21_mod_core(s21_big_decimal value_1, s21_big_decimal value_2,
                 s21_big_decimal *result);
int s21_big_decimal_to_decimal(s21_big_decimal value, s21_decimal *result);
void to_big_decimal(s21_decimal small, s21_big_decimal *big);
void equalize(s21_big_decimal *number1, s21_big_decimal *number2);
int mul_bit(s21_big_decimal number1, s21_big_decimal number2,
            s21_big_decimal *result);
int sum_bit(s21_big_decimal number1, s21_big_decimal number2,
            s21_big_decimal *result);
int big_get_bit(s21_big_decimal number, int index);
void big_set_bit(s21_big_decimal *number, int index, int value);
int big_get_len(s21_big_decimal number);
int big_get_scale(s21_big_decimal number);
void big_set_scale(s21_big_decimal *number, int scale);
int shift_left_compare(s21_big_decimal *number);
int first_mismatch(s21_big_decimal number1, s21_big_decimal number2);
int big_get_sign(s21_big_decimal number);
void big_set_sign(s21_big_decimal *number, int value);

// Другие функции
int s21_floor(s21_decimal value, s21_decimal *result);
int s21_round(s21_decimal value, s21_decimal *result);
int s21_truncate(s21_decimal value, s21_decimal *result);
int s21_negate(s21_decimal value, s21_decimal *result);
