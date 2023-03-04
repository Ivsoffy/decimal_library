#include "s21_decimal.h"

// // errokele
// int s21_wrt(s21_decimal result) {
//   for (int j = 0; j < 4; j++) {
//     for (int i = 0; i < 32; i++) {
//       if ((result.bits[3 - j] >> (31 - i)) & 1) {
//         printf("1");
//       } else {
//         printf("0");
//       }
//     }
//     printf("\n");
//   }
//   printf("\n");
//   return 0;
// }

// int main (){
//     s21_decimal src1 = {0};
//     s21_decimal src2 = {0};
//     // scanf("%u", &a.bits[0]);
//     // scanf("%u", &b.bits[0]);
//     //             00000000000001011001000011001100
//   src1.bits[0] = 0;
//   src1.bits[1] = 0;
//   src1.bits[2] = 1;
//   src1.bits[3] = 0b10000000000000000000000000000000;
  
//   src2.bits[0] = 0;
//   src2.bits[1] = 0;
//   src2.bits[2] = 2;
//   src2.bits[3] = 0b00000000000000000000000000000000;
//   // 0b000000000001011001000011001011

// // 00000000000000000000000000000000
// // 11111111111111111111111111111111
// // 11111111111000001011111000001110
// // 00010111010101110100100010001000s


    
    
//     // b.bits[3] = b.bits[3] | (1 << 31);
//     //  a.bits[3] = a.bits[3] | (1 << 31);
//    //  b.bits[3] = 4;
//     //b.bits[6] = b.bits[6] | (1 << 17);
//     s21_decimal * result = malloc(sizeof(int) * 4);
//     // s21_big_decimal_to_decimal(b, result);
//     printf("inf?  %u\n", s21_div(src1, src2, result));
//     // wrt(*result);
//     // printf("inf? a - b %u\n",s21_minus(a, b, result));
//     // wrt(*result)

//     // printf("inf? a /b  %u\n",s21_div_core(a, b, result ));
//     // wrt(*result);
//     // printf("inf? a * b %u\n",s21_mul (a, b, result));
//     for (int i = 0; i < 3; i++){
//       // result->bits[i] =b.bits[i] ;
//     }
//      s21_wrt(*result);
//     return 0;

// }

void set_sign(int sign, s21_decimal *bit) {
  if (sign) {
    unsigned int mask = POSITIVE_MASK;
    bit->bits[3] = mask | bit->bits[3];
  } else {
    unsigned int mask = POSITIVE_MASK;
    mask = ~mask;
    bit->bits[3] = mask & bit->bits[3];
  }
}

void get_sign(int *sign, s21_decimal bit) {
  if (POSITIVE_MASK & bit.bits[3]) {
    *sign = 1;
  } else {
    *sign = 0;
  }
}

void set_scale(int coef, s21_decimal *bit) {
  bit->bits[3] = bit->bits[3] & 0xFF00FFFF;
  bit->bits[3] = bit->bits[3] | (coef * 65536);  // 2^16
}

void get_scale(int *coef, s21_decimal bit) {
  *coef = (bit.bits[3] << 1) >> 17;  // 2^24s
}

int get_bit(int num, s21_decimal *value) {
  return (value->bits[num / 32] >> (num % 32)) & 1u;
}

int s21_from_int_to_decimal(int src, s21_decimal *dst) {
  int status = 0;
  for (int i = 0; i < 4; i++) {
    dst->bits[i] = 0;
  }
  (src >= 0) ? set_sign(0, dst) : set_sign(1, dst);
  if (src < 0) src *= -1;
  dst->bits[0] = src;
  return status;
}

int s21_from_decimal_to_int(s21_decimal src, int *dst) {
  int status = 0;
  int sign = 0;
  s21_decimal tmp = {0};
  s21_truncate(src, &tmp);
  get_sign(&sign, tmp);
  if ((tmp.bits[2] == 0) && (tmp.bits[1] == 0)) {
    *dst = tmp.bits[0];
    if (sign) {
      *dst *= -1;
    }
  } else {
    *dst = 0;
    status = 1;
  }
  return status;
}

int s21_from_float_to_decimal(float src, s21_decimal *dst) {
  int status = 0;
  for (int i = 0; i < 4; i++) {
    dst->bits[i] = 0;
  }
  if ((src != src) | (fabsf(src) > MAX_VALUE) | (fabsf(src) == INFINITY) |
      (fabsf(src) < MIN_VALUE)) {
    if (fabsf(src) != 0) {
      status = 1;
    }
  } else {
    char str[100] = {0};
    snprintf(str, 100, "%.6e", src);
    int sign = (src > 0 ? 0 : 1);
    src = fabsf(src);
    char numb[100];
    int j = 0;
    for (size_t i = 0; i < strlen(str); i++) {
      if (str[i] == 'e') {
        numb[j] = '\0';
        break;
      }
      if (str[i] != '.') {
        numb[j] = str[i];
        j++;
      }
    }
    int scale = 0;
    int pred_float = 0;
    sscanf(numb, "%d", &pred_float);
    sscanf(str, "%*d.%*de%d", &scale);
    scale++;
    dst->bits[0] = abs(pred_float);
    scale = (scale >= 7) ? 0 : 7 - scale;
    set_sign(sign, dst);
    set_scale(scale, dst);
  }
  return status;
}

int s21_from_decimal_to_float(s21_decimal src, float *dst) {
  int status = 0;
  int scale = 0;
  get_scale(&scale, src);
  if (dst == NULL || scale > 28) {
    status = 1;
  } else {
    double result = 0;
    for (int i = 0; i < 96; i++) {
      if (get_bit(i, &src) != 0) {
        result += pow(2, i);
      }
    }
    while (scale != 0) {
      result /= 10;
      scale--;
    }
    *dst = 0;
    *dst = (float)result;
    int sign;
    get_sign(&sign, src);
    if (sign == 1) {
      *dst *= -1;
    }
  }
  return status;
}

int s21_truncate(s21_decimal value, s21_decimal *result) {  // cp result
  int status = 0;
    int scale = 0;
  s21_decimal result1 = {0};
  // int tmp = 0;
  copy_decimal(&result1, value);
  // get_scale(&scale, value);
  // tmp = (value.bits[3] << 1) >> 17;
  scale = (value.bits[3] << 1) >> 17;
  for ( int i = 0; i < scale; i++) {
    mul_X(value, &result1);
    copy_decimal(&value, result1);
  }
  if (!scale) {
    copy_decimal(&result1, value);
  }
  copy_decimal(result, result1);
  set_scale(0, result);
  return status;
}

void copy_decimal(s21_decimal *cp, s21_decimal src) {
  for (int i = 0; i <= 3; i++) {
    cp->bits[i] = src.bits[i];
  }
}

void mul_XX(s21_big_decimal src, s21_big_decimal *rez) {
  unsigned long long modd = 0;
  unsigned long long buff = 0;
  for (int i = 0; i < 6; i++) {
    rez->bits[i] = 0;
  }
  rez->bits[6] = src.bits[6];
  for (int i = 5; i >= 0; i--) {
    modd = (unsigned long long)(buff + src.bits[i]) % 10;
    rez->bits[i] = (unsigned long long)(buff + src.bits[i]) / 10;
    buff = modd << 32;
  }
}

void mul_X(s21_decimal src, s21_decimal *rez) {
  unsigned long long modd = 0;
  unsigned long long buff = 0;
  int scale = 0;
  for (int i = 0; i < 3; i++) {
    rez->bits[i] = 0;
  }
  rez->bits[3] = src.bits[3];
  for (int i = 2; i >= 0; i--) {
    modd = (buff + src.bits[i]) % 10;
    rez->bits[i] = (buff + src.bits[i]) / 10;
    buff = modd << 32;
  }
  get_scale(&scale, src);
  set_scale(scale + 1, rez);
}

int s21_negate(s21_decimal value, s21_decimal *result) {
  int status = 0, sign;
  copy_decimal(result, value);
  get_sign(&sign, *result);  // error status
  if (sign)
    set_sign(0, result);
  else
    set_sign(1, result);
  return status;
}

int s21_round(s21_decimal value, s21_decimal *result) {
  int status = 0;
  int scale = 0, sign = 0;
  int flag = 0;
  for (int i = 0; i < 3; i++) {
    if (value.bits[i] != 0) flag = 1;
  }
  get_sign(&sign, value);
  if (flag) {
    int i = 0;
    s21_decimal result1 = {0};
    s21_decimal one = {0};
    s21_decimal ten = {0};
    s21_decimal mantissa_src = {0};
    s21_from_int_to_decimal(10, &ten);
    s21_from_int_to_decimal(1, &one);
    copy_decimal(&result1, value);
    get_scale(&scale, value);

    set_sign(0, &value);
    s21_truncate(result1, &result1);
    set_sign(0, &result1);
    int flag = 0;
    for (int i = 0; i < 3; i++) {
      if (result1.bits[i] != 0) flag = 1;
    }
    if (flag) {
      s21_mod(value, result1, &mantissa_src);
      get_scale(&scale, value);
      set_scale(scale, &mantissa_src);
    } else {
      for (int i = 0; i < 4; i++) {
        mantissa_src.bits[i] = 0;  //
      }
    }
    for (int i = 0; i < scale; i++) {
      mul_X(value, &result1);
      copy_decimal(&value, result1);
    }
    set_scale(0, &value);
    s21_decimal polovinka;
    for (int i = 0; i < 4; i++) {
      polovinka.bits[i] = 0;
    }
    polovinka.bits[0] = 5;
    set_scale(1, &polovinka);
    if (s21_is_greater_or_equal(mantissa_src, polovinka)) {
      s21_add(one, value, &result1);
      copy_decimal(&value, result1);
    }
    if (!i) {
      copy_decimal(&result1, value);
    }
    copy_decimal(result, result1);
    result->bits[3] = 0;
    set_sign(sign, result);
    set_scale(0, &value);
  } else {
    for (int i = 0; i < 4; i++) {
      result->bits[i] = 0;
    }
    set_sign(sign, result);
  }
  return status;
}

int s21_floor(s21_decimal value, s21_decimal *result) {
  int status = 0, sign = 0;
  s21_decimal value_cp = {0};
  s21_truncate(value, &value_cp);
  get_sign(&sign, value);
  if (s21_is_not_equal(value, value_cp)) {
    if (sign) {
      s21_decimal one;
      s21_from_int_to_decimal(1, &one);
      s21_sub(value_cp, one, result);
    } else {
      copy_decimal(result, value_cp);
    }
  } else {
    copy_decimal(result, value_cp);
  }
  return status;
}
// harodonf

int s21_add(s21_decimal value_1, s21_decimal value_2, s21_decimal *resuls) {
  s21_big_decimal *a = malloc(sizeof(unsigned int) * 7);
  s21_big_decimal *b = malloc(sizeof(unsigned int) * 7);
  int flag = 0;
  s21_big_decimal *resuls_tmp = malloc(sizeof(unsigned int) * 7);
  int inf = 0;
  for (int i = 0; i < 7; i++) {
    a->bits[i] = 0;
    resuls_tmp->bits[i] = 0;
    b->bits[i] = 0;
  }
  if ((value_1.bits[3] >> 31) ^ (value_2.bits[3] >> 31)) {
    if ((value_1.bits[3] >> 31)) {
      value_1.bits[3] = (value_1.bits[3] << 1) >> 1;
      flag++;
    } else {
      value_2.bits[3] = (value_2.bits[3] << 1) >> 1;
    }
    s21_sub(value_1, value_2, resuls);
    if (flag) {
      if (resuls->bits[3] >> 31) {
        resuls->bits[3] = (resuls->bits[3] << 1) >> 1;
      } else {
        resuls->bits[3] = resuls->bits[3] | (1 << 31);
      }
    }
  } else {
    scale_norming(value_1, value_2, a, b);
    s21_add_core(*a, *b, resuls_tmp);
    resuls_tmp->bits[6] = (((a->bits[6] << 1) >> 1)) | resuls_tmp->bits[6];
    inf = s21_big_decimal_to_decimal(*resuls_tmp, resuls);
    if ((value_1.bits[3] >> 31) & 1) {
      resuls->bits[3] = resuls->bits[3] | (1 << 31);
    }
  }
  if ((resuls->bits[3] >> 31) && inf) {
    inf = 2;
  }
  free(a);
  free(b);
  free(resuls_tmp);
  return inf;
}

int s21_sub(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
  s21_big_decimal *a = malloc(sizeof(unsigned int) * 7);
  s21_big_decimal *b = malloc(sizeof(unsigned int) * 7);
  s21_big_decimal *resuls_tmp = malloc(sizeof(unsigned int) * 7);
  int flag = 0;
  int inf = 0;
  for (int i = 0; i < 7; i++) {
    a->bits[i] = 0;
    b->bits[i] = 0;
    resuls_tmp->bits[i] = 0;
  }
  if ((value_1.bits[3] >> 31) ^ (value_2.bits[3] >> 31)) {
    if ((value_1.bits[3] >> 31)) {
      value_1.bits[3] = (value_1.bits[3] << 1) >> 1;
    } else {
      value_2.bits[3] = (value_2.bits[3] << 1) >> 1;
      flag++;
    }
    s21_add(value_1, value_2, result);
    if (!flag) {
      if (result->bits[3] >> 31) {
        result->bits[3] = (result->bits[3] << 1) >> 1;
      } else {
        result->bits[3] = result->bits[3] | (1 << 31);
      }
    }
  } else {
    scale_norming(value_1, value_2, a, b);

    if (s21_is_greater_or_equal_big(*a, *b)) {
      s21_minus_core(*a, *b, resuls_tmp);
    } else {
      s21_minus_core(*b, *a, resuls_tmp);
    }
    inf = s21_big_decimal_to_decimal(*resuls_tmp, result);
    if (s21_is_greater(value_2, value_1)) {
      result->bits[3] = result->bits[3] | (1 << 31);
    }
  }
  if ((result->bits[3] >> 31) && inf) {
    inf = 2;
  }
  free(resuls_tmp);
  free(a);
  free(b);
  return inf;
}
int s21_is_greater_or_equal_big(s21_big_decimal value_1,
                                s21_big_decimal value_2) {  //

  int tmp = 0;
  int flag = 1;
  // wrt(value_1)
  for (int i = 0; i < 6; i++) {
    if (value_1.bits[i] == value_2.bits[i]) {
      tmp = 1;

    } else {
      flag = 0;
      tmp = 0;
      break;
    }
  }
  if (flag == 0) {
    for (int i = 0; i < 1; i++) {
      if (value_1.bits[5] < value_2.bits[5]) {
        tmp = 0;
        break;
      } else if (value_1.bits[5] > value_2.bits[5]) {
        tmp = 1;
        break;
      } else if (value_1.bits[4] < value_2.bits[4]) {
        tmp = 0;
        break;
      } else if (value_1.bits[4] > value_2.bits[4]) {
        tmp = 1;
        break;
      } else if (value_1.bits[3] < value_2.bits[3]) {
        tmp = 0;
        break;
      } else if (value_1.bits[3] > value_2.bits[3]) {
        tmp = 1;
        break;
      } else if (value_1.bits[2] < value_2.bits[2]) {
        tmp = 0;
        break;
      } else if (value_1.bits[2] > value_2.bits[2]) {
        tmp = 1;
        break;

      } else if (value_1.bits[1] < value_2.bits[1]) {
        tmp = 0;
        break;
      } else if (value_1.bits[1] > value_2.bits[1]) {
        tmp = 1;
        break;
      } else if (value_1.bits[0] < value_2.bits[0]) {
        tmp = 0;
        break;
      } else if (value_1.bits[0] > value_2.bits[0]) {
        tmp = 1;
        break;
      }
    }
  }
  return tmp;
}
void scale_norming(s21_decimal value_1, s21_decimal value_2, s21_big_decimal *a,
                   s21_big_decimal *b) {
  unsigned int total_scale = 0;
  unsigned int max_scale = 0;
  unsigned int scale_1 = 0;
  for (int i = 16; i < 24; i++) {
    if ((value_1.bits[3] >> i) & 1) {
      scale_1 = (value_1.bits[3] | (1 << i));
    }
  }
  scale_1 = scale_1 << 1;
  scale_1 = scale_1 >> 17;

  unsigned int scale_2 = 0;
  for (int i = 16; i < 24; i++) {
    if ((value_2.bits[3] >> i) & 1) {
      scale_2 = (value_2.bits[3] | (1 << i));
    }
  }
  scale_2 = scale_2 << 1;
  scale_2 = scale_2 >> 17;

  for (int i = 0; i < 3; i++) {
    a->bits[i] = value_1.bits[i];
    b->bits[i] = value_2.bits[i];
  }

  if (scale_1 < scale_2) {
    max_scale = scale_2;
    total_scale = scale_2 - scale_1;
    shift_10(total_scale, a);

  } else if (scale_1 != 0 || scale_2 != 0) {
    max_scale = scale_1;
    total_scale = scale_1 - scale_2;
    shift_10(total_scale, b);
  }
  a->bits[6] = a->bits[6] | (max_scale << (16));
  b->bits[6] = b->bits[6] | (max_scale << (16));
}

s21_big_decimal shift(s21_big_decimal value_1, int shift) {
  int shift_tmp = 0;
  s21_big_decimal tmp = {0};
  s21_big_decimal tmp1 = {0};
  shift_tmp = shift;
  if (shift == 33) {
    tmp = tmp1;
  }
  if (shift != 0) {
    for (int i = 0; i < 6; i++) {
      tmp.bits[i] = value_1.bits[i];
      if (shift >= 32) {
        shift = 0;
      } else {
        tmp1.bits[i] = (tmp.bits[i] << shift);
      }
    }
    shift = shift_tmp;
    for (int j = 0; j < 5; j++) {
      for (int i = 0; i < shift; i++) {
        if ((tmp.bits[j] >> (31 - i)) & 1) {
          tmp1.bits[j + 1] = (tmp1.bits[j + 1] | (1 << (shift - 1 - i)));
        }
      }
    }
  }
  return (shift != 0) ? tmp1 : value_1;
}

s21_big_decimal s21_minus_core(s21_big_decimal value_1, s21_big_decimal value_2,
                               s21_big_decimal *result) {
  int over_value = 0;

  for (int i = 0; i < 6; i++) {
    result->bits[i] = 0;
  }
  result->bits[6] = (value_1.bits[6] << 1) >> 1;
  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 32; i++) {
      if ((((value_1.bits[j] >> i) & (value_2.bits[j] >> i)) & 1) &&
          over_value == 1) {
        over_value = 1;
        result->bits[j] = result->bits[j] | (1 << i);
      }

      if (((((value_1.bits[j] >> i) & 1) && (!((value_2.bits[j] >> i) & 1)))) &&
          over_value == 0) {
        over_value = 0;
        result->bits[j] = result->bits[j] | (1 << i);
      }

      if (((((value_1.bits[j] >> i) & 1) && (!((value_2.bits[j] >> i) & 1)))) &&
          over_value == 1) {
        over_value = 0;
      }

      if ((((!((value_1.bits[j] >> i) & 1)) && ((value_2.bits[j] >> i) & 1))) &&
          over_value == 1) {
        over_value = 1;
      }

      if ((((!((value_1.bits[j] >> i) & 1)) && ((value_2.bits[j] >> i) & 1))) &&
          over_value == 0) {
        over_value = 1;
        result->bits[j] = result->bits[j] | (1 << i);
      }

      if ((!((value_1.bits[j] >> i) & 1) && !((value_2.bits[j] >> i) & 1)) &&
          over_value == 1) {
        over_value = 1;
        result->bits[j] = (result->bits[j] | (1 << i));
      }
    }
  }

  return *result;
}

int shift_10(int total_scale, s21_big_decimal *result) {
  s21_big_decimal a_s = {0};
  s21_big_decimal b_s = {0};
  s21_big_decimal *tmp = malloc(sizeof(unsigned int) * 7);
  for(int i = 0; i < 7; i++){
    tmp->bits[i] = 0;
  }

  for (int i = 0; i < 3; i++) {
    tmp->bits[i] = result->bits[i];
  }

  for (int i = 0; i < total_scale; i++) {
    a_s = shift(*result, 3);
    b_s = shift(*result, 1);
    s21_add_core(a_s, b_s, result);
  }

  free(tmp);
  return 0;
}

s21_big_decimal s21_add_core(s21_big_decimal value_1, s21_big_decimal value_2,
                             s21_big_decimal *result) {
  int over_value = 0;

  for (int i = 0; i < 6; i++) {
    result->bits[i] = 0;
  }

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 32; i++) {
      if ((((value_1.bits[j] >> i) & (value_2.bits[j] >> i)) & 1) &&
          over_value == 1) {
        over_value = 1;
        result->bits[j] = result->bits[j] | (1 << i);
      }
      if ((((value_1.bits[j] >> i) & (value_2.bits[j] >> i)) & 1) &&
          over_value == 0) {
        over_value = 1;
      }
      if ((((value_1.bits[j] >> i) ^ (value_2.bits[j] >> i)) & 1) &&
          over_value == 1) {
        over_value = 1;
      }
      if ((((value_1.bits[j] >> i) ^ (value_2.bits[j] >> i)) & 1) &&
          over_value == 0) {
        over_value = 0;
        result->bits[j] = result->bits[j] | (1 << i);
      }
      if ((!((value_1.bits[j] >> i) & 1) && !((value_2.bits[j] >> i) & 1)) &&
          over_value == 1) {
        over_value = 0;
        result->bits[j] = (result->bits[j] | (1 << i));
      }
    }
  }
  return *result;
}

int s21_mul(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
  s21_big_decimal *a = malloc(sizeof(unsigned int) * 7);
  s21_big_decimal *b = malloc(sizeof(unsigned int) * 7);
  s21_big_decimal *result_tmp = malloc(sizeof(unsigned int) * 7);
  for (int i = 0; i < 7; i++) {
    a->bits[i] = 0;
    b->bits[i] = 0;
    result_tmp->bits[i] = 0;
  }
  int inf = 0;

  s21_big_decimal tmp = {0};
  s21_big_decimal value_1_tmp = {0};
  for (int i = 0; i < 3; i++) {
    a->bits[i] = value_1.bits[i];
    b->bits[i] = value_2.bits[i];
  }
  for (int i = 0; i < 6; i++) {
    result_tmp->bits[i] = 0;
  }

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 32; i++) {
      if ((b->bits[j] >> i) & 1) {
        value_1_tmp = shift(*a, (i + 32 * j));
        tmp = s21_add_core(tmp, value_1_tmp, result_tmp);
      }
    }
  }
  int scale1, scale2;
  get_scale(&scale1, value_1);
  get_scale(&scale2, value_2);
  tmp.bits[6] = (scale1 + scale2) << 16;
  inf = s21_big_decimal_to_decimal(tmp, result);
  // printf("%d\n", result->bits[3]);
  free(a);
  free(b);
  free(result_tmp);
  if ((value_1.bits[3] >> 31) ^ (value_2.bits[3] >> 31)) {
    result->bits[3] = result->bits[3] | (1 << 31);
    if (inf == 1) {
      inf = 2;
    }
  }
  return inf;
}

int s21_mod_core(s21_big_decimal value_1, s21_big_decimal value_2,
                 s21_big_decimal *result) {
  int inf = 0;
  if (value_2.bits[0] == 0 && value_2.bits[1] == 0 && value_2.bits[2] == 0) {
    inf = 1;
  } else {
    s21_big_decimal tmp = {0};
    s21_big_decimal count = {0};
    s21_big_decimal count_tmp = {0};

    count.bits[0] = 1;

    s21_big_decimal one = {0};
    s21_big_decimal value_1_tmp = {0};
    s21_big_decimal value_2_tmp = {0};
    s21_big_decimal value_2_ = {0};
    int flag = 1;
    one.bits[0] = 1;

    for (int i = 0; i < 6; i++) {
      value_2_tmp.bits[i] = value_2.bits[i];
      value_2_.bits[i] = value_2.bits[i];
      value_1_tmp.bits[i] = value_1.bits[i];
    }
    for (int i = 0; i < 6; i++) {
      result->bits[i] = 0;
    }

    while (flag) {
      tmp = shift(value_2_tmp, 1);

      if (s21_is_greater_or_equal_big(value_1_tmp, tmp)) {
        value_2_tmp = shift(value_2_tmp, 1);
        count = shift(count, 1);
        continue;
      } else {
        if (s21_is_greater_or_equal_big(value_1_tmp, value_2_tmp)) {
          value_1_tmp = s21_minus_core(value_1_tmp, value_2_tmp, result);
        }
        if (s21_is_greater_or_equal_big(value_1_tmp, value_2_)) {
          value_2_tmp = value_2_;
          tmp = value_2_;
          count_tmp = s21_add_core(count, count_tmp, result);
          count = one;
          continue;
        }
      }
      flag = 0;
      count = s21_add_core(count, count_tmp, result);
    }
    for (int i = 0; i < 6; i++) {
      result->bits[i] = value_1_tmp.bits[i];
    }
  }
  return inf;
}

int s21_div(s21_decimal value_1, s21_decimal value_2, s21_decimal *result_tmp) {
  int inf = 0;
  int sign1 = 0, sign2 = 0;
  get_sign(&sign1, value_1);
  get_sign(&sign2, value_2);

  if (0 == value_2.bits[0] && 0 == value_2.bits[1] && 0 == value_2.bits[2]) {
    inf = 3;
  } else {
    s21_big_decimal *a = malloc(sizeof(unsigned int) * 7);
    s21_big_decimal *b = malloc(sizeof(unsigned int) * 7);
    s21_big_decimal *result = malloc(sizeof(unsigned int) * 7);
    for(int i = 0; i < 7; i++){
      a->bits[i] = 0;
      b->bits[i] = 0;
      result->bits[i] = 0;
    }
    scale_norming(value_1, value_2, a, b);
    s21_big_decimal tmp = {0};
    s21_big_decimal count = {0};
    s21_big_decimal count_tmp = {0};
    s21_big_decimal count_post_point = {0};
    s21_big_decimal count_after_point = {0};
    long int qwe = 0;

    s21_big_decimal null = {0};

    count.bits[0] = 1;

    s21_big_decimal one = {0};

    s21_big_decimal value_1_tmp = {0};
    s21_big_decimal value_2_tmp = {0};
    s21_big_decimal value_2_ = {0};
    int flag = 1;
    one.bits[0] = 1;

    for (int i = 0; i < 6; i++) {
      value_2_tmp.bits[i] = b->bits[i];
      value_2_.bits[i] = b->bits[i];
      value_1_tmp.bits[i] = a->bits[i];
    }
    for (int i = 0; i < 6; i++) {
      result->bits[i] = 0;
    }
    int c_flag = 0;
    int i = 0;
    int suka = 0;
    int scale = 0;
    int scale_sum = 0;

    for (int w = 0; w < 2; w++) {
      while (flag) {
        tmp = shift(value_2_tmp, 1);
        qwe++;
        if (s21_is_greater_or_equal_big(value_1_tmp, tmp)) {
          value_2_tmp = shift(value_2_tmp, 1);
          count = shift(count, 1);
          continue;
        } else {
          if (s21_is_greater_or_equal_big(value_1_tmp, value_2_)) {
            if (s21_is_greater_or_equal_big(value_1_tmp, value_2_tmp)) {
              value_1_tmp = s21_minus_core(value_1_tmp, value_2_tmp, result);
            }
            value_2_tmp = value_2_;
            tmp = value_2_;

            if (c_flag) {
              if (suka) {
                for (int j = 0; j < scale; j++) {
                  count_tmp = s21_add_core(shift(count_tmp, 1),
                                           shift(count_tmp, 3), result);
                  scale_sum++;
                  suka = 0;
                }
                scale = 0;
              }
              count_tmp = s21_add_core(count, count_tmp, result);
              count = one;
              value_2_tmp = value_2_;
            }

            if (!c_flag) {
              count_tmp = s21_add_core(count, count_tmp, result);
            }
            count = one;
            continue;
          } else if (c_flag) {
            value_1_tmp = s21_add_core(shift(value_1_tmp, 1),
                                       shift(value_1_tmp, 3), result);
            scale++;

            suka = 1;
            i++;
            if (i < 28) {
              continue;
            }
          }
          flag = 0;
        }
      }
      if (!c_flag) {
        count_after_point = s21_add_core(count, count_tmp, result);
        value_2_tmp = value_2_;
        c_flag = 1;
        flag = 1;
        count = one;
        count_tmp = null;
      } else {
        count_post_point = s21_add_core(null, count_tmp, result);
      }
    }
    count_after_point = s21_minus_core(count_after_point, one, result);
    for (int i = 0; i < scale_sum; i++) {
      count_after_point = s21_add_core(shift(count_after_point, 1),
                                       shift(count_after_point, 3), result);
    }
    count_post_point =
        s21_add_core(count_post_point, count_after_point, result);
    for (int i = 0; i < 6; i++) {
      result->bits[i] = count_post_point.bits[i];
    }

    inf = s21_big_decimal_to_decimal(*result, result_tmp);
    if (((value_1.bits[3] >> 31) ^ (value_2.bits[3] >> 31)) && inf == 1) {
      inf = 2;
    }
    result_tmp->bits[3] = (result_tmp->bits[3] | (scale_sum << 16));
    result_tmp->bits[3] = result_tmp->bits[3] | ((sign1 ^ sign2) << 31);
    free(a);
    free(b);
    free(result);
  }
  return inf;
}

int s21_big_decimal_to_decimal(s21_big_decimal value, s21_decimal *result) {
  result->bits[3] = value.bits[6];
  value.bits[6] = value.bits[6] << 1;
  value.bits[6] = value.bits[6] >> 1;

  s21_big_decimal ten = {0};
  s21_big_decimal one = {0};
  int inf = 0;
  one.bits[0] = 1;
  ten.bits[0] = 10;

  s21_big_decimal tmp = {0};
  s21_big_decimal ost = {0};


  tmp = value;
  int scale = 0;
  scale = value.bits[6] >> 16;
  while (scale != 0) {
    if (tmp.bits[3] == 0 && tmp.bits[4] == 0 && tmp.bits[5] == 0) {
      for (int i = 0; i < 3; i++) {
        result->bits[i] = tmp.bits[i];
      }
      break;
    } else {
      scale--;
      tmp.bits[6] = scale << 16;
      s21_mod_core(tmp, ten, &ost);
      if (ost.bits[0] == 5) {
        mul_XX(tmp, &tmp);
        if (tmp.bits[0] % 2) {
          tmp = s21_add_core(tmp, one, &tmp);
        }
      } else {
        s21_mod_core(tmp, ten, &ost);
        if (ost.bits[0] < 5) {
          mul_XX(tmp, &tmp);
        }
        if (ost.bits[0] > 5) {
          mul_XX(tmp, &tmp);
          tmp = s21_add_core(tmp, one, &tmp);
        }
      }
    }

  }
  if (tmp.bits[3] == 0 && tmp.bits[4] == 0 && tmp.bits[5] == 0) {
    for (int i = 0; i < 3; i++) {
      result->bits[i] = tmp.bits[i];
    }

    result->bits[3] = tmp.bits[6];
    if (value.bits[6] >> 31) {
      result->bits[3] = result->bits[3] | (1 << 31);
    }
  } else {
    inf = 1;
  }

  return inf;
}

// vileplme

int s21_is_less(s21_decimal number1, s21_decimal number2) {
  int rv = 0;
  if (checkMinusZero(number1, number2)) {
    rv = 1;
  } else if (compare(number1, number2) < 0) {
    rv = 1;
  }
  return rv;
}

int checkMinusZero(s21_decimal number1, s21_decimal number2) {
  int rv = 0;
  int sign1 = 0;
  int sign2 = 0;
  get_sign(&sign1, number1);
  get_sign(&sign2, number2);
  if (number1.bits[0] == 0 && number1.bits[1] == 0 && number1.bits[2] == 0 &&
      sign1 == 1) {
    if (number2.bits[0] == 0 && number2.bits[1] == 0 && number2.bits[2] == 0 &&
        sign2 == 0) {
      rv = 1;
    }
  }
  return rv;
}

int s21_is_less_or_equal(s21_decimal number1, s21_decimal number2) {
  int rv = 0;
  if (compare(number1, number2) <= 0) {
    rv = 1;
  }
  return rv;
}

int s21_is_greater(s21_decimal number1, s21_decimal number2) {
  int rv = 0;
  if (compare(number1, number2) > 0) {
    rv = 1;
  }
  return rv;
}

int s21_is_greater_or_equal(s21_decimal number1, s21_decimal number2) {
  int rv = 0;
  if (compare(number1, number2) >= 0) {
    rv = 1;
  }
  return rv;
}

int s21_is_equal(s21_decimal number1, s21_decimal number2) {
  int rv = 0;
  if (compare(number1, number2) == 0) {
    rv = 1;
  }
  return rv;
}

int s21_is_not_equal(s21_decimal number1, s21_decimal number2) {
  int rv = 0;
  if (compare(number1, number2) != 0) {
    rv = 1;
  }
  return rv;
}

int compare(s21_decimal number1, s21_decimal number2) {
  int rv = 0;
  s21_big_decimal copy_number1 = {0};
  s21_big_decimal copy_number2 = {0};
  s21_big_decimal zero = {0};
  to_big_decimal(number1, &copy_number1);
  to_big_decimal(number2, &copy_number2);
  equalize(&copy_number1, &copy_number2);
  int sign1 = big_get_sign(copy_number1);
  int sign2 = big_get_sign(copy_number2);
  int answer = first_mismatch(copy_number1, copy_number2);
  if (memcmp(&copy_number1, &zero, 6 * sizeof(int)) == 0 &&
      memcmp(&copy_number2, &zero, 6 * sizeof(int)) == 0)
    rv = 0;
  else if (!sign1 && !sign2) {
    rv = answer;
  } else if (sign1 && sign2)
    rv = -answer;
  else if (sign1)
    rv = -1;
  else
    rv = 1;
  return rv;
}

void to_big_decimal(s21_decimal small, s21_big_decimal *big) {
  s21_big_decimal zero = {{0, 0, 0, 0, 0, 0, 0}};
  memcpy(big, &zero, 7 * sizeof(int));
  big->bits[6] = small.bits[3];
  big->bits[5] = small.bits[0];
  big->bits[4] = small.bits[1];
  big->bits[3] = small.bits[2];
}

void equalize(s21_big_decimal *number1, s21_big_decimal *number2) {
  s21_big_decimal temp = {{0, 0, 0, 0, 0, 10, 0}};
  int scale1 = big_get_scale(*number1);
  int scale2 = big_get_scale(*number2);
  while (scale1 > scale2) {
    mul_bit(*number2, temp, number2);
    scale2++;
  }
  while (scale1 < scale2) {
    mul_bit(*number1, temp, number1);
    scale1++;
  }
  big_set_scale(number1, scale1);
  big_set_scale(number2, scale2);
}

int mul_bit(s21_big_decimal number1, s21_big_decimal number2,
            s21_big_decimal *result) {
  int rv = 0;
  int length = BITS_IN_BIG - big_get_len(number2) - 1;
  s21_big_decimal zero = {0};
  for (int i = BITS_IN_BIG - 1; i >= length; i--) {
    if (big_get_bit(number2, i) == 1) {
      rv = sum_bit(zero, number1, &zero);
    }
    if (rv == 0 && i > 0) {
      rv = shift_left_compare(&number1);
    }
    if (rv != 0) {
      break;
    }
  }
  if (rv == 0) {
    memcpy(result, &zero, 6 * sizeof(int));
  }
  return rv;
}

int sum_bit(s21_big_decimal number1, s21_big_decimal number2,
            s21_big_decimal *result) {
  int rv = 0;
  int bit = 0;
  int save = 0;
  for (int i = BITS_IN_BIG - 1; i >= 0; i--) {
    bit = big_get_bit(number1, i) + big_get_bit(number2, i) + save;
    save = bit / 2;
    bit %= 2;
    big_set_bit(result, i, bit);
  }
  if (save > 0) {
    rv = 1;
  }
  return rv;
}

int big_get_bit(s21_big_decimal number, int index) {
  int rv = (number.bits[index / 32] >> (31 - (index % 32))) & 1;
  if (index >= BITS_IN_BIG) {
    rv = -1;
  }
  return rv;
}

void big_set_bit(s21_big_decimal *number, int index, int value) {
  if (index < BITS_IN_BIG) {
    if (value == 1) {
      number->bits[index / 32] |= 1 << (31 - (index % 32));
    } else {
      number->bits[index / 32] &= ~(1 << (31 - (index % 32)));
    }
  }
}

int big_get_len(s21_big_decimal number) {
  int length = 0;
  int bit = 0;
  while (length < BITS_IN_BIG) {
    bit = big_get_bit(number, length);
    length++;
    if (bit != 0) {
      length--;
      break;
    }
  }
  return BITS_IN_BIG - length;
}

int big_get_scale(s21_big_decimal number) {
  return (number.bits[6] >> 16) & 31;
}

void big_set_scale(s21_big_decimal *number, int scale) {
  int sign = big_get_sign(*number);
  number->bits[6] = scale << 16;
  big_set_sign(number, sign);
}

int shift_left_compare(s21_big_decimal *number) {
  int rv = 0;
  int bit1 = big_get_bit(*number, 0);
  int bit2 = big_get_bit(*number, 1 * BITS_IN_INT);
  int bit3 = big_get_bit(*number, 2 * BITS_IN_INT);
  int bit4 = big_get_bit(*number, 3 * BITS_IN_INT);
  int bit5 = big_get_bit(*number, 4 * BITS_IN_INT);
  int bit6 = big_get_bit(*number, 5 * BITS_IN_INT);
  if (bit1 == 0) {
    number->bits[0] <<= 1;
    number->bits[2] <<= 1;
    number->bits[1] <<= 1;
    number->bits[3] <<= 1;
    number->bits[4] <<= 1;
    number->bits[5] <<= 1;
    big_set_bit(number, 1 * BITS_IN_INT - 1, bit2);
    big_set_bit(number, 2 * BITS_IN_INT - 1, bit3);
    big_set_bit(number, 3 * BITS_IN_INT - 1, bit4);
    big_set_bit(number, 4 * BITS_IN_INT - 1, bit5);
    big_set_bit(number, 5 * BITS_IN_INT - 1, bit6);
  } else {
    rv = 1;
  }
  return rv;
}

int first_mismatch(s21_big_decimal number1, s21_big_decimal number2) {
  int rv = 0;
  for (int i = 0; i < BITS_IN_BIG; i++) {
    rv = big_get_bit(number1, i) - big_get_bit(number2, i);
    if (rv != 0) {
      break;
    }
  }
  return rv;
}

int big_get_sign(s21_big_decimal number) {
  int rv = 0;
  if (POSITIVE_MASK & number.bits[6]) {
    rv = 1;
  } else {
    rv = 0;
  }
  return rv;
}

void big_set_sign(s21_big_decimal *number, int value) {
  int sign = big_get_sign(*number);
  if ((sign == 0 && value != 0) || (sign != 0 && value == 0)) {
    number->bits[6] = number->bits[6] ^ (1 << 31);
  }
}

int s21_mod(s21_decimal value_1, s21_decimal value_2, s21_decimal *result) {
  int status = 0;
  s21_big_decimal result0 = {0};
  s21_big_decimal *a = malloc(sizeof(int) * 7);
  s21_big_decimal *b = malloc(sizeof(int) * 7);
  for (int i = 0; i < 7; i++) {
    a->bits[i] = 0;
    b->bits[i] = 0;
  }
  to_big_decimal(*result, &result0);
  scale_norming(value_1, value_2, a, b);
  s21_mod_core(*a, *b, &result0);
  if (s21_big_decimal_to_decimal(result0, result)) {
    status = 1;
  } else {
    status = 0;
  }
  if (value_1.bits[3] >> 31) {
    result->bits[3] = result->bits[3] | (1 << 31);
  }
  if (value_2.bits[0] == 0 && value_2.bits[1] == 0 && value_2.bits[2] == 0) {
    status = 3;
  }
  free(a);
  free(b);
  return status;
}
