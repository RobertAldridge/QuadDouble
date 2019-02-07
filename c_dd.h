
// c_dd.h

void c_dd_add(const double* a, const double* b, double* c);

void c_dd_add_d_dd(double a, const double* b, double* c);

void c_dd_add_dd_d(const double* a, double b, double* c);

void c_dd_sub(const double* a, const double* b, double* c);

void c_dd_sub_d_dd(double a, const double* b, double* c);

void c_dd_sub_dd_d(const double* a, double b, double* c);

void c_dd_mul(const double* a, const double* b, double* c);

void c_dd_mul_d_dd(double a, const double* b, double* c);

void c_dd_mul_dd_d(const double* a, double b, double* c);

void c_dd_div(const double* a, const double* b, double* c);

void c_dd_div_d_dd(double a, const double* b, double* c);

void c_dd_div_dd_d(const double* a, double b, double* c);

void c_dd_copy(const double* a, double* b);

void c_dd_copy_d(double a, double* b);

void c_dd_sqrt(const double* a, double* b);

void c_dd_sqr(const double* a, double* b);

void c_dd_abs(const double* a, double* b);

void c_dd_npwr(const double* a, int b, double* c);

void c_dd_nroot(const double* a, int b, double* c);

void c_dd_nint(const double* a, double* b);

void c_dd_aint(const double* a, double* b);

void c_dd_floor(const double* a, double* b);

void c_dd_ceil(const double* a, double* b);

void c_dd_exp(const double* a, double* b);

void c_dd_log(const double* a, double* b);

void c_dd_log10(const double* a, double* b);

void c_dd_sin(const double* a, double* b);

void c_dd_cos(const double* a, double* b);

void c_dd_tan(const double* a, double* b);

void c_dd_asin(const double* a, double* b);

void c_dd_acos(const double* a, double* b);

void c_dd_atan(const double* a, double* b);

void c_dd_atan2(const double* a, const double* b, double* c);

void c_dd_sinh(const double* a, double* b);

void c_dd_cosh(const double* a, double* b);

void c_dd_tanh(const double* a, double* b);

void c_dd_asinh(const double* a, double* b);

void c_dd_acosh(const double* a, double* b);

void c_dd_atanh(const double* a, double* b);

void c_dd_sincos(const double* a, double* s, double* c);

void c_dd_sincosh(const double* a, double* s, double* c);

void c_dd_read(const char* s, double* a);

void c_dd_swrite(const double* a, int precision, char* s, int len);

void c_dd_write(const double* a);

void c_dd_neg(const double* a, double* b);

void c_dd_rand(double* a);

void c_dd_comp(const double* a, const double* b, int* result);

void c_dd_comp_dd_d(const double* a, double b, int* result);

void c_dd_comp_d_dd(double a, const double* b, int* result);

void c_dd_pi(double* a);
