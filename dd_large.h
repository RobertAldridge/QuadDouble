
// dd_large.h

#undef isnan

#undef isfinite

#undef isinf

#undef max

#undef min

struct dd_real
{
  double x[2];

  dd_real(double hi, double lo)
  {
    x[0] = hi;

    x[1] = lo;
  }

  dd_real()
  {
    x[0] = 0.0;

    x[1] = 0.0;
  }

  dd_real(double h)
  {
    x[0] = h;

    x[1] = 0.0;
  }

  dd_real(int h)
  {
    x[0] = (static_cast<double>(h) );

    x[1] = 0.0;
  }

  dd_real (const char* s);

  explicit dd_real(const double* d)
  {
    x[0] = d[0];

    x[1] = d[1];
  }

  static void error(const char* msg);

  double _hi() const
  {
    return x[0];
  }

  double _lo() const
  {
    return x[1];
  }

  static const dd_real _2pi;

  static const dd_real _pi;

  static const dd_real _3pi4;

  static const dd_real _pi2;

  static const dd_real _pi4;

  static const dd_real _e;

  static const dd_real _log2;

  static const dd_real _log10;

  static const dd_real _nan;

  static const dd_real _inf;

  static const double _eps;

  static const double _min_normalized;

  static const dd_real _max;

  static const dd_real _safe_max;

  static const int _ndigits;

  bool isnan() const
  {
    return std::isnan(x[0] ) || std::isnan(x[1] );
  }

  bool isfinite() const
  {
    return std::isfinite(x[0] );
  }

  bool isinf() const
  {
    return std::isinf(x[0] );
  }

  static dd_real add(double a, double b);

  static dd_real ieee_add(const dd_real& a, const dd_real& b);

  static dd_real sloppy_add(const dd_real& a, const dd_real& b);

  dd_real& operator+=(double a);

  dd_real& operator+=(const dd_real& a);

  static dd_real sub(double a, double b);

  dd_real& operator-=(double a);

  dd_real& operator-=(const dd_real& a);

  dd_real operator-() const;

  static dd_real mul(double a, double b);

  dd_real& operator*=(double a);

  dd_real& operator*=(const dd_real& a);

  static dd_real div(double a, double b);

  static dd_real sloppy_div(const dd_real& a, const dd_real& b);

  static dd_real accurate_div(const dd_real& a, const dd_real& b);

  dd_real& operator/=(double a);

  dd_real& operator/=(const dd_real& a);

  dd_real& operator=(double a);

  dd_real& operator=(const char* s);

  dd_real operator^(int n);

  static dd_real sqr(double d);

  static dd_real sqrt(double a);

  bool is_zero() const;

  bool is_one() const;

  bool is_positive() const;

  bool is_negative() const;

  static dd_real rand(void);

  void to_digits(char* s, int& expn, int precision = _ndigits) const;

  void write(char* s, int len, int precision = _ndigits, bool showpos = false, bool uppercase = false) const;

  std::string to_string(int precision = _ndigits, int width = 0, std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), bool showpos = false, bool uppercase = false, char fill = ' ') const;

  int read(const char* s, dd_real& a);

  void dump(const std::string& name = "", std::ostream& os = std::cerr) const;

  void dump_bits(const std::string& name = "", std::ostream& os = std::cerr) const;

  static dd_real debug_rand();
};

namespace std
{
  template<> class numeric_limits<dd_real> : public numeric_limits<double>
  {
  public:

    static double epsilon()
    {
      return dd_real::_eps;
    }

    static dd_real max()
    {
      return dd_real::_max;
    }

    static dd_real safe_max()
    {
      return dd_real::_safe_max;
    }

    static double min()
    {
      return dd_real::_min_normalized;
    }

    static const int digits = 104;

    static const int digits10 = 31;
  };
}

dd_real ddrand(void);

dd_real sqrt(const dd_real& a);

dd_real polyeval(const dd_real* c, int n, const dd_real& x);

dd_real polyroot(const dd_real* c, int n, const dd_real& x0, int max_iter = 32, double thresh = 0.0);

bool isnan(const dd_real& a)
{
  return a.isnan();
}

bool isfinite(const dd_real& a)
{
  return a.isfinite();
}

bool isinf(const dd_real& a)
{
  return a.isinf();
}

dd_real mul_pwr2(const dd_real& dd, double d);

dd_real operator+(const dd_real& a, double b);

dd_real operator+(double a, const dd_real& b);

dd_real operator+(const dd_real& a, const dd_real& b);

dd_real operator-(const dd_real& a, double b);

dd_real operator-(double a, const dd_real& b);

dd_real operator-(const dd_real& a, const dd_real& b);

dd_real operator*(const dd_real& a, double b);

dd_real operator*(double a, const dd_real& b);

dd_real operator*(const dd_real& a, const dd_real& b);

dd_real operator/(const dd_real& a, double b);

dd_real operator/(double a, const dd_real& b);

dd_real operator/(const dd_real& a, const dd_real& b);

dd_real inv(const dd_real& a);

dd_real rem(const dd_real& a, const dd_real& b);

dd_real drem(const dd_real& a, const dd_real& b);

dd_real divrem(const dd_real& a, const dd_real& b, dd_real& r);

dd_real pow(const dd_real& a, int n);

dd_real pow(const dd_real& a, const dd_real& b);

dd_real npwr(const dd_real& a, int n);

dd_real sqr(const dd_real& a);

dd_real sqrt(const dd_real& a);

dd_real nroot(const dd_real& a, int n);

bool operator==(const dd_real& a, double b);

bool operator==(double a, const dd_real& b);

bool operator==(const dd_real& a, const dd_real& b);

bool operator<=(const dd_real& a, double b);

bool operator<=(double a, const dd_real& b);

bool operator<=(const dd_real& a, const dd_real& b);

bool operator>=(const dd_real& a, double b);

bool operator>=(double a, const dd_real& b);

bool operator>=(const dd_real& a, const dd_real& b);

bool operator<(const dd_real& a, double b);

bool operator<(double a, const dd_real& b);

bool operator<(const dd_real& a, const dd_real& b);

bool operator>(const dd_real& a, double b);

bool operator>(double a, const dd_real& b);

bool operator>(const dd_real& a, const dd_real& b);

bool operator!=(const dd_real& a, double b);

bool operator!=(double a, const dd_real& b);

bool operator!=(const dd_real& a, const dd_real& b);

dd_real nint(const dd_real& a);

dd_real floor(const dd_real& a);

dd_real ceil(const dd_real& a);

dd_real aint(const dd_real& a);

dd_real ddrand(void);

double to_double(const dd_real& a);

int to_int(const dd_real& a);

dd_real exp(const dd_real& a);

dd_real ldexp(const dd_real& a, int exp);

dd_real log(const dd_real& a);

dd_real log10(const dd_real& a);

dd_real sin(const dd_real& a);

dd_real cos(const dd_real& a);

dd_real tan(const dd_real& a);

void sincos(const dd_real& a, dd_real& sin_a, dd_real& cos_a);

dd_real asin(const dd_real& a);

dd_real acos(const dd_real& a);

dd_real atan(const dd_real& a);

dd_real atan2(const dd_real& y, const dd_real& x);

dd_real sinh(const dd_real& a);

dd_real cosh(const dd_real& a);

dd_real tanh(const dd_real& a);

void sincosh(const dd_real& a, dd_real& sinh_a, dd_real& cosh_a);

dd_real asinh(const dd_real& a);

dd_real acosh(const dd_real& a);

dd_real atanh(const dd_real& a);

dd_real fabs(const dd_real& a);

dd_real abs(const dd_real& a);

dd_real fmod(const dd_real& a, const dd_real& b);

std::ostream& operator<<(std::ostream& s, const dd_real& a);

std::istream& operator>>(std::istream& s, dd_real& a);
