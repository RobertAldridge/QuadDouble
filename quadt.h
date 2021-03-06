
// quadt.h

template<class T> class quadt
{
public:

  quadt(double eps);

  ~quadt();

  template<class F> int integrate_u(const F& f, double tol, T& result, double& err);

  template<class F> int integrate(const F& f, T a, T b, double tol, T& result, double& err);

private:

  int max_level;

  double initial_width;

  double final_width;

  int table_size;

  double eps;

  T* weights;

  T* points;

  template<class F> class UnitFunction
  {
    F f;

    T offset;

    T h;

  public:

    UnitFunction(const F& f, const T& a, const T& b) : f(f)
    {
      offset = 0.5 * (a + b);

      h = (b - a) * 0.5;
    }

    T operator()(T x) const
    {
      return f(offset + h * x) * h;
    }
  };

  void init_table();
};

template<class T> quadt<T>::quadt(double eps)
{
  max_level = 11;

  initial_width = 0.5;

  final_width = std::ldexp(initial_width, -max_level + 1);

  table_size = static_cast<int>(2.0 * 7.0 / final_width);

  this->eps = eps;

  init_table();
}

template<class T> quadt<T>::~quadt()
{
  delete[] weights;

  delete[] points;
}

template<class T> void quadt<T>::init_table()
{
  double h = initial_width * 2.0;

  double dt = 0.0;

  double t = 0.0;

  int i = 0;

  T sinh_t;

  T cosh_t;

  T sinh_s;

  T cosh_s;

  T x;

  T w;

  weights = new T[table_size];

  points = new T[table_size];

  for(int level = 1; level <= max_level; level++, h *= 0.5)
  {
    t = h * 0.5;

    dt = (level == 1) ? t : h;

    for(/*nop*/; /*nop*/; t += dt)
    {
      sincosh(T(t), sinh_t, cosh_t);

      sincosh(sinh_t, sinh_s, cosh_s);

      x = sinh_s / cosh_s;

/*w = (cosh_t / cosh_s) / cosh_s;*/

      w = (cosh_t / sqr(cosh_s) );

      if(x == 1.0 || w < eps)
      {
        weights[i++] = 0.0;

        break;
      }

      points[i] = x;

      weights[i] = w;

      i++;
    }
  }

}

template<class T> template<class F> int quadt<T>::integrate_u(const F& f, double tol, T& result, double& err)
{
  T r1;

  T r2;

  T r3;

  T s;

  T x;

  T w;

  int level = 0;

  double h = initial_width;

  bool conv = false;

  int i = 0;

  r1 = 0.0;

  r2 = 0.0;

  r3 = 0.0;

  s = f(T(0.0) );

  for(level = 1; level <= max_level; level++, h *= 0.5)
  {
    for(/*nop*/; /*nop*/; /*nop*/)
    {
      x = points[i];

      w = weights[i];

      i++;

      if(w == 0.0)
      {
        break;
      }

      s += w * (f(x) + f(-x) );
    }

    r1 = s * h;

    if(level > 2)
    {
      double e1 = 0.0;

      double e2 = 0.0;

      double d1 = 0.0;

      double d2 = 0.0;

      e1 = abs(to_double(r1 - r2) );

      if(e1 == 0.0)
      {
        err = eps;
      }
      else
      {
        e2 = abs(to_double(r1 - r3) );

        d1 = log(e1);

        d2 = log(e2);

        err = exp(d1 * d1 / d2);
      }

      std::cout << " level = " << level << std::endl;

      std::cout << " r = " << r1 << std::endl;

      std::cout << " err = " << err << std::endl;

      if(err < abs(r1) * tol)
      {
        conv = true;

        break;
      }
    }

    r2 = r1;

    r3 = r2;
  }

  if(level > max_level)
  {
    puts("Level exhausted.");
  }

  result = r1;

  if( !conv)
  {
    return -1;
  }

  return 0;
}

template<class T> template<class F> int quadt<T>::integrate(const F& f, T a, T b, double tol, T& result, double& err)
{
  if(a == -1.0 && b == 1.0)
  {
    return integrate_u(f, tol, result, err);
  }
  else
  {
    UnitFunction<F> unit_f(f, a, b);

    return integrate_u<UnitFunction<F> >(unit_f, tol, result, err);
  }
}
