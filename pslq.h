
// pslq.h

int minimum(int lhs, int rhs)
{
  return (lhs <= rhs)? lhs : rhs;
}

int maximum(int lhs, int rhs)
{
  return (lhs >= rhs)? lhs : rhs;
}

template<class T> void swap(T* lhs, T* rhs, int index)
{
  T temporary = lhs[index];

  lhs[index] = rhs[index];

  rhs[index] = temporary;
}

template<class T> T square(T lhs)
{
  return lhs * lhs;
}

template<class T> T** new_matrix(int nr_rows, int nr_cols, T diag = 0.0, T elem = 0.0)
{
  T** m = new T*[nr_rows];

  int i = 0;

  int j = 0;

  for(i = 0; i < nr_rows; i++)
  {
    m[i] = new T[nr_cols];
  }

  for(i = 0; i < nr_rows; i++)
  {
    for(j = 0; j < nr_cols; j++)
    {
      m[i][j] = (i == j) ? diag : elem;
    }
  }

  return m;
}

template<class T> T* new_vector(int n, T elem = 0.0)
{
  T* v = new T[n];

  for(int i = 0; i < n; i++)
  {
    v[i] = elem;
  }

  return v;
}

template<class T> void delete_matrix(T** m, int nr_rows)
{
  for(int i = 0; i < nr_rows; i++)
  {
    delete[] m[i];
  }

  delete[] m;
}

template<class T> void delete_vector(T* v)
{
  delete[] v;
}

static const double gam = 1.2;

template<class T> int pslq(const T* x, int n, T* r, double eps, int max_itr)
{
  T* s = new_vector<T>(n);

  T* y = new_vector<T>(n);

  T** a = new_matrix<T>(n, n, 1.0);

  T** b = new_matrix<T>(n, n, 1.0);

  T** h = new_matrix<T>(n, n - 1, 0.0);

  T t;

  double teps = 16.0 * eps;

  int i = 0;

  int j = 0;

  int k = 0;

  int err = 0;

  t = x[n-1] * x[n-1];

  s[n-1] = abs(x[n-1] );

  for(i = n-2; i >= 0; i--)
  {
    t += x[i] * x[i];

    s[i] = sqrt(t);
  }

  t = s[0];

  for(i = 0; i < n; i++)
  {
    y[i] = x[i] / t;
  }

  s[0] = 1.0;

  for(i = 1; i < n; i++)
  {
    s[i] /= t;
  }

  for(i = 0; i < n; i++)
  {
    for(j = 0; j <= minimum<T>(i, n - 2); j++)
    {
      h[i][j] = (i == j) ? s[j + 1] / s[j] : - y[i] * y[j] / (s[j] * s[j + 1] );
    }
  }

  for(i = 1; i < n; i++)
  {
    for(j = i - 1; j >= 0; j--)
    {
      t = nint(h[i][j] / h[j][j] );

      y[j] += t * y[i];

      for(k = 0; k <= j; k++)
      {
        h[i][k] -= t * h[j][k];
      }

      for(k = 0; k < n; k++)
      {
        a[i][k] -= t * a[j][k];

        b[k][j] += t * b[k][i];
      }
    }
  }

  int m = 0;

  int itr = 0;

  bool done = false;

  while( !done)
  {
    itr++;

    T m_val = -1.0;

    T g = gam;

    m = -1;

    for(i = 0; i < n - 1; i++, g *= gam)
    {
      t = abs(g * h[i][i] );

      if(t > m_val)
      {
        m_val = t;

        m = i;
      }
    }

    if(m < 0)
    {
      err = 1;

      break;
    }

    swap<T>(y/*[m]*/, y + 1/*[m + 1]*/, m);

    for(i = 0; i < n; i++)
    {
      swap<T>(a[m]/*[m][i]*/, a[m + 1]/*[m + 1][i]*/, i);
    }

    for(i = 0; i < n - 1; i++)
    {
      swap<T>(h[m]/*[m][i]*/, h[m + 1]/*[m + 1][i]*/, i);
    }

    for(i = 0; i < n; i++)
    {
      swap<T>(b[i]/*[i][m]*/, b[i] + 1/*[i][m + 1]*/, m);
    }

    if(m < n-2)
    {
      T t0;

      T t1;

      T t2;

      T t3;

      T t4;

      t0 = sqrt(square<T>(h[m][m] ) + square<T>(h[m][m + 1] ) );

      t1 = h[m][m] / t0;

      t2 = h[m][m + 1] / t0;

      for(i = m; i < n; i++)
      {
        t3 = h[i][m];

        t4 = h[i][m + 1];

        h[i][m] = t1 * t3 + t2 * t4;

        h[i][m + 1] = t1 * t4 - t2 * t3;
      }
    }

    for(i = m + 1; i < n; i++)
    {
      for(j = minimum<T>(i - 1, m + 1); j >= 0; j--)
      {
        t = nint(h[i][j] / h[j][j] );

        y[j] += t * y[i];

        for(k = 0; k <= j; k++)
        {
          h[i][k] -= t * h[j][k];
        }

        for(k = 0; k < n; k++)
        {
          a[i][k] -= t * a[j][k];

          b[k][j] += t * b[k][i];
        }
      }
    }

    m_val = -1.0e308;

    for(j = 0; j < n - 1; j++)
    {
      t = abs(h[j][j] );

      if(t > m_val)
      {
        m_val = t;
      }
    }

    for(i = 0; i < n; i++)
    {
      t = abs(y[i] );

      if(t < teps)
      {
        m = i;

        done = true;

        break;
      }
    }

    if(itr > max_itr)
    {
      done = true;

      err = -1;
    }
  }

  if(err == 0)
  {
    for(i = 0; i < n; i++)
    {
      r[i] = b[i][m];
    }
  }

  delete_matrix<T>(h, n);

  delete_matrix<T>(a, n);

  delete_matrix<T>(b, n);

  delete_vector<T>(y);

  delete_vector<T>(s);

  return err;
}
