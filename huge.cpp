
// huge.cpp

#include "include.h"

static bool flag_test_dd = false;
static bool flag_test_qd = false;
static bool flag_verbose = false;

static bool print_result(bool result)
{
  if(result)
  {
    cout << "Test passed." << endl;
  }
  else
  {
    cout << "Test FAILED." << endl;
  }

  return result;
}

static void print_usage()
{
  cout << "qd_test [-h] [--dd] [--qd] [--all] [-v]" << endl;
  cout << "  Tests output of large numbers." << endl;
  cout << endl;
  cout << "-h --help    Prints this usage message." << endl;
  cout << "--dd         Perform tests with double-double types." << endl;
  cout << "--qd         Perform tests with quad-double types." << endl;
  cout << "               This is the default." << endl;
  cout << "--all        Perform both double-double and quad-double tests." << endl;
  cout << "-v --verbose Print detailed information for each test." << endl;
}

bool check(string str, string true_str)
{
  bool pass = (str == true_str);

  if( !pass)
  {
    cout << "     fail: " << str << endl;
    cout << "should be: " << true_str << endl;
  }
  else if(flag_verbose)
  {
    cout << "     pass: " << str << endl;
  }

  return pass;
}

template <class T> bool test_huge()
{
  bool pass = true;

  int digits = T::_ndigits - 1;

  T x = T::_pi* T("1.0e290");

  string pi_str = T::_pi.to_string(digits, 0, std::ios_base::fixed);

  if(flag_verbose)
  {
    cout << pi_str << endl;
  }

  for(int i = 0; i < 18; i++, x *= 10.0)
  {
    std::ostringstream os;

    os << pi_str << "e+" << (290 + i);

    pass &= check(x.to_string(digits), os.str() );
  }

  x = -T::_pi* T("1.0e290");

  pi_str = "-" + pi_str;

  for(int i = 0; i < 18; i++, x *= 10.0)
  {
    std::ostringstream os;

    os << pi_str << "e+" << (290 + i);

    pass &= check(x.to_string(digits), os.str() );
  }

  return pass;
}

template <class T> bool test_max(string true_str)
{
  bool pass = true;

  int digits = T::_ndigits - 1;

  pass &= check(T::_max.to_string(digits), true_str);

  pass &= check( (-T::_max).to_string(digits), "-" + true_str);

  return pass;
}

int main3(int argc, const char* argv[] )
{
  bool pass = true;

//unsigned int old_cw = 0;
//fpu_fix_start(&old_cw);

  for(int i = 1; i < argc; i++)
  {
    string arg(argv[i] );

    bool result = false;

    if(arg == "-h" || arg == "--help")
    {
      result = true;

      print_usage();

      exit(0);
    }

    if(arg == "--dd")
    {
      result = true;

      flag_test_dd = true;
    }

    if(arg == "--qd")
    {
      result = true;

      flag_test_qd = true;
    }

    if(arg == "--all")
    {
      result = true;

      flag_test_dd = flag_test_qd = true;
    }

    if(arg == "-v" || arg == "--verbose")
    {
      result = true;

      flag_verbose = true;
    }

    if( !result)
    {
      cerr << "Unknown flag '" << arg << "'." << endl;
    }
  }

  if( !flag_test_dd && !flag_test_qd)
  {
    flag_test_dd = true;
    flag_test_qd = true;
  }

  cout << "Testing output of huge numbers..." << endl;

  if(flag_test_dd)
  {
    cout << endl;
    cout << "Testing dd_real ..." << endl;

    pass &= test_huge<dd_real>();
    pass &= test_max<dd_real>("1.797693134862315807937289714053e+308");

    print_result(pass);
  }

  if(flag_test_qd)
  {
    cout << endl;
    cout << "Testing qd_real ..." << endl;

    pass &= test_huge<qd_real>();
    pass &= test_max<qd_real>("1.7976931348623158079372897140530286112296785259868571699620069e+308");

    print_result(pass);
  }

//fpu_fix_end(&old_cw);

  return pass ? 0 : 1;
}
