#include <R.h>
#include <Rinternals.h>

// Add two scalars.
SEXP add_(SEXP x_, SEXP y_) {
  double x = asReal(x_);
  double y = asReal(y_);

  double sum = x + y;

  return ScalarReal(sum);
}

// Multiply two scalars.
SEXP mult_(SEXP x_, SEXP y_) {
  double x = asReal(x_);
  double y = asReal(y_);

  double mult = x * y;

  return ScalarReal(mult);
}
