/**
 * $Rev$ $Date$
 *
 * Copyright &copy; 2012-2017, Tech-X Corporation, Boulder, CO.
 * See LICENSE file (EclipseLicense.txt) for conditions of use.
 *
 * Determine whether the compiler knows std::abs<double>.
 */

#include <cmath>

int main(int argc, char** argv) {
  double a = 0;
  double b = std::abs(a);
  return 0;
}

