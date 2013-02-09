#include <iostream>
#include "bicycle.h"
#include "whipple.h"

int main() {
  Bicycle b;
  std::cout << b;
  Whipple w;
  b.set_parameters_from_whipple(w);
  std::cout << b;
}
