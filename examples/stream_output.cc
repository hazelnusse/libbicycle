#include <iostream>
#include "bicycle.h"
#include "whipple.h"

using namespace bicycle;

int main() {
  Bicycle b;
  std::cout << b;
  Whipple w;
  b.set_parameters(w);
  std::cout << b;
}
