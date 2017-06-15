#include <iostream>
#include <cstdlib>

extern"C" {
  double fact_ (int l);
}

int main()
{
  int l = 5;
  double y;
  y = fact_(l);
  std::cout << "Testing fact: l = " << l << "y = " << y;
  return (0);
  
}
