#include "incoming.h"
#include "target.h"
#include "outgoing.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{

  target   t;
  incoming inc;
  outgoing o;
  
  t.set_mass(5.5);
  inc.set_mvi(2.2,5.5,.1);

  std::cout << o.p1_mass() << std::endl;
  
  return 0;
  
}
