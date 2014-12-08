#include "incoming.h"
#include "target.h"
#include "outpair.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{

  target   t;
  incoming inc;
  outpair  op;
  
  t.set_mass(5.5);
  inc.set_mvi(2.2,5.5,.1);

  std::cout << op.p1vx() << std::endl;
  
  return 0;
  
}
