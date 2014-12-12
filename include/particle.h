#ifndef PARTICLE_H
#define PARTICLE_H

#include "kinvec.h"
#include <vector>

class particle {

private:

  std::vector<kinvec> _kinvecs;
  double              _mass;
  double              _radius;

public:

  particle(const double m, const double r) : _mass(m), _radius(r) {}
  virtual ~particle() {}

  void add_step(const kinvec a) { _kinvecs.push_back(a); }

  double radius() const { return _radius; }
  double mass()   const { return _mass;   }

  const std::vector<kinvec> kinvecs() const { return _kinvecs; }
  
};

#endif
