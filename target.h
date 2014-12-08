#ifndef TARGET_H
#define TARGET_H

class target {

private:

  double _mass;
  
public:

  target() {}
  virtual ~target() {}

  void set_mass(const double m) { _mass = m; }
  
  double mass() const { return _mass; }

};

#endif
