#ifndef INCOMING_H
#define INCOMING_H

class incoming {

private:

  double _mass;
  double _velocity;
  double _impact_param;
  
public:

  incoming() {}
  virtual ~incoming() {}

  void set_mvi(const double m, const double v, const double i)
  {
    _mass         = m;
    _velocity     = v;
    _impact_param = i;
  }

  double mass()         const { return _mass;         }
  double velocity()     const { return _velocity;     }
  double impact_param() const { return _impact_param; }
  
};

#endif
