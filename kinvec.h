#ifndef KINVEC_H
#define KINVEC_H

class kinvec {

private:

  double _x;
  double _y;
  double _vx;
  double _vy;

public:

  kinvec() {}
  virtual ~kinvec() {}

  void set_params(const double x,
		  const double y,
		  const double vx,
		  const double vy)
  {
    _x = x; _y = y; _vx = vx, _vy = vy;
  }

  const double x()  const { return _x;  }
  const double y()  const { return _y;  }
  const double vx() const { return _vx; }
  const double vy() const { return _vy; }
  
};

#endif
