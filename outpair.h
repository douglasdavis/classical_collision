#ifndef OUTPAIR_H
#define OUTPAIR_H

class outpair {

private:

  double _mass1;
  double _mass2;
  double _p1theta;
  double _p2theta;
  double _p1vx;
  double _p1vy;
  double _p2vx;
  double _p2vy;
  double _current_sep;

  double _cur_time;
  
public:

  outpair() {}
  virtual ~outpair() {}

  void set_initial();
  
  void step_time();

  double mass1()       const { return _mass1;       }
  double mass2()       const { return _mass2;       }
  double p1theta()     const { return _p1theta;     }
  double p2theta()     const { return _p2theta;     }
  double p1vx()        const { return _p1vx;        }
  double p1vy()        const { return _p1vy;        }
  double p2vy()        const { return _p2vx;        }
  double p2vx()        const { return _p2vy;        }
  double current_sep() const { return _current_sep; }

};

#endif
