// outgoing.h
// Douglas Davis, Matthew Epland
// holds information about pair of scattered particles

#ifndef OUTGOING_H
#define OUTGOING_H

#include <vector>
#include <utility>

class outgoing {

private:

  double _p1_mass;
  double _p2_mass;

  // pairs of x,y coordinates
  std::vector<std::pair<double,double> > _p1_trackPts;
  std::vector<std::pair<double,double> > _p2_trackPts;
  
public:

  outgoing() {}
  virtual ~outgoing() {}

  // setters
  void set_p1_mass(const double m) { _p1_mass = m; }
  void set_p2_mass(const double m) { _p2_mass = m; }

  void add_p1_trackPt(const double x, const double y) { _p1_trackPts.push_back(std::make_pair(x,y)); }
  void add_p2_trackPt(const double x, const double y) { _p2_trackPts.push_back(std::make_pair(x,y)); }

  // getters
  double p1_mass() const { return _p1_mass; }
  double p2_mass() const { return _p2_mass; }

  const std::vector<std::pair<double,double> > p1_trackPts() const { return _p1_trackPts; }
  const std::vector<std::pair<double,double> > p2_trackPts() const { return _p2_trackPts; }
  
};

#endif
