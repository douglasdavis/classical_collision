// Matt Epland, Douglas Davis
// PHY 761 Fall 2014 Duke University Physics Dept.
// C++ main for classical mechanics project

#include "particle.h"
#include <iostream>
#include <vector>
#include <cmath>

int main(int argc, char *argv[])
 {
  double projectile_mass   = 1.0;
  double projectile_radius = 1.0;
  double target_mass       = 1.0;
  double target_radius     = 1.0;

  particle p1(target_mass,target_radius);
  particle p2(projectile_mass,projectile_radius);

  kinvec initial_vec1;
  kinvec initial_vec2;
  
  initial_vec1.set_params(1,2,3,4);
  initial_vec2.set_params(2,5,6,8);

  p1.add_step(initial_vec1);
  p2.add_step(initial_vec2);
  
  double separation;      // |position1 - position2|
  double Fx, Fy, Fmag;    // force components, magnitude
  double K = 1;           // spring constant
  double Lambda = 1;      // force power law
  
  for ( auto i = 1; i < 10; ++i ) {

    // First we calculate the serpation of the particles determined from the
    // previous time step for calculations in the current step.
    // We also define force componenets for calculations in the current step.

    // sqrt[ (x1-x2)^2 + (y1-y2)^2 ]
    separation = std::sqrt(std::pow(p1.kinvecs().at(i-1).x() - p2.kinvecs().at(i-1).x(),2) + std::pow(p1.kinvecs().at(i-1).y() - p2.kinvecs().at(i-1).y(),2));
    // (x2 - x1 ) / separation
    Fx = ( p2.kinvecs().at(i-1).x() - p1.kinvecs().at(i-1).x() ) / separation;
    // (y2 - y1 ) / separation
    Fy = ( p2.kinvecs().at(i-1).y() - p1.kinvecs().at(i-1).y() ) / separation;
    // ( k * sqrt[ r1 + r2 - separation ]^(lambda) = k*displacement_from_equilbrium raised to lambda
    Fmag = K*std::pow(((p1.radius() + p2.radius()) - separation),Lambda);
    
    kinvec temp_vec1;
    kinvec temp_vec2;

    temp_vec1.set_params(p1.kinvecs().at(i-1).x()+1,
			 p1.kinvecs().at(i-1).y()+1,
			 p1.kinvecs().at(i-1).vx()+1,
			 p1.kinvecs().at(i-1).vy()+1);

    temp_vec2.set_params(p2.kinvecs().at(i-1).x()+i,
			 p2.kinvecs().at(i-1).y()+i,
			 p2.kinvecs().at(i-1).vx()+i,
			 p2.kinvecs().at(i-1).vy()+i);
    
    p1.add_step(temp_vec1);
    p2.add_step(temp_vec2);
    
    std::cout << separation << std::endl;
    
  }

  for ( auto const& step : p1.kinvecs() ) {
    std::cout << step.x() << " " << step.vx() << std::endl;
  }
  
  return 0;
}
