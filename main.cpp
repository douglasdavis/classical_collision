// Matt Epland, Douglas Davis
// PHY 761 Fall 2014 Duke University Physics Dept.
// C++ main for classical mechanics project

#include "particle.h"
#include <iostream>
#include <vector>
#include <cmath>

#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEllipse.h"

#include "boost/program_options.hpp"

int main(int argc, char *argv[])
{
  namespace po = boost::program_options;
  po::options_description desc("options");
  desc.add_options()
    ("help,h","Print help message")
    ("p-radius", po::value<double>()->default_value(1.000),"projectile radius")
    ("t-radius", po::value<double>()->default_value(3.000),"target radius")
    ("p-mass",   po::value<double>()->default_value(1.000),"projectile mass")
    ("t-mass",   po::value<double>()->default_value(5.000),"target mass")
    ("p-vinit",  po::value<double>()->default_value(5.000),"projectile initial velocity")
    ("delta-t,t",po::value<double>()->default_value(0.001),"time step")
    ("impact,s", po::value<double>()->default_value(1.000),"impact parameter")
    ("spring,k", po::value<double>()->default_value(1.000),"spring constant")
    ("eta,e",    po::value<double>()->default_value(0.000),"energy disipation")
    ("lambda,l", po::value<double>()->default_value(1.000),"power law");

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);

  if ( vm.count("help") ) {
    std::cout << desc << std::endl;
    return 0;
  }
    
  double delta_t           = vm["delta-t"].as<double>();

  double projectile_mass   = vm["p-mass"].as<double>();
  double projectile_radius = vm["p-radius"].as<double>();
  double projectile_vx0    = vm["p-vinit"].as<double>(); // initial velocity of project (complete in x dir)
  
  double target_mass       = vm["t-mass"].as<double>();
  double target_radius     = vm["t-radius"].as<double>();
  double target_x0         = 1.00; // initial x position of the target
  double target_y0         = 0.00; // initial y position of the target
  
  double impact_paramter   = vm["impact"].as<double>();

  double K      = vm["spring"].as<double>();               // spring constant
  double Lambda = vm["lambda"].as<double>();               // force power law
  double eta    = vm["eta"].as<double>();               // energy disapation parameter (0 = none) (fraction loss/time step)
  
  particle p1(target_mass,target_radius);
  particle p2(projectile_mass,projectile_radius);

  // ____________________________________________________________
  
  kinvec initial_kinvec1;
  kinvec initial_kinvec2;

  double projectile_x0 = target_x0 - std::cos(std::asin(impact_paramter/(target_radius+projectile_radius)))*(target_radius+projectile_radius);
  double projectile_y0 = impact_paramter;
  
  initial_kinvec1.set_params(target_x0,target_y0,0,0); // target initial params (target as 0 velocity)
  initial_kinvec2.set_params(projectile_x0,projectile_y0,projectile_vx0,0.0); // projectile initial params (initial vy is 0)

  p1.add_step(initial_kinvec1);
  p2.add_step(initial_kinvec2);

  std::vector<double> sep_vector;
  std::vector<double> time_vector;
  std::vector<double> Fmag_vector;
  std::vector<double> CM_x_vector;
  std::vector<double> CM_y_vector;
  std::vector<double> CM_vx_vector;
  std::vector<double> CM_vy_vector;
  CM_x_vector.push_back((initial_kinvec1.x()*target_mass + initial_kinvec2.x()*projectile_mass)/(projectile_mass+target_mass));
  CM_y_vector.push_back((initial_kinvec1.y()*target_mass + initial_kinvec2.y()*projectile_mass)/(projectile_mass+target_mass));
  CM_vx_vector.push_back((initial_kinvec1.vx()*target_mass + initial_kinvec2.vx()*projectile_mass)/(projectile_mass+target_mass));
  CM_vy_vector.push_back((initial_kinvec1.vy()*target_mass + initial_kinvec2.vy()*projectile_mass)/(projectile_mass+target_mass));
  sep_vector.push_back(projectile_radius + target_radius);
  time_vector.push_back(0.0);

  double separation;      // |position1 - position2|
  double Fx, Fy, Fmag;    // force components, magnitude

  // we want to know what the minimum separation is and
  // and what iteration in the loop its at to draw the ellipses later
  double min_sep = 2*(projectile_radius + target_radius);
  double min_sep_time = 0;
  int    min_sep_i = 0;
  
  int i = 1;  
  do {
    // First we calculate the serpation of the particles determined from the
    // previous time step for calculations in the current step.
    // We also define force componenets for calculations in the current step.

    // sqrt[ (x1-x2)^2 + (y1-y2)^2 ]
    separation = std::sqrt(std::pow(p1.kinvecs().at(i-1).x() - p2.kinvecs().at(i-1).x(),2) + std::pow(p1.kinvecs().at(i-1).y() - p2.kinvecs().at(i-1).y(),2));

    if ( separation < min_sep ) {
      min_sep = separation;
      min_sep_i = i;
      min_sep_time = delta_t*i;
    }
    
    // (x2 - x1 ) / separation
    Fx = ( p2.kinvecs().at(i-1).x() - p1.kinvecs().at(i-1).x() ) / separation;
    // (y2 - y1 ) / separation
    Fy = ( p2.kinvecs().at(i-1).y() - p1.kinvecs().at(i-1).y() ) / separation;
    // ( k * sqrt[ r1 + r2 - separation ]^(lambda) = k*displacement_from_equilbrium raised to lambda
    Fmag = K*std::pow(((p1.radius() + p2.radius()) - separation),Lambda);
    
    kinvec new_kv1;
    kinvec new_kv2;

    // ___ target
    new_kv1.set_x(p1.kinvecs().at(i-1).x() + p1.kinvecs().at(i-1).vx()*delta_t + 0.5*Fmag*(-1)*Fx*delta_t*delta_t/target_mass);
    new_kv1.set_y(p1.kinvecs().at(i-1).y() + p1.kinvecs().at(i-1).vy()*delta_t + 0.5*Fmag*(-1)*Fy*delta_t*delta_t/target_mass);

    new_kv1.set_vx((1-eta)*p1.kinvecs().at(i-1).vx() + Fmag*(-1)*Fx*delta_t/target_mass);
    new_kv1.set_vy((1-eta)*p1.kinvecs().at(i-1).vy() + Fmag*(-1)*Fy*delta_t/target_mass);

    // ___ projectile
    new_kv2.set_x(p2.kinvecs().at(i-1).x() + p2.kinvecs().at(i-1).vx()*delta_t + 0.5*Fmag*Fx*delta_t*delta_t/projectile_mass);
    new_kv2.set_y(p2.kinvecs().at(i-1).y() + p2.kinvecs().at(i-1).vy()*delta_t + 0.5*Fmag*Fy*delta_t*delta_t/projectile_mass);

    new_kv2.set_vx((1-eta)*p2.kinvecs().at(i-1).vx() + Fmag*Fx*delta_t/projectile_mass);
    new_kv2.set_vy((1-eta)*p2.kinvecs().at(i-1).vy() + Fmag*Fy*delta_t/projectile_mass);


    // __ center of mass
    CM_x_vector.push_back((p1.kinvecs().at(i-1).x()*target_mass + p2.kinvecs().at(i-1).x()*projectile_mass)/(projectile_mass+target_mass));
    CM_y_vector.push_back((p1.kinvecs().at(i-1).y()*target_mass + p2.kinvecs().at(i-1).y()*projectile_mass)/(projectile_mass+target_mass));
   
    CM_vx_vector.push_back((p1.kinvecs().at(i-1).vx()*target_mass + p2.kinvecs().at(i-1).vx()*projectile_mass)/(projectile_mass+target_mass));
    CM_vy_vector.push_back((p1.kinvecs().at(i-1).vy()*target_mass + p2.kinvecs().at(i-1).vy()*projectile_mass)/(projectile_mass+target_mass));
 
    p1.add_step(new_kv1);
    p2.add_step(new_kv2);
    
    //    std::cout << separation << " " << projectile_radius + target_radius << std::endl;
    sep_vector.push_back(separation);
    time_vector.push_back(delta_t*i);
    Fmag_vector.push_back(Fmag);   
     
    i++;
  } while ( separation <= ( projectile_radius + target_radius ) ) ;

  TEllipse *min_sep_ellipse1 = new TEllipse(p1.kinvecs().at(min_sep_i).x(),p1.kinvecs().at(min_sep_i).y(),target_radius,target_radius);
  TEllipse *min_sep_ellipse2 = new TEllipse(p2.kinvecs().at(min_sep_i).x(),p2.kinvecs().at(min_sep_i).y(),projectile_radius,projectile_radius);
  TEllipse *end_ellipse1     = new TEllipse(p1.kinvecs().at(p1.kinvecs().size()-1).x(),p1.kinvecs().at(p1.kinvecs().size()-1).y(),target_radius,target_radius);
  TEllipse *end_ellipse2     = new TEllipse(p2.kinvecs().at(p2.kinvecs().size()-1).x(),p2.kinvecs().at(p2.kinvecs().size()-1).y(),projectile_radius,projectile_radius);
  min_sep_ellipse1->SetFillColorAlpha(kWhite,0);
  min_sep_ellipse2->SetFillColorAlpha(kWhite,0);
  min_sep_ellipse1->SetLineColor(kBlack);
  min_sep_ellipse2->SetLineColor(kBlue);
  end_ellipse1->SetLineColor(kBlack);
  end_ellipse2->SetLineColor(kBlue);
  
  /*
    std::cout << " ****\t******\t*****\t******\t** " << std::endl;
    std::cout << " * x \t * y \t * vx \t * vy \t * " << std::endl;
    for ( auto const& step : p1.kinvecs() ) {
    std::cout << " * " << step.x() << " \t * " << step.y() << " \t * " << step.vx() << " \t * " << step.vy() << " \t * " << std::endl;
    }

    TApplication tapp("tapp",&argc,argv);

  */

  // Seperation vs time graph

  TCanvas* c1 = new TCanvas("c1"," ",600,600);
  gPad->SetTopMargin(.13);
  gPad->SetBottomMargin(.13);
  gPad->SetRightMargin(.13);
  gPad->SetLeftMargin(.13);
  
  TGraph *sep_graph = new TGraph(sep_vector.size(),&time_vector[0],&sep_vector[0]);
  sep_graph->GetXaxis()->SetTitle("t");
  sep_graph->GetYaxis()->SetTitle("#||{#vec{x_{1}}-#vec{x_{2}}}");
  sep_graph->SetMarkerStyle(7);
  sep_graph->Draw("AP");
  sep_graph->SetTitle("Radial Seperation vs Time");

  c1->Print("out/sep_vs_t.pdf", "Portrait pdf");

  // x vs y graph
  TCanvas* c2 = new TCanvas("c2"," ",600,600);

  TGraph *xy1 = new TGraph();
  xy1->SetName("Target");
  TGraph *xy2 = new TGraph();
  xy2->SetName("Projectile");
  TGraph *CM_xy_graph = new TGraph(CM_x_vector.size(),&CM_x_vector[0],&CM_y_vector[0]);
  CM_xy_graph->SetName("Center of Mass");

  for ( auto j = 0; j < p1.kinvecs().size(); ++j ) {
    xy1->SetPoint(j,p1.kinvecs().at(j).x(),p1.kinvecs().at(j).y());
    xy2->SetPoint(j,p2.kinvecs().at(j).x(),p2.kinvecs().at(j).y());
  }
 

  TMultiGraph *mgxy = new TMultiGraph();
  xy1->SetMarkerStyle(7);
  xy1->SetMarkerColor(kBlack);
  xy2->SetMarkerStyle(7);
  xy2->SetMarkerColor(kBlue);
  CM_xy_graph->SetMarkerStyle(7);
  CM_xy_graph->SetMarkerColor(kGreen);
  mgxy->Add(xy1);
  mgxy->Add(xy2);
  mgxy->Add(CM_xy_graph);
  mgxy->Draw("AP");

  mgxy->SetTitle("Center of Mass Trajectories");
  mgxy->GetXaxis()->SetTitle("x");
  mgxy->GetYaxis()->SetTitle("y");
  mgxy->GetYaxis()->SetTitleOffset(1.25); // This plot looked bad it without it but not the others for some reason... ?

  TLegend* mgxyleg = new TLegend(0.15,0.75,0.4,0.9);
  mgxyleg->AddEntry(xy2,"Projectile","p");
  mgxyleg->AddEntry(xy1,"Target","p");
  mgxyleg->AddEntry(CM_xy_graph,"Center of Mass","p");
  mgxyleg->Draw("same");
  c2->Print("out/x_vs_y.pdf", "Portrait pdf");

  TCanvas *c2we = new TCanvas("c2we","c2we",600,600);
  mgxy->Draw("AP");
  mgxyleg->Draw("same");
  end_ellipse1->SetFillColorAlpha(kWhite,0);
  end_ellipse2->SetFillColorAlpha(kWhite,0);
  end_ellipse1->Draw("same");
  end_ellipse2->Draw("same");

  
  c2we->Modified(); c2we->Update();
  mgxy->GetXaxis()->SetLimits(-2.2*target_radius,2.2*target_radius);
  c2we->Modified(); c2we->Update();
  mgxy->GetYaxis()->SetRangeUser(-2.2*target_radius,2.2*target_radius);
  gPad->Modified();
  c2we->Print("out/x_vs_y_with_ellipse.pdf", "Portrait pdf");
  
  TCanvas *c2weMS = new TCanvas("c2weMS","c2we",600,600);
  mgxy->Draw("AP");
  mgxyleg->Draw("same");
  min_sep_ellipse1->Draw("same");
  min_sep_ellipse2->Draw("same");

  
  c2weMS->Modified(); c2weMS->Update();
  mgxy->GetXaxis()->SetLimits(-2.2*target_radius,2.2*target_radius);
  c2weMS->Modified(); c2weMS->Update();
  mgxy->GetYaxis()->SetRangeUser(-2.2*target_radius,2.2*target_radius);
  gPad->Modified();
  c2weMS->Print("out/x_vs_y_with_ellipse_MS.pdf", "Portrait pdf");



  // vx vs vy graph
  TCanvas* c3 = new TCanvas("c3"," ",400,350);

  TGraph *vxy1 = new TGraph();
  vxy1->SetName("Target");
  TGraph *vxy2 = new TGraph();
  vxy2->SetName("Projectile");
  TGraph *CM_vxy_graph = new TGraph(CM_vx_vector.size(),&CM_vx_vector[0],&CM_vy_vector[0]);
  CM_vxy_graph->SetName("Center of Mass");

  for ( auto j = 0; j < p1.kinvecs().size(); ++j ) {
    vxy1->SetPoint(j,p1.kinvecs().at(j).vx(),p1.kinvecs().at(j).vy());
    vxy2->SetPoint(j,p2.kinvecs().at(j).vx(),p2.kinvecs().at(j).vy());
  }
 

  TMultiGraph *mgvxy = new TMultiGraph();
  vxy1->SetMarkerStyle(7);
  vxy1->SetMarkerColor(kBlack);
  vxy2->SetMarkerStyle(7);
  vxy2->SetMarkerColor(kBlue);
  CM_vxy_graph->SetMarkerStyle(7);
  CM_vxy_graph->SetMarkerColor(kGreen);
  mgvxy->Add(vxy1);
  mgvxy->Add(vxy2);
  mgvxy->Add(CM_vxy_graph);
  mgvxy->Draw("AP");

  mgvxy->SetTitle("Center of Mass Velocities");
  mgvxy->GetXaxis()->SetTitle("v_{x}");
  mgvxy->GetYaxis()->SetTitle("v_{y}");

  TLegend* mgvxyleg = new TLegend(0.15,0.75,0.4,0.9);
  mgvxyleg->AddEntry(vxy2,"Projectile","p");
  mgvxyleg->AddEntry(vxy1,"Target","p");
  mgvxyleg->AddEntry(CM_vxy_graph,"Center of Mass","p");
  mgvxyleg->Draw();

  gPad->Modified();
  c3->Print("out/vx_vs_vy.pdf", "Portrait pdf");



  // vx vs t graph
  TCanvas* c4 = new TCanvas("c4"," ",600,600);

  TGraph *vx1 = new TGraph();
  vx1->SetName("Target");
  TGraph *vx2 = new TGraph();
  vx2->SetName("Projectile");
  TGraph *CM_vx_graph = new TGraph(CM_vx_vector.size(),&time_vector[0],&CM_vx_vector[0]);
  CM_vx_graph->SetName("Center of Mass");

  for ( auto j = 0; j < p1.kinvecs().size(); ++j ) {
    vx1->SetPoint(j,time_vector.at(j),p1.kinvecs().at(j).vx());
    vx2->SetPoint(j,time_vector.at(j),p2.kinvecs().at(j).vx());
  }
 
  TMultiGraph *mgvx = new TMultiGraph();
  vx1->SetMarkerStyle(7);
  vx1->SetMarkerColor(kBlack);
  vx2->SetMarkerStyle(7);
  vx2->SetMarkerColor(kBlue);
  CM_vx_graph->SetMarkerStyle(7);
  CM_vx_graph->SetMarkerColor(kGreen);
  mgvx->Add(vx1);
  mgvx->Add(vx2);
  mgvx->Add(CM_vx_graph);
  mgvx->Draw("AP");

  mgvx->SetTitle("X Velocities vs t");
  mgvx->GetXaxis()->SetTitle("t");
  mgvx->GetYaxis()->SetTitle("v_{x}");

  TLegend* mgvxleg = new TLegend(0.73,0.5-0.075,0.98,0.5+0.075);
  mgvxleg->AddEntry(vx2,"Projectile","p");
  mgvxleg->AddEntry(vx1,"Target","p");
  mgvxleg->AddEntry(CM_vx_graph,"Center of Mass","p");
  mgvxleg->Draw();

  gPad->Modified();
  c4->Print("out/vx_vs_t.pdf", "Portrait pdf");



  // vy vs t graph
  TCanvas* c5 = new TCanvas("c5"," ",600,600);

  TGraph *vy1 = new TGraph();
  vy1->SetName("Target");
  TGraph *vy2 = new TGraph();
  vy2->SetName("Projectile");
  TGraph *CM_vy_graph = new TGraph(CM_vy_vector.size(),&time_vector[0],&CM_vy_vector[0]);
  CM_vy_graph->SetName("Center of Mass");

  for ( auto j = 0; j < p1.kinvecs().size(); ++j ) {
    vy1->SetPoint(j,time_vector.at(j),p1.kinvecs().at(j).vy());
    vy2->SetPoint(j,time_vector.at(j),p2.kinvecs().at(j).vy());
  }
 
  TMultiGraph *mgvy = new TMultiGraph();
  vy1->SetMarkerStyle(7);
  vy1->SetMarkerColor(kBlack);
  vy2->SetMarkerStyle(7);
  vy2->SetMarkerColor(kBlue);
  CM_vy_graph->SetMarkerStyle(7);
  CM_vy_graph->SetMarkerColor(kGreen);
  mgvy->Add(vy1);
  mgvy->Add(vy2);
  mgvy->Add(CM_vy_graph);
  mgvy->Draw("AP");

  mgvy->SetTitle("X Velocities vs t");
  mgvy->GetXaxis()->SetTitle("t");
  mgvy->GetYaxis()->SetTitle("v_{y}");

  TLegend* mgvyleg = new TLegend(0.73,0.5-0.075,0.98,0.5+0.075);
  mgvyleg->AddEntry(vy2,"Projectile","p");
  mgvyleg->AddEntry(vy1,"Target","p");
  mgvyleg->AddEntry(CM_vy_graph,"Center of Mass","p");
  mgvyleg->Draw();

  gPad->Modified();
  c5->Print("out/vy_vs_t.pdf", "Portrait pdf");



  // Force vs time graph

  TCanvas* c6 = new TCanvas("c6"," ",600,600);

  TGraph *Fmag_graph = new TGraph(Fmag_vector.size(),&time_vector[0],&Fmag_vector[0]);
  Fmag_graph->GetXaxis()->SetTitle("t");
  Fmag_graph->GetYaxis()->SetTitle("#||{F}");
  Fmag_graph->SetMarkerStyle(7);
  Fmag_graph->Draw("AP");
  Fmag_graph->SetTitle("Force Magnitude vs Time");

  c6->Print("out/Fmag_vs_t.pdf", "Portrait pdf");

  //  tapp.Run();
  return 0;
}
