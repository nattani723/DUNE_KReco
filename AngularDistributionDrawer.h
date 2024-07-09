/*
 * Header file for Angular Distribution Drawer tool class
 * Plots angular distribution of hits and draw as histograms
 * Cheated PDG information also avialable
 */

#ifndef ANGULAR_DISTRIBUTION_DRAWER
#define ANGULAR_DISTRIBUTION_DRAWER 1

//#include "HitSplitAlg_module.h"
#include "ParticleDirectionFinder.h"
#include <cmath>
#include <TH2D.h>
#include <TFile.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>


namespace kaon_reconstruction
{


  class AngularDistributionDrawer
  {

  public:

    typedef std::map<int, std::map<int, std::map<int, double>>> AngularDistribution3DCheatPDGMap;

    AngularDistributionDrawer(const ParticleDirectionFinder& particle_direction_finder);

    void runDrawer(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DCheatPDGMap& angular_distribution_map_cheated_pdg, std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, std::map<recob::Hit,int>& hit_pdg_map, TCanvas* &c);


  private:

    void fill_angular_distribution_map_cheated_pdg(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DCheatPDGMap& angular_distribution_map_cheated_pdg, std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, std::map<recob::Hit,int>& hit_pdg_map) const;

    void draw_hist_angular_distribution_map_cheated_pdg(std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, TString outfile_name, TCanvas* &c) const;

    float m_theta_bin_size;
    float m_phi_bin_size;
    int m_num_bin_theta;
    int m_num_bin_phi;

    //TCanvas * c = new TCanvas("c", "c", 800, 600);
  };


  //------------------------


  AngularDistributionDrawer::AngularDistributionDrawer(const ParticleDirectionFinder& particle_direction_finder) :
    m_theta_bin_size(particle_direction_finder.get_theta_bin_size()),
    m_phi_bin_size(particle_direction_finder.get_phi_bin_size()),
    m_num_bin_theta(static_cast<int>(M_PI / particle_direction_finder.get_theta_bin_size())),
    m_num_bin_phi(static_cast<int>(2*M_PI / particle_direction_finder.get_phi_bin_size()))
  {
  }


  //------------------------------------------------------------------------------------------------------------------------------------------   

    void AngularDistributionDrawer::runDrawer(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DCheatPDGMap& angular_distribution_map_cheated_pdg, std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, std::map<recob::Hit,int>& hit_pdg_map, TCanvas* &c)
  {

    this->fill_angular_distribution_map_cheated_pdg(sp_list_roi, k_end, angular_distribution_map_cheated_pdg, h_angular_distribution_cheated_pdg, spacepointToHitMap, hit_pdg_map);

    this->draw_hist_angular_distribution_map_cheated_pdg(h_angular_distribution_cheated_pdg, "test.root", c);

  } 

  //------------------------------------------------------------------------------------------------------------------------------------------      

  void AngularDistributionDrawer::fill_angular_distribution_map_cheated_pdg(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DCheatPDGMap& angular_distribution_map_cheated_pdg, std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap, std::map<recob::Hit,int>& hit_pdg_map) const
  {

    //const TVector3 x_axis(1.,0.,0.);

    for(auto it_sp = sp_list_roi.begin(); it_sp != sp_list_roi.end(); ++it_sp){

      const TVector3 hit_position = (*it_sp)->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;

      double theta = displacement_vector.Theta();
      double phi = displacement_vector.Phi();
      int theta_factor = (int)(std::floor(theta / m_theta_bin_size));
      int phi_factor = (int)(std::floor(phi / m_phi_bin_size));

      // retrieve PDG info
      art::Ptr<recob::Hit> phit = spacepointToHitMap.at(*it_sp);
      int pdg = hit_pdg_map[*phit];

      angular_distribution_map_cheated_pdg[pdg][theta_factor][phi_factor] += TMath::Sin(theta);;
    }

    if(angular_distribution_map_cheated_pdg.empty()) return;

    for(const auto& angular_distribution_map : angular_distribution_map_cheated_pdg) {
      TH2D* h = new TH2D("", "", m_num_bin_theta, 0, M_PI, m_num_bin_phi, -M_PI, M_PI);
      for(const auto& theta_entry : angular_distribution_map.second) {
	for(const auto& phi_entry : theta_entry.second) {
	  h->Fill(m_theta_bin_size * theta_entry.first, m_phi_bin_size * phi_entry.first, phi_entry.second);
	}
      }
      h_angular_distribution_cheated_pdg[angular_distribution_map.first] = h;
    }

  }

  //------------------------------------------------------------------------------------------------------------------------------------------                       

  void AngularDistributionDrawer::draw_hist_angular_distribution_map_cheated_pdg(std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, TString outfile_name, TCanvas* &c) const
  {

    if(h_angular_distribution_cheated_pdg.empty()) return;

    TFile outfile(outfile_name, "update");
    auto h_stack = new THStack("h_stack", "");

    c->Clear();
    c->SetFillStyle(1001);

    TLegend * leg = new TLegend(0.1, 0.0, 0.9, 1.0, "");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.2);  
    leg->SetNColumns(2);

    const double Single_PadSplit = 0.85;
    TPad *p_plot = new TPad("p_plot","p_plot",0,0,1,Single_PadSplit);
    TPad *p_legend = new TPad("p_legend","p_legend",0,Single_PadSplit,1,1);
    p_legend->SetBottomMargin(0);
    p_legend->SetTopMargin(0.1);
    p_plot->SetTopMargin(0.01);

    std::set<int> added_pids;

    for(const auto& it_pdg : h_angular_distribution_cheated_pdg) {
      int pdg_code = it_pdg.first;
      TH2D* histogram = it_pdg.second;
      histogram->SetFillStyle(1001);

      int fill_color = kGray; // Default color

      switch (pdg_code) {
      case 321:
	fill_color = kBlue;
	if (added_pids.find(pdg_code) == added_pids.end()) {
	  leg->AddEntry(histogram, "K^{+}", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      case -13:
	fill_color = kCyan;
	if (added_pids.find(pdg_code) == added_pids.end()) {
	  leg->AddEntry(histogram, "#mu^{+}", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      case 211:
	fill_color = kMagenta;
	if (added_pids.find(pdg_code) == added_pids.end()) {
	  leg->AddEntry(histogram, "#pi^{+}", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      case -11:
      case 11:
	fill_color = kGreen + 2;
	if (added_pids.find(11) == added_pids.end() && added_pids.find(-11) == added_pids.end()) {
	  leg->AddEntry(histogram, "Shower", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      case 2212:
      case 2112:
	fill_color = kRed;
	if (added_pids.find(2212) == added_pids.end() && added_pids.find(2112) == added_pids.end() ) {
	  leg->AddEntry(histogram, "Nucleon", "f");
	  added_pids.insert(pdg_code);
	}
	break;
      default:
	if (added_pids.find(999) == added_pids.end()) {
	  leg->AddEntry(histogram, "Others", "f");
	  added_pids.insert(999);
	}
	break; // Use default color
      }

    //TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");
    //TCanvas * c = new TCanvas("c", "Canvas", 800, 600);
    //TCanvas * c("c", "Canvas", 800, 600);
    //c->SetFillStyle(1001);

      /*
    for(const auto& it_pdg : h_angular_distribution_cheated_pdg) {

      int pdg_code = it_pdg.first;
      TH2D* histogram = it_pdg.second;
      histogram->SetFillStyle(1001);

      int fill_color = kBlack; // Default color
      switch (pdg_code) {
      case 321: fill_color = kBlue; break;
      case -13: fill_color = kCyan; break;
      case 211: fill_color = kMagenta; break;
      case -11:
      case 11: fill_color = kGreen + 2; break;
      case 2212:
      case 2112: fill_color = kRed; break;
      default: break; // Use default color
      }
      */

      histogram->SetFillColor(fill_color);
      histogram->SetFillStyle(1001);
      h_stack->Add(histogram);
    }

    p_legend->Draw();
    p_legend->cd();
    leg->Draw();
    c->cd();
    p_plot->Draw();
    p_plot->cd();
    h_stack->Draw("hist lego3 0");

    h_stack->GetXaxis()->SetTitle("#Phi [rad]");
    h_stack->GetYaxis()->SetTitle("#Theta [rad]");
    h_stack->GetXaxis()->SetTitleOffset(1.5);
    h_stack->GetYaxis()->SetTitleOffset(1.5);

    c->Modified();
    c->Write();

    //h_stack->Draw("hist lego3 0");
    //leg->Draw("same");
    //c->Write();
  }


} //namespace kaon_reconstruction

#endif 
