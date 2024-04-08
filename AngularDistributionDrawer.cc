#include "AngularDistributionDrawer.h"

namespace kaon_reconstruction
{
  AngularDistributionDrawer::AngularDistributionDrawer(const ParticleDirectionFinder particle_direction_finder) :
    m_num_bin_theta(static_cast<int>(M_PI / particle_direction_finder.get_theta_bin_size())),
    m_num_bin_phi(static_cast<int>(2*M_PI / particle_direction_finder.get_phi_bin_size()))
  {
  }
  
  
  void AngularDistributionDrawer::fill_angular_distribution_map_cheated_pdg(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DCheatPDGMap& angular_distribution_map_cheated_pdg, std::map<int, TH2D*> &h_angular_distribution_cheated_pdg) const
  {

    const TVector3 x_axis(1.,0.,0.);

    for(auto it_sp = sp_list_roi.begin(); it_sp != sp_list_roi.end(); ++it_sp){

      const TVector3 hit_position = (*it_sp)->XYZ();
      const TVector3 displacement_vector = hit_position - k_end;

      double theta = distance_vector.Theta();
      double phi = distance_vector.Phi();
      int theta_factor = (int)(std::floor(theta / m_theta_bin_size));
      int phi_factor = (int)(std::floor(phi / m_phi_bin_size));

      // retrieve PDG info
      art::Ptr<recob::Hit> hit = fSpacePointsToHits.at(*spIter);
      truthHitMatcher(hit, particletmp);
      int pdg = particletmp->PdgCode();

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

  void AngularDistributionDrawer::draw_hist_angular_distribution_map_cheated_pdg(std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, TString outfile_name) const
  {

    if(h_angular_distribution_cheated_pdg.empty()) return;

    TFile outfile(outfile_name, "update");
    auto h_stack = new THStack("h_stack", "");
    TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");
    TCanvas * c("c", "Canvas", 800, 600);
    c->SetFillStyle(1001);
    
    for(const auto& it_pdg : h_angular_distribution_cheated_pdg) {
      
      int pdg_code = it_pdg.first;
      TH2D* histogram = it_pdg.second;
      histogram->SetFillStyle(1001);
      
      int fill_color = kBlack; // Default color
      switch (pdg_code) {
      case 321: fillColor = kBlue; break;
      case -13: fillColor = kCyan; break;
      case 211: fillColor = kMagenta; break;
      case -11:
      case 11: fillColor = kGreen + 2; break;
      case 2212:
      case 2112: fillColor = kRed; break;
      default: break; // Use default color
      }

      histogram->SetFillColor(fill_color);
      h_stack->Add(histogram);
      
    }
    
    h_stack->Draw("hist lego3 0");
    leg->Draw();
    c->Write();
 
  }


}// namespace kaon_reconstruction 
