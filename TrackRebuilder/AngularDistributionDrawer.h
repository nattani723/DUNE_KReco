/*

 * Header file for Angular Distribution Drawer tool class
 * Plots angular distribution of hits and draw as histograms
 * Cheated PDG information also avialable

 */

#ifndef ANGULAR_DISTRIBUTION_DRAWER
#define ANGULAR_DISTRIBUTION_DRAWER 1

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
    AngularDistributionDrawer(const ParticleDirectionFinder& particle_direction_finder);

  private:

    typedef std::map<int, std::map<int, std::map<int, double>>> AngularDistribution3DCheatPDGMap;

    void fill_angular_distribution_map_cheated_pdg(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistribution3DCheatPDGMap& angular_distribution_map_cheated_pdg, std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, const std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>>& spacepointToHitMap) const;

    void draw_hist_angular_distribution_map_cheated_pdg(std::map<int, TH2D*> &h_angular_distribution_cheated_pdg, TString outfile_name) const;

    float m_theta_bin_size;
    float m_phi_bin_size;
    int m_num_bin_theta;
    int m_num_bin_phi;

  };


} //namespace kaon_reconstruction

#endif 


