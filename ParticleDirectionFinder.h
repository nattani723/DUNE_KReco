/*

 * Header file for particle direction finder tool class.
 * Plots angular distribution of hits insider the region of interest,
 * and find peaks as candidates for daughter particle direction.

 */

#ifndef PARTICLE_DIRECTION_FINDER
#define PARTICLE_DIRECTION_FINDER

#include <bits/c++config.h>
#undef _GLIBCXX14_CONSTEXPR
#define _GLIBCXX14_CONSTEXPR

//#include "HitSplitAlg_module.h"
#include "TMath.h"
#include <algorithm>
#include <functional>

#include "Api/PandoraApi.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArContent.h"


//using namespace pandora;
//using namespace std;

namespace kaon_reconstruction
{

  class ParticleDirectionFinder
  {
    
  public:
    ParticleDirectionFinder();
    
    //Main functionality
    
  private:
    
    typdef std::vector<art::Ptr<recob::SpacePoint>> SPList;
    typdef std::map<int, std::map<int, double>> AngularDistributionMap3D;
    
  /*
   * @brief  Collect spacepoints inside the region of interest (ROI)
   *
   * @param sp_list: spacepoint list of K+ daughter candidate hits
   * @param k_end: vector of K+ track candidate's end
   * @param sp_list_roi: spacepoint list inside the ROI
   *
   */

  //void collect_sp_in_roi(const SPList& sp_list, const TVector3 k_end, SPList& sp_list_roi) const;
  void collect_sp_in_roi(const SPList& sp_list, const TVector3& k_end, double& region_of_interest, SPList& sp_list_roi) const;

  /*
   * @brief  Fill map of angular distribution for spacepoints inside the region of interest (ROI)
   *
   * @param sp_list_roi: spacepoint list inside the ROI
   * @param k_end: vector of K+ track candidate's end
   * @param angular_distribution_map: output map 
   *
   */

  //void FillAngularDistributionMap(const std::vector<art::Ptr<recob::SpacePoint>>& SPListROI, const TVector3 KEnd, AngularDistributionMap3D& AngularDistributionMap) const;
  void fill_angular_distribution_map(const std::vector<art::Ptr<recob::SpacePoint>>& sp_list_roi, const TVector3 k_end, AngularDistributionMap3D& angular_distribution_map) const;

  /*
   * @brief  Smooth out the angular distribution map
   *
   * @param angular_distribution_map: output map
   *
   */

  void smooth_angular_distribution_map(AngularDistributionMap3D& angular_distribution_map) const;

  /*
   * @brief  Obtain the directions from angular distribution peaks
   * 
   * @param angular_distribution_map: output map 
   * @param sort_peak_direction_map: vector of obtained peak directions sorted by peak height from high to low
   * 
   */

  //void retrieve_peak_directions(const angular_distribution_map_3d& angular_distribution_map, std::vector<TVector2>& peak_direction_vectors) const;
  void retrieve_peak_directions(const angular_distribution_map_3d& angular_distribution_map, std::map<double, TVector3, std::greater<>>& sort_peak_direction_map) const;

  /*
   * @brief  Remove some peak candidates depending on their height or opnening angles
   * 
   * @param sort_peak_direction_map: vector of obtained peak directions sorted by peak height from high to low
   * @param peak_direction_vector: vector of peak direction vector
   * 
   */

  void refine_peak_directions(const std::map<double, TVector3, std::greater<>>& sort_peak_direction_map, vector<TVector3> &peak_direction_vector) const;


  //int m_peak_search_region;
  float m_theta_bin_size;
  float m_phi_bin_size;
  int m_smoothing_window;
  int m_peak_search_window;
  double m_peak_open_angle;
  double m_min_peak_height;
  int m_max_num_peak;

};


} // end of namespace kaon_reconstruction

#endif // #ifndef PARTICLE_DIRECTION_FINDER




    //------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                                                   
    void HitSplitAlg::drawHistAngularDistributionMap3DCheat( std::map<int, TH2D*> &h_angular_distribution_pfparticle_cheat_3D,
						TString outfile_name,
						TCanvas* &c){
      
      TFile outfile(outfile_name, "update");
      auto hs = new THStack("hs", "");
      TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "");
      c->SetFillStyle(1001);
      
      if(h_angular_distribution_pfparticle_cheat_3D.size()==0) return;
      
      std::map<int, TH2D*>::iterator it;
      //cout << "h_angular_distribution_pfparticle_cheat_3D.size(): " << h_angular_distribution_pfparticle_cheat_3D.size() << endl;

      for(it=h_angular_distribution_pfparticle_cheat_3D.begin(); it!=h_angular_distribution_pfparticle_cheat_3D.end(); it++){

	it->second->SetFillStyle(1001);
	//cout << "it->first: " << it->first << ", it->second->GetEntries(): " << it->second->GetEntries() << endl;

	if(it->first == 321){
	  it->second->SetFillColor(kBlue);
	  //leg->AddEntry(it->second, "K+");
	}
	else if(it->first == -13){
	  it->second->SetFillColor(kCyan);
	  //leg->AddEntry(it->second, "mu+");
	}
	else if(it->first == 211){
	  it->second->SetFillColor(kMagenta);
	  //leg->AddEntry(it->second, "pi+");
	}
	else if(it->first == -11 || it->first == 11){
	  it->second->SetFillColor(kGreen+2);
	  //leg->AddEntry(it->second, "shower");
	}
	else if(it->first == 2212 || it->first == 2112){
	  it->second->SetFillColor(kRed);
	  //leg->AddEntry(it->second, "nucleon");
	}
	else{
	  it->second->SetFillColor(kBlack);
	}
	hs->Add(it->second);

      }

      /*
      for(auto const& h_org: h_angular_distribution_pfparticle_cheat_3D){

	TH2D h = h_org;	
	h.second.SetFillStyle(1001);
	
	if(h.first == 321) h.second.SetFillColor(kBlue);
	else if(h.first == -13) h.second.SetFillColor(kCyan);
	else if(h.first == 211) h.second.SetFillColor(kMagenta);
	else if(h.first == -11 || h.first == 11) h.second.SetFillColor(kGreen+2);
	else if(h.first == 2212 || 2112) h.second.SetFillColor(kRed);
	else h.second.SetFillColor(kBlack);
	hs->Add(h.second);
	
      }
      */
      hs->Draw("hist""lego3 0");
      //hs->Draw("hist nostack");                                                                                                                                                                                                                                                                                      
      leg->Draw();
      c->Write();
      
      delete hs;
      delete leg;
      
      return;
      
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------   
   
#endif
