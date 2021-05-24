#include "TH2.h"
#include "TCanvas.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "JMTucker/Tools/interface/ExtValue.h"
#include "JMTucker/Tools/interface/PairwiseHistos.h"
#include "JMTucker/Tools/interface/Utilities.h"
#include "JMTucker/MFVNeutralinoFormats/interface/Event.h"
#include "JMTucker/MFVNeutralinoFormats/interface/VertexAux.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
class MFVVertexHistos : public edm::EDAnalyzer {
public:
	explicit MFVVertexHistos(const edm::ParameterSet&);
	void analyze(const edm::Event&, const edm::EventSetup&);

private:
	const edm::EDGetTokenT<MFVEvent> mevent_token;
	const edm::EDGetTokenT<double> weight_token;
	const edm::EDGetTokenT<MFVVertexAuxCollection> vertex_token;
	const int max_ntrackplots;
	const bool do_scatterplots;

	enum sv_index { sv_all, sv_num_indices };
	static const char* sv_index_names[sv_num_indices];

	void fill(TH1F** hs, const int, const double val, const double weight) const { hs[sv_all]->Fill(val, weight); }
	void fill(TH2F** hs, const int, const double val, const double val2, const double weight) const { hs[sv_all]->Fill(val, val2, weight); }

	Measurement1D miss_dist(const reco::Vertex& sv, const AlgebraicVector3& ref, const AlgebraicVector3& mom) {
		// miss distance is magnitude of (jet direction (= n) cross (tv - sv) ( = d))
		// to calculate uncertainty, use |n X d|^2 = (|n||d|)^2 - (n . d)^2
		AlgebraicVector3 n = ROOT::Math::Unit(mom);
		AlgebraicVector3 d(sv.x() - ref(0),
			sv.y() - ref(1),
			sv.z() - ref(2));
		AlgebraicVector3 n_cross_d = ROOT::Math::Cross(n, d);
		double n_dot_d = ROOT::Math::Dot(n, d);
		double val = ROOT::Math::Mag(n_cross_d);
		AlgebraicVector3 jac(2 * d(0) - 2 * n_dot_d * n(0),
			2 * d(1) - 2 * n_dot_d * n(1),
			2 * d(2) - 2 * n_dot_d * n(2));
		return Measurement1D(val, sqrt(ROOT::Math::Similarity(jac, sv.covariance())) / 1 / val); // modified err from 2->1 of sv
	}

	Measurement1D miss_dist_2D(const reco::Vertex & sv, const AlgebraicVector3 & ref, const AlgebraicVector3 & mom) {
		// miss distance is magnitude of (jet direction (= n) cross (tv - sv) ( = d))
		// to calculate uncertainty, use |n X d|^2 = (|n||d|)^2 - (n . d)^2
		AlgebraicVector3 n = ROOT::Math::Unit(mom);
		n(2) = 0.0;
		AlgebraicVector3 d(sv.x() - ref(0),
			sv.y() - ref(1),
			sv.z() - ref(2));
		d(2) = 0.0;
		AlgebraicVector3 n_cross_d = ROOT::Math::Cross(n, d);
		double n_dot_d = ROOT::Math::Dot(n, d);
		double val = ROOT::Math::Mag(n_cross_d);
		AlgebraicVector3 jac(2 * d(0) - 2 * n_dot_d * n(0),
			2 * d(1) - 2 * n_dot_d * n(1),
			2 * d(2) - 2 * n_dot_d * n(2));
		return Measurement1D(val, sqrt(ROOT::Math::Similarity(jac, sv.covariance())) / 1 / val); // modified err from 2->1 of sv and need sv to be modified for sig
	}
	
	TH1F* h_nsv;
	TH1F* h_output_shared_jet_or_not;
	

	
};

const char* MFVVertexHistos::sv_index_names[MFVVertexHistos::sv_num_indices] = { "all" };

MFVVertexHistos::MFVVertexHistos(const edm::ParameterSet & cfg)
	: mevent_token(consumes<MFVEvent>(cfg.getParameter<edm::InputTag>("mevent_src"))),
	weight_token(consumes<double>(cfg.getParameter<edm::InputTag>("weight_src"))),
	vertex_token(consumes<MFVVertexAuxCollection>(cfg.getParameter<edm::InputTag>("vertex_src"))),
	max_ntrackplots(cfg.getParameter<int>("max_ntrackplots")),
	do_scatterplots(cfg.getParameter<bool>("do_scatterplots"))
{
	edm::Service<TFileService> fs;

	h_nsv = fs->make<TH1F>("h_nsv", ";# of secondary vertices;arb. units", 15, 0, 15);
	h_output_shared_jet_or_not = fs->make<TH1F>("h_output_shared_jet_or_not", ";SV tracks share jet?", 2, 0, 2);

	}

void MFVVertexHistos::analyze(const edm::Event & event, const edm::EventSetup&) {

	unsigned run = event.id().run();
	unsigned lumi = event.luminosityBlock();
	unsigned long long evt = event.id().event();
	//if ( 42000 <= evt && evt <= 42115 ){
	
	edm::Handle<MFVEvent> mevent;
	event.getByToken(mevent_token, mevent);
	edm::Handle<double> weight;
	event.getByToken(weight_token, weight);
	const double w = *weight;
	
	const double bsx = mevent->bsx;
	const double bsy = mevent->bsy;
	const double bsz = mevent->bsz;
	const math::XYZPoint bs(bsx, bsy, bsz);
	const math::XYZPoint pv(mevent->pvx, mevent->pvy, mevent->pvz);

	
	edm::Handle<MFVVertexAuxCollection> auxes;
	event.getByToken(vertex_token, auxes);
	const int nsv = int(auxes->size());

	h_nsv->Fill(nsv,w);
	
	
	

	//////////////////////////////////////////////////////////////////////
	std::vector<std::vector<int> > sv_track_which_idx;
	std::vector<std::vector<int> > sv_track_which_jet;
	for (int isv = 0; isv < nsv; ++isv) {			//loop over vertices
		const MFVVertexAux& aux = auxes->at(isv);
		const int ntracks = aux.ntracks();
        
		
		//double phin = atan2(aux.y - bsy, aux.x - bsx);
		std::vector<int> track_which_idx;
		std::vector<int> track_which_jet;
		std::vector<double> absdeltaphi_sv_jets;
		for (int i = 0; i < ntracks; ++i) {		   //loop over tracks associated with a vertex
			double match_threshold = 1.3;
			int jet_index = 255;
			for (unsigned j = 0; j < mevent->jet_track_which_jet.size(); ++j) {		   //loop over all jets to take ones associated with each track
				double a = fabs(aux.track_pt(i) - fabs(mevent->jet_track_qpt[j])) + 1;
				double b = fabs(aux.track_eta[i] - mevent->jet_track_eta[j]) + 1;
				double c = fabs(aux.track_phi[i] - mevent->jet_track_phi[j]) + 1;
				if (a * b * c < match_threshold) {
					match_threshold = a * b * c;
					jet_index = mevent->jet_track_which_jet[j];
				}
			}
			if (jet_index != 255) {
				track_which_idx.push_back(i);                   // get one-to-one track_idx : track_idx
				track_which_jet.push_back((int)jet_index);	  // get one-to-one track_idx : jet_Index
			}
		}
		sv_track_which_jet.push_back(track_which_jet);
		sv_track_which_idx.push_back(track_which_idx);	   //[[0,1,2,3],[0,1]]


	

               }

	
	

	if (nsv >= 2) {
		const MFVVertexAux& sv0 = auxes->at(0);
		const MFVVertexAux& sv1 = auxes->at(1);
		double svdist2d = mag(sv0.x - sv1.x, sv0.y - sv1.y);
		double svdist3d = mag(sv0.x - sv1.x, sv0.y - sv1.y, sv0.z - sv1.z);
		double phi0 = atan2(sv0.y - bsy, sv0.x - bsx);
		double phi1 = atan2(sv1.y - bsy, sv1.x - bsx);
		double eta0 = atan2(sv0.y - bsy, sv0.z - bsz);
		double eta1 = atan2(sv1.y - bsy, sv1.z - bsz);

		/*
		std::vector<int> sv0_track_which_idx(int(sv0.ntracks()));
		int idx0 = 0;
		std::generate(sv0_track_which_idx.begin(), sv0_track_which_idx.end(), [&] { return idx0++; });
		std::vector<int> sv1_track_which_idx(int(sv1.ntracks()));
		int idx1 = 0;
		std::generate(sv1_track_which_idx.begin(), sv1_track_which_idx.end(), [&] { return idx1++; });
		*/	
		

		bool shared_jet = std::find_first_of(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()) != sv_track_which_jet[0].end();
		h_output_shared_jet_or_not->Fill(shared_jet);

		int nsharedjets = 0;
		std::vector<double> nsharedjet_phis;
                std::vector<double> nsharedjet_ets;
                std::vector<double> nsharedjet_etas;
				std::vector<int> nsharedjet_idxs;
		std::vector<std::vector<int>> sv_track_which_jet_copy = sv_track_which_jet;

		std::vector<int> nsharedjet_tracks_sv0;                                                                                                                                             
		std::vector<int> nsharedjet_tracks_sv1;
		std::vector<std::vector<int> >sv0_sharedjet_which_idx;                                                                                                                              
		std::vector<std::vector<int> >sv1_sharedjet_which_idx;

		std::vector<int> sv0_track_which_jet = sv_track_which_jet[0];                                                                                                    
		std::vector<int> sv0_track_which_idx = sv_track_which_idx[0];                                                                                                    
		std::vector<int> sv0_track_which_temp_idx;                                                                                                                                                                                                                                                                                                                              
		std::vector<int> sv1_track_which_jet = sv_track_which_jet[1];                                                                                                   
		std::vector<int> sv1_track_which_idx = sv_track_which_idx[1];                                                                                                   
		std::vector<int> sv1_track_which_temp_idx;

		while (std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end()) != sv_track_which_jet_copy[0].end()) {

			nsharedjets++;
			std::vector<int> sv0_track_which_idx_copy = sv0_track_which_idx;
			std::cout << " sv0's size: " << sv0_track_which_idx_copy.size() << std::endl;
			std::vector<int> sv1_track_which_idx_copy = sv1_track_which_idx;
			std::cout << " sv1's size: " << sv1_track_which_idx_copy.size() << std::endl;
			std::vector<int>::iterator it = std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end());
			int idx = std::distance(sv_track_which_jet_copy[0].begin(), it);
			int jet_index = sv_track_which_jet_copy[0].at(idx);
			/*
			std::vector<int>::iterator itr = std::find(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index);
			if (itr != sv_track_which_jet[0].cend()) {
				int j = std::distance(sv_track_which_jet[0].begin(), itr);
				nsharedjet_phis.push_back(mevent->jet_track_phi[j]);
                nsharedjet_ets.push_back(fabs(mevent->jet_track_qpt[j]));
                nsharedjet_etas.push_back(mevent->jet_track_eta[j]);
			}
			*/

			nsharedjet_phis.push_back(mevent->jet_phi[jet_index]);
			nsharedjet_ets.push_back(mevent->jet_energy[jet_index]);
			nsharedjet_etas.push_back(mevent->jet_eta[jet_index]);
			nsharedjet_idxs.push_back(jet_index);

			sv_track_which_jet_copy[0].erase(std::remove(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), jet_index), sv_track_which_jet_copy[0].end());
			sv_track_which_jet_copy[1].erase(std::remove(sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end(), jet_index), sv_track_which_jet_copy[1].end());
		
			// start counting shared tracks of sv0 for each shared jet                                                                                                                          
			nsharedjet_tracks_sv0.push_back(std::count(sv0_track_which_jet.begin(), sv0_track_which_jet.end(), jet_index));        
			std::cout << "# of sv0's shared-tracks = " << nsharedjet_tracks_sv0.back() << std::endl;
			std::multimap<int, size_t> sv0_m;                                                                                                                                                                       
			for (size_t k = 0; k < sv0_track_which_jet.size(); k++) if (sv0_track_which_jet[k] == jet_index) { sv0_m.insert({ sv0_track_which_jet[k], k }); }                                                                                                                                                                                                                                                                                                                                                                                                           
			for (auto it = sv0_m.begin(); it != sv0_m.end(); )                                                                                                                                  
			{                                                                                                                                                                                           
				auto p = sv0_m.equal_range(it->first);                                                                                                                                                                                                                                                                                                                                  
				while (p.first != p.second)                                                                                                                                                         
				{                                                                                                                                                                                           
					sv0_track_which_temp_idx.push_back(sv0_track_which_idx[p.first++->second]); 
					std::cout << "with jet index: " << jet_index << "idx is appended to a sv0 temp list: " << sv0_track_which_temp_idx.back() << std::endl;
					//[..] -> [..,1]
				}                                                                                                                                                                                   
				it = p.second;                                                                                                                                                                                                                                                                                                                                                  
			}                                                                                                                                                                                                                                                                                                                                                                       
			sv0_sharedjet_which_idx.push_back(sv0_track_which_temp_idx);                                                                                                                        
			for (size_t k = 0; k < sv0_track_which_temp_idx.size(); k++) {                                                                                                                              
				int track_index = sv0_track_which_temp_idx[k];                                                                                                                                      
				sv0_track_which_idx_copy.erase(std::remove(sv0_track_which_idx_copy.begin(), sv0_track_which_idx_copy.end(), track_index), sv0_track_which_idx_copy.end());                                                                                                                                                                                                                                                                                                                                                                                         
			}                                                                                                                                                                                                                                                                                                                                                                       
			                                                                                                                                                                                                                                                                                              
			sv0_track_which_temp_idx = {};                                                                                                                                                                                                                                                                                                                                          
			// start counting shared tracks of sv1 for each shared jet                                                                                                                          
			nsharedjet_tracks_sv1.push_back(std::count(sv1_track_which_jet.begin(), sv1_track_which_jet.end(), jet_index));   
			std::cout << "# of sv1's shared-tracks = " << nsharedjet_tracks_sv1.back() << std::endl;
			std::multimap<int, size_t> sv1_m;                                                                                                                                                   
			for (size_t k = 0; k < sv1_track_which_jet.size(); k++) if (sv1_track_which_jet[k] == jet_index) { sv1_m.insert({ sv1_track_which_jet[k], k }); }                                                                                                                                                                                                                                                                                                                                                                                                           
			for (auto it = sv1_m.begin(); it != sv1_m.end(); )                                                                                                                                  
			{                                                                                                                                                                                           
				auto p = sv1_m.equal_range(it->first);                                                                                                                                                                                                                                                                                                                                  
				while (p.first != p.second)                                                                                                                                                         
				{                                                                                                                                                                                           
					sv1_track_which_temp_idx.push_back(sv1_track_which_idx[p.first++->second]);  
					std::cout << "with jet index: " << jet_index << "idx is appended to a sv1 temp list: " << sv1_track_which_temp_idx.back() << std::endl;
				}                                                                                                                                                                                   
				it = p.second;                                                                                                                                                                                                                                                                                                                                                  
			}

			sv1_sharedjet_which_idx.push_back(sv1_track_which_temp_idx);                                                                                                                        
			for (size_t k = 0; k < sv1_track_which_temp_idx.size(); k++) { int track_index = sv1_track_which_temp_idx[k];                                                                                                                                      
			sv1_track_which_idx_copy.erase(std::remove(sv1_track_which_idx_copy.begin(), sv1_track_which_idx_copy.end(), track_index), sv1_track_which_idx_copy.end()); }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
			sv1_track_which_temp_idx = {}; 
		
		}


		if (shared_jet) {

			std::cout << "shared-jet event id: " << " run: " << run << " lumi: " << lumi << " event: " << evt << std::endl;
			std::cout << "the number of shared-jets is " << nsharedjets << std::endl;
			std::cout << "sv0's phi = " << phi0 << " and " << "sv1's phi = " << phi1 << std::endl;
			std::vector<int> sv0_sum_pt_track_which_idx = sv0_track_which_idx;
			std::vector<int> sv1_sum_pt_track_which_idx = sv1_track_which_idx;
			for (int i = 0; i < nsharedjets; ++i) {
				std::cout << i+1 << " shared jet's phi: " << nsharedjet_phis[i] << " shared jet's eta " << nsharedjet_etas[i] << " shared jet's eT " << nsharedjet_ets[i] << std::endl;
				double sum_pt_i_sv0 = 0;                                                                                                                                                            
				std::vector<int> sv0_i_sharedjet_which_idx = sv0_sharedjet_which_idx[i];       //consider [1,3,5,6]
				std::cout << "# of sv0's shared-tracks = " << sv0_i_sharedjet_which_idx.size() << std::endl;
				for (unsigned int j = 0; j < sv0_i_sharedjet_which_idx.size(); j++) { 		//used to be nsharedjet_tracks_sv0[i]
					int idx = sv0_i_sharedjet_which_idx[j]; 
					std::cout << "with jet index: " << nsharedjet_idxs[i] << "idx is appended to a sv0 temp list: " << idx << std::endl;
					sum_pt_i_sv0 = sum_pt_i_sv0 + sv0.track_pt(idx); 
					AlgebraicVector3 mom_tk(sv0.track_px[idx], sv0.track_py[idx], sv0.track_pz[idx]);
					AlgebraicVector3 ref_tk(sv0.track_vx[idx], sv0.track_vy[idx], sv0.track_vz[idx]);
					Measurement1D tkvtx_dist = miss_dist(sv0, ref_tk, mom_tk);
					std::cout << "  " << j + 1 << " shared track's phi: " << sv0.track_phi[idx] << " shared track's pt: " << sv0.track_pt(idx) << " shared track's sig_dxy" << tkvtx_dist.significance() << std::endl;
				}                                                                                                                                                                                   
				double sum_pt_i_sv1 = 0;                                                                                                                                                            
				std::vector<int> sv1_i_sharedjet_which_idx = sv1_sharedjet_which_idx[i];     
				std::cout << "# of sv1's shared-tracks = " << sv1_i_sharedjet_which_idx.size() << std::endl;
				for (unsigned int j = 0; j < sv1_i_sharedjet_which_idx.size(); j++) {	   //used to be nsharedjet_tracks_sv1[i]
					int idx = sv1_i_sharedjet_which_idx[j];  
					std::cout << "with jet index: " << nsharedjet_idxs[i] << "idx is appended to a sv1 temp list: " << idx << std::endl;
					sum_pt_i_sv1 = sum_pt_i_sv1 + sv1.track_pt(idx); 
					AlgebraicVector3 mom_tk(sv1.track_px[idx], sv1.track_py[idx], sv1.track_pz[idx]);
					AlgebraicVector3 ref_tk(sv1.track_vx[idx], sv1.track_vy[idx], sv1.track_vz[idx]);
					Measurement1D tkvtx_dist = miss_dist(sv1, ref_tk, mom_tk);
					std::cout << "  " << j + 1 << " shared track's phi: " << sv1.track_phi[idx] << " shared track's pt: " << sv1.track_pt(idx) << " shared track's sig_dxy" << tkvtx_dist.significance() << std::endl;

				}

				if (sum_pt_i_sv0 >= sum_pt_i_sv1) {
					std::cout << i+1 << ": sv0 is selected with the number of shared tracks of " << nsharedjet_tracks_sv0[i] << std::endl;
                    std::cout << i+1 << ": sv1 is non-selected with the number of shared tracks of " << nsharedjet_tracks_sv1[i] << std::endl;                                           
				}
				else {
					std::cout << i+1 << ": sv1 is selected with the number of shared tracks of " << nsharedjet_tracks_sv1[i] << std::endl;
                    std::cout << i+1 << ": sv0 is non-selected with the number of shared tracks of " << nsharedjet_tracks_sv0[i] << std::endl;                                       
				}
			}
			
		  
		}
		else {

          if (41000 <= evt && evt <= 42115){
			std::cout << "shared-jet event id: " << "run: " << run << "lumi: " << lumi << "event: " << evt << std::endl;
			std::cout << "the number of shared-jets is " << nsharedjets << std::endl;
			std::cout << "sv0's phi = " << phi0 << " and " << "sv1's phi = " << phi1 << std::endl;
			std::vector<int> sv0_sum_pt_track_which_idx = sv0_track_which_idx;
			std::vector<int> sv1_sum_pt_track_which_idx = sv1_track_which_idx;
          }
		
        }

		
	
    
     }
	
	
  
}

DEFINE_FWK_MODULE(MFVVertexHistos);

