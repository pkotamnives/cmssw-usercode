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
	

	TH1F* h_ratio_ntracks_large_nsv2;

	
	TH1F* h_absdeltaphi_jet_shared_tracks_large_nsv2;
	TH1F* h_pt_minor_shared_tracks_large_nsv2;
	TH1F* h_pt_major_shared_tracks_large_nsv2;
	TH1F* h_dxy_err_minor_shared_tracks_large_nsv2;
	TH1F* h_dxy_err_major_shared_tracks_large_nsv2;
	TH1F* h_absdeltaphi_jet_major_sv_large_nsv2;
	
	TH1F* h_ratio_nsharedjets_large_nsv2;



	TH1F * h_lspdist2d_nsv2_shared_jets;
	TH1F * h_lspdist3d_nsv2_shared_jets;
	TH1F * h_absdeltaphi01_genlsp_nsv2_shared_jets;
	TH1F * h_lspdist2d_nsv2_no_shared_jets;
	TH1F * h_lspdist3d_nsv2_no_shared_jets;
	TH1F * h_absdeltaphi01_genlsp_nsv2_no_shared_jets;

	TH1F * h_nsharedjets_nsv2_shared_jets;



	TH1F * h_nsharedjets_large_nsv2_shared_jets;
	TH1F * h_svdist2d_large_absdeltaphi01_nsv2_shared_jets;
	TH1F * h_svdist3d_large_absdeltaphi01_nsv2_shared_jets;
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

	h_absdeltaphi_jet_shared_tracks_large_nsv2 = fs->make<TH1F>("h_absdeltaphi_jet_shared_tracks_large_nsv2", "nsv = 2, absdeltaphi01 > 0.5, track ratios >= 4; abs(delta(shared tracks, phi of a shared jet));arb. units", 316, 0, 3.16);
	h_ratio_ntracks_large_nsv2 = fs->make<TH1F>("h_ratio_ntracks_large_nsv2", "nsv = 2, absdeltaphi01 > 0.5;ratios of shared tracks (>=1);arb. units", 50, 0, 10);
	h_pt_minor_shared_tracks_large_nsv2 = fs->make<TH1F>("h_pt_minor_shared_tracks_large_nsv2", "nsv = 2, absdeltaphi01 > 0.5, track ratios >= 4; minor track p_{T} (GeV);arb. units", 100, 0, 100);
	h_pt_major_shared_tracks_large_nsv2 = fs->make<TH1F>("h_pt_major_shared_tracks_large_nsv2", "nsv = 2, absdeltaphi01 > 0.5, track ratios >= 4; major track p_{T} (GeV);arb. units", 100, 0, 100);
	h_dxy_err_minor_shared_tracks_large_nsv2 = fs->make<TH1F>("h_dxy_err_minor_shared_tracks_large_nsv2", "nsv = 2, absdeltaphi01 > 0.5, track ratios >= 4; minor track d_{xy} err (cm);arb. units", 100, 0, 0.1);
	h_dxy_err_major_shared_tracks_large_nsv2 = fs->make<TH1F>("h_dxy_err_major_shared_tracks_large_nsv2", "nsv = 2, absdeltaphi01 > 0.5, track ratios >= 4; major track d_{xy} err (cm);arb. units", 100, 0, 0.1);
	h_absdeltaphi_jet_major_sv_large_nsv2 = fs->make<TH1F>("h_absdeltaphi_jet_major_sv_large_nsv2", "nsv = 2, absdeltaphi01 > 0.5, track ratios >= 4; abs(delta(phi of major SV, phi of a shared jet));arb. units", 316, 0, 3.16);


    h_ratio_nsharedjets_large_nsv2 = fs->make<TH1F>("h_ratio_nsharedjets_large_nsv2", "nsv = 2, absdeltaphi01 > 0.5, semi-Fig2 && Fig2;ratios of shared jets (>=1);arb. units", 50, 0, 10);

	
	
	h_lspdist2d_nsv2_shared_jets = fs->make<TH1F>("h_lspdist2d_nsv2_shared_jets", "nsv = 2;dist2d(gen vtx #0, #1) (cm)", 200, 0, 2);
	h_lspdist3d_nsv2_shared_jets = fs->make<TH1F>("h_lspdist3d_nsv2_shared_jets", " nsv = 2;dist3d(gen vtx #0, #1) (cm)", 200, 0, 2);
	h_absdeltaphi01_genlsp_nsv2_shared_jets = fs->make<TH1F>("h_absdeltaphi01_genlsp_nsv2_shared_jets", "nsv = 2;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);
	h_lspdist2d_nsv2_no_shared_jets = fs->make<TH1F>("h_lspdist2d_nsv2_no_shared_jets", "nsv = 2;dist2d(gen vtx #0, #1) (cm)", 200, 0, 2);
	h_lspdist3d_nsv2_no_shared_jets = fs->make<TH1F>("h_lspdist3d_nsv2_no_shared_jets", " nsv = 2;dist3d(gen vtx #0, #1) (cm)", 200, 0, 2);
	h_absdeltaphi01_genlsp_nsv2_no_shared_jets = fs->make<TH1F>("h_absdeltaphi01_genlsp_nsv2_no_shared_jets", "nsv = 2;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);

    h_nsharedjets_nsv2_shared_jets = fs->make<TH1F>("h_nsharedjets_nsv2_shared_jets", "nsv = 2;# of shared jets;arb. units", 10, 0, 10);		
        	
	h_nsharedjets_large_nsv2_shared_jets = fs->make<TH1F>("h_nsharedjets_large_nsv2_shared_jets", "nsv = 2, absdeltaphi01 > 0.5;# of shared jets;arb. units", 10, 0, 10);
	h_svdist2d_large_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist2d_large_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 > 0.5 ;dist2d(sv #0, #1) (cm);arb. units", 500, 0, 1);
	h_svdist3d_large_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist3d_large_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 > 0.5 ;dist3d(sv #0, #1) (cm);arb. units", 500, 0, 1);
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
		sv_track_which_idx.push_back(track_which_idx);


	

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

              

		bool shared_jet = std::find_first_of(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()) != sv_track_which_jet[0].end();
		if (shared_jet) {

			int nsharedjets = 1;
			std::vector<int> nsharedjet_jet_index;
			std::vector<std::vector<int> > sv_track_which_jet_copy;
			sv_track_which_jet_copy = sv_track_which_jet;
			std::vector<int> nsharedjet_tracks_sv0;
			std::vector<int> nsharedjet_tracks_sv1;
			std::vector<std::vector<int> >sv0_sharedjet_which_idx;
			std::vector<std::vector<int> >sv1_sharedjet_which_idx;
			std::vector<int>::iterator it = std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end());
			int idx = std::distance(sv_track_which_jet_copy[0].begin(), it);
			int jet_index = sv_track_which_jet_copy[0].at(idx);
			nsharedjet_jet_index.push_back(jet_index);
			sv_track_which_jet_copy[0].erase(std::remove(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), jet_index), sv_track_which_jet_copy[0].end());
			sv_track_which_jet_copy[1].erase(std::remove(sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end(), jet_index), sv_track_which_jet_copy[1].end());
			nsharedjet_tracks_sv0.push_back(std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index));

			std::vector<int> sv0_track_which_jet = sv_track_which_jet[0];
			std::vector<int> sv0_track_which_idx = sv_track_which_idx[0];
			std::vector<int> sv0_track_which_temp_idx;
			std::multimap<int, size_t> sv0_m_nshj1;

			for (size_t k = 0; k < sv0_track_which_jet.size(); k++) if (sv0_track_which_jet[k] == jet_index) { sv0_m_nshj1.insert({ sv0_track_which_jet[k], k }); }

			for (auto it = sv0_m_nshj1.begin(); it != sv0_m_nshj1.end(); )
			{
				auto p = sv0_m_nshj1.equal_range(it->first);

				while (p.first != p.second)
				{
					sv0_track_which_temp_idx.push_back(sv0_track_which_idx[p.first++->second]);
				}
				it = p.second;

			}

			sv0_sharedjet_which_idx.push_back(sv0_track_which_temp_idx);                                                                                                                        sv0_track_which_temp_idx = {};

			nsharedjet_tracks_sv1.push_back(std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index));

			std::vector<int> sv1_track_which_jet = sv_track_which_jet[1];
			std::vector<int> sv1_track_which_idx = sv_track_which_idx[1];
			std::vector<int> sv1_track_which_temp_idx;
			std::multimap<int, size_t> sv1_m_nshj1;

			for (size_t k = 0; k < sv1_track_which_jet.size(); k++) if (sv1_track_which_jet[k] == jet_index) { sv1_m_nshj1.insert({ sv1_track_which_jet[k], k }); }

			for (auto it = sv1_m_nshj1.begin(); it != sv1_m_nshj1.end(); )
			{
				auto p = sv1_m_nshj1.equal_range(it->first);

				while (p.first != p.second)
				{
					sv1_track_which_temp_idx.push_back(sv1_track_which_idx[p.first++->second]);
				}
				it = p.second;

			}
			sv1_sharedjet_which_idx.push_back(sv1_track_which_temp_idx);                                                                                                                        
			sv1_track_which_temp_idx = {};
			// std::cout << "shared-jet #" << nsharedjets << " has shared-jet ntracks from sv#0 =" << std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index) << ", shared-jet ntracks from sv#1 =" << std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index) << std::endl;

			while (std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end()) != sv_track_which_jet_copy[0].end()) {
				nsharedjets++;
				it = std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end());
				idx = std::distance(sv_track_which_jet_copy[0].begin(), it);
				jet_index = sv_track_which_jet_copy[0].at(idx);
				nsharedjet_jet_index.push_back(jet_index);
				sv_track_which_jet_copy[0].erase(std::remove(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), jet_index), sv_track_which_jet_copy[0].end());
				sv_track_which_jet_copy[1].erase(std::remove(sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end(), jet_index), sv_track_which_jet_copy[1].end());
				nsharedjet_tracks_sv0.push_back(std::count(sv0_track_which_jet.begin(), sv0_track_which_jet.end(), jet_index));

				std::multimap<int, size_t> sv0_m;
				for (size_t k = 0; k < sv0_track_which_jet.size(); k++) if (sv0_track_which_jet[k] == jet_index) { sv0_m.insert({ sv0_track_which_jet[k], k }); }

				for (auto it = sv0_m.begin(); it != sv0_m.end(); )
				{
					auto p = sv0_m.equal_range(it->first);

					while (p.first != p.second)
					{
						sv0_track_which_temp_idx.push_back(sv0_track_which_idx[p.first++->second]);
					}
					it = p.second;

				}
				//std::cout << sv0_track_which_temp_idx.size() << "==" << nsharedjet_tracks_sv0.back() << std::endl;
				sv0_sharedjet_which_idx.push_back(sv0_track_which_temp_idx);                                                                                                                        
				sv0_track_which_temp_idx = {};

				nsharedjet_tracks_sv1.push_back(std::count(sv1_track_which_jet.begin(), sv1_track_which_jet.end(), jet_index));
				std::multimap<int, size_t> sv1_m;
				for (size_t k = 0; k < sv1_track_which_jet.size(); k++) if (sv1_track_which_jet[k] == jet_index) { sv1_m.insert({ sv1_track_which_jet[k], k }); }

				for (auto it = sv1_m.begin(); it != sv1_m.end(); )
				{
					auto p = sv1_m.equal_range(it->first);

					while (p.first != p.second)
					{
						sv1_track_which_temp_idx.push_back(sv1_track_which_idx[p.first++->second]);
					}
					it = p.second;

				}
				//std::cout << sv1_track_which_temp_idx.size() << "==" << nsharedjet_tracks_sv1.back() << std::endl;
				sv1_sharedjet_which_idx.push_back(sv1_track_which_temp_idx);                                                                                                                
				sv1_track_which_temp_idx = {};
				//  std::cout << "shared-jet #" << nsharedjets << " has shared-jet ntracks from sv#0 =" << std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index) << ", shared-jet ntracks from sv#1 =" << std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index) << std::endl;

			
			//std::cout << "# of shared jets are " << nsharedjets << std::endl;
			//std::cout << "sv#0 has ntracks="<< sv_track_which_jet[0].size() << ", sv#1 has ntracks="<< sv_track_which_jet[1].size() << std::endl;
			
		        }	
			
			if (nsv == 2) {	   //start nsv=2
				h_nsharedjets_nsv2_shared_jets->Fill(nsharedjets,w);
				std::vector<double> absdeltaphi_sv0_shared_jets;
				std::vector<double> absdeltaphi_sv1_shared_jets;
				std::vector<double> max_dphi_sv0_sv1;
                                
                                
				h_lspdist2d_nsv2_shared_jets->Fill(mevent->lspdist2d(), w);
				h_lspdist3d_nsv2_shared_jets->Fill(mevent->lspdist3d(), w);
				h_absdeltaphi01_genlsp_nsv2_shared_jets->Fill(std::abs(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1])), w);

				if (fabs(reco::deltaPhi(phi0, phi1)) > 0.5) {	//start no split vertex 

					h_nsharedjets_large_nsv2_shared_jets->Fill(nsharedjets, w);
					h_svdist2d_large_absdeltaphi01_nsv2_shared_jets->Fill(svdist2d, w);
					h_svdist3d_large_absdeltaphi01_nsv2_shared_jets->Fill(svdist3d, w);
					
					int nsharedjets_sv0 = 0;
					int nsharedjets_sv1 = 0;
                    for (int i = 0; i < nsharedjets; i++) {				//start nsharedjet loop
                                                
						int jet_index = nsharedjet_jet_index[i];

						double dphi_sv0_sharedjet = double(fabs(reco::deltaPhi(phi0, mevent->jet_phi[jet_index])));
						double dphi_sv1_sharedjet = double(fabs(reco::deltaPhi(phi1, mevent->jet_phi[jet_index])));


						double ratio_ntracks_nsv2;
						if (nsharedjet_tracks_sv0[i] > nsharedjet_tracks_sv1[i]) { ratio_ntracks_nsv2 = nsharedjet_tracks_sv0[i] / nsharedjet_tracks_sv1[i]; }
						else { ratio_ntracks_nsv2 = nsharedjet_tracks_sv1[i] / nsharedjet_tracks_sv0[i]; }
						h_ratio_ntracks_large_nsv2->Fill(ratio_ntracks_nsv2);

						if (ratio_ntracks_nsv2 >= 4) {		// study sample with all sh jets

						
							if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {		 // major vertex is SV1


								h_absdeltaphi_jet_major_sv_large_nsv2->Fill(dphi_sv1_sharedjet, w);


								std::vector<int> sv1_nsharedjets1_which_idx = sv1_sharedjet_which_idx[i];
								for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {		  // major tracks
									int track_idx = sv1_nsharedjets1_which_idx[j];
									
									double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx])));
									h_absdeltaphi_jet_shared_tracks_large_nsv2->Fill(absdelta_jet_track, w);
									h_pt_major_shared_tracks_large_nsv2->Fill(sv1.track_pt(track_idx), w);
									h_dxy_err_major_shared_tracks_large_nsv2->Fill(sv1.track_dxy_err(track_idx), w);
																		
								}
								
								std::vector<int> sv0_nsharedjets1_which_idx = sv0_sharedjet_which_idx[i];
								for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {		   // minor tracks
									int track_idx = sv0_nsharedjets1_which_idx[j];
									double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx])));                                                                     
									h_absdeltaphi_jet_shared_tracks_large_nsv2->Fill(absdelta_jet_track, w);
									h_pt_minor_shared_tracks_large_nsv2->Fill(sv0.track_pt(track_idx), w);
									h_dxy_err_minor_shared_tracks_large_nsv2->Fill(sv0.track_dxy_err(track_idx), w);
									
								}
								



							}
							
							if (nsharedjet_tracks_sv0[i] > nsharedjet_tracks_sv1[i]) {	  // major vertex is SV0


								h_absdeltaphi_jet_major_sv_large_nsv2->Fill(dphi_sv0_sharedjet, w);
								
								std::vector<int> sv0_nsharedjets1_which_idx = sv0_sharedjet_which_idx[i];
								for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {	// major tracks
									int track_idx = sv0_nsharedjets1_which_idx[j];
									double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx])));
									h_absdeltaphi_jet_shared_tracks_large_nsv2->Fill(absdelta_jet_track, w);
									h_pt_major_shared_tracks_large_nsv2->Fill(sv0.track_pt(track_idx), w);
									h_dxy_err_major_shared_tracks_large_nsv2->Fill(sv0.track_dxy_err(track_idx), w);
									
								}


								std::vector<int> sv1_nsharedjets1_which_idx = sv1_sharedjet_which_idx[i];
								for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {	 // minor tracks
									int track_idx = sv1_nsharedjets1_which_idx[j];
									double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx])));
									h_absdeltaphi_jet_shared_tracks_large_nsv2->Fill(absdelta_jet_track, w);
									h_pt_minor_shared_tracks_large_nsv2->Fill(sv1.track_pt(track_idx), w);
									h_dxy_err_minor_shared_tracks_large_nsv2->Fill(sv1.track_dxy_err(track_idx), w);
								}
								

							}

						}

						

					
					   }	//end nsharedjets loop
					   
					   

					                     
                       
					   double ratio_nsharedjets;
					   if (nsharedjets_sv0 != 0 && nsharedjets_sv1 != 0) {
						   if (nsharedjets_sv0 > nsharedjets_sv1) { ratio_nsharedjets = nsharedjets_sv0 / nsharedjets_sv1; }
						   else { ratio_nsharedjets = nsharedjets_sv1/ nsharedjets_sv0; }
						   h_ratio_nsharedjets_large_nsv2->Fill(ratio_nsharedjets);
					   }

					   else if (nsharedjets_sv0 != 0 && nsharedjets_sv1 == 0) {
						   h_ratio_nsharedjets_large_nsv2->Fill(double(0));
					   }

					   else if (nsharedjets_sv0 == 0 && nsharedjets_sv1 != 0) {
						   h_ratio_nsharedjets_large_nsv2->Fill(double(0));
					   }

											   


					

                                           

					}   //end no split vertex 
				
				}	//end nsv=2
		
                         


		

           }
		else {	   // no shared jets

		   h_lspdist2d_nsv2_no_shared_jets->Fill(mevent->lspdist2d(), w);
		   h_lspdist3d_nsv2_no_shared_jets->Fill(mevent->lspdist3d(), w);
		   h_absdeltaphi01_genlsp_nsv2_no_shared_jets->Fill(std::abs(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1])), w);
		   if (nsv == 2) {
			   h_nsharedjets_nsv2_shared_jets->Fill((int)0, w);
			   if ((reco::deltaPhi(phi0, phi1)) > 0.5) { h_nsharedjets_large_nsv2_shared_jets->Fill((int)0, w); }
		   }
		}
	}
 //  }
  
}

DEFINE_FWK_MODULE(MFVVertexHistos);
