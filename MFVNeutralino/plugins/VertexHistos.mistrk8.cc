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
	
	TH1F* h_ratio_diff_pT_sum_sv_nsv2_no_shj;
	TH1F* h_ratio_diff_pT_sum_sv_nsv2_large_no_shj;
	TH1F* h_ratio_diff_pT_sum_major_minor_sv_nsv2_no_shj;
	TH1F* h_ratio_diff_pT_sum_major_minor_sv_nsv2_large_no_shj;
	
	TH1F* h_sv_njets_nsv1;
	TH1F* h_sv_njets_large_nsv2_no_shj;
	

	TH1F* h_ratio_ntracks_large_nsv2_shared_jets;
	TH1F* h_sv_njets_large_nsv2_shj;
	TH1F* h_sv_njets_no_less_sum_pt_shared_tracks_large_nsv2;
	TH1F* h_sv_njets_no_less_ratio_sum_pt_shared_tracks_large_nsv2;
	//efficiency plots for all shared jets 
	TH2F* h_2D_sv_tracks_large_nsv2;
	TH2F* h_2D_sv_tracks_shared_jets_large_nsv2;
	TH2F* h_2D_sv_tracks_no_shared_tracks_large_nsv2;
	TH2F* h_2D_sv_tracks_no_minor_shared_tracks_large_nsv2;

	//efficiency plots for all nsharedjets=1 
	
	TH2F* h_2D_sv_tracks_shared_jets_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2;
	TH2F* h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj1;	  //repetitive 
	TH2F* h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj2;
	TH2F* h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_more_nshj3;

	TH2F* h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj1;	  //repetitive 
	TH2F* h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj2;
	TH2F* h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_more_nshj3;

	TH2F* h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_large_nsv2;
	TH2F* h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_large_nsv2;
	TH2F* h_2D_sv_tracks_no_less_avg_pt_shared_tracks_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_large_nsv2;

	TH1F* h_diff_pT_avg_dPhi_shj_sv_large_nsv2;
	TH1F* h_diff_ratio_pT_avg_dPhi_shj_sv_large_nsv2;
	TH1F* h_diff_pT_sum_dPhi_shj_sv_large_nsv2;
	TH1F* h_diff_ratio_pT_sum_dPhi_shj_sv_large_nsv2;

	TH2F* h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2;

	TH2F* h_2D_poor_sv_tracks_shared_jets_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2;
	TH2F* h_2D_poor_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2;
	


	TH1F* h_less_dphi_more_shared_tracks_large_nsv2_all_fig2;
	TH1F* h_more_ratio_sum_major_sv_nsv2_all_fig2;
	TH1F* h_more_ratio_avg_major_sv_nsv2_all_fig2;
	TH1F* h_more_pT_avg_major_sv_nsv2_all_fig2;
	TH1F* h_more_pT_sig_avg_major_sv_nsv2_all_fig2;

	TH1F* h_diff_pT_avg_sv0_sv1_sv_nsv2;
	TH1F* h_diff_ratio_pT_avg_sv0_sv1_sv_nsv2;
	TH1F* h_diff_pT_sum_sv0_sv1_sv_nsv2;
	TH1F* h_diff_ratio_pT_sum_sv0_sv1_sv_nsv2;
	TH1F* h_ratio_diff_pT_sum_sv0_sv1_sv_nsv2;
	TH1F* h_ratio_diff_pT_avg_sv0_sv1_sv_nsv2;

	
	TH1F* h_vertex_chi2dof_absdeltaphi0_large_nsv2_fig2;
	TH1F* h_absdeltaphi_large_jet_shared_tracks_nsv2_fig2;
	TH1F* h_pt_poor_shared_tracks_large_nsv2_fig2;
	TH1F* h_pt_good_shared_tracks_large_nsv2_fig2;
	TH1F* h_dxy_err_poor_shared_tracks_large_nsv2_fig2;
	TH1F* h_dxy_err_good_shared_tracks_large_nsv2_fig2;
	TH1F* h_dxy_sig_poor_shared_tracks_large_nsv2_fig2;
	TH1F* h_dxy_sig_good_shared_tracks_large_nsv2_fig2;
	TH1F* h_ratio_ntracks_large_nsv2_shared_jets_fig2;
	TH1F* h_dPhi_sv_large_nsv2_fig2;																						 
	TH1F* h_dEta_sv_large_nsv2_fig2;
	TH1F* h_dPhi_good_tracks_its_sv_large_fig2;
	TH1F* h_dPhi_poor_tracks_its_sv_large_fig2;
	TH1F* h_dPhi_good_tracks_other_sv_large_fig2;
	TH1F* h_dPhi_poor_tracks_other_sv_large_fig2;
	TH2F* h_2D_dEta_dPhi_sv_large_nsv2_fig2;
	TH2F* h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_fig2;
	TH2F* h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_fig2;
	TH2F* h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_fig2;
	TH2F* h_2D_poor_dxy_err_good_dxy_err_large_sv_fig2;
	TH2F* h_2D_poor_dz_err_good_dz_err_large_sv_fig2;
	TH2F* h_2D_poor_dxy_sig_good_dxy_sig_large_sv_fig2;
	TH2F* h_2D_poor_dz_sig_good_dz_sig_large_sv_fig2;
	TH2F* h_2D_poor_pt_good_pt_large_sv_fig2;
	TH2F* h_2D_poor_avg_pt_good_avg_pt_large_sv_fig2;
	TH2F* h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_fig2;
	TH2F* h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_fig2;
	TH2F* h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_fig2;
	TH2F* h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_fig2;
	TH2F* h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_fig2;
	TH1F* h_pt_shared_jets_large_nsv2_fig2;
	TH1F* h_more_ratio_sum_major_sv_nsv2_fig2;
	TH1F* h_more_ratio_avg_major_sv_nsv2_fig2;
	TH1F* h_more_pT_avg_major_sv_nsv2_fig2;
	TH1F* h_more_pT_sig_avg_major_sv_nsv2_fig2;


	
	TH1F* h_vertex_chi2dof_absdeltaphi0_large_nsv2_semi_fig2;
	TH1F* h_absdeltaphi_large_jet_shared_tracks_nsv2_semi_fig2;
	TH1F* h_pt_poor_shared_tracks_large_nsv2_semi_fig2;
	TH1F* h_pt_good_shared_tracks_large_nsv2_semi_fig2;
	TH1F* h_dxy_err_poor_shared_tracks_large_nsv2_semi_fig2;
	TH1F* h_dxy_err_good_shared_tracks_large_nsv2_semi_fig2;
	TH1F* h_dxy_sig_poor_shared_tracks_large_nsv2_semi_fig2;
	TH1F* h_dxy_sig_good_shared_tracks_large_nsv2_semi_fig2;
	TH1F* h_ratio_ntracks_large_nsv2_shared_jets_semi_fig2;
	TH1F* h_dPhi_sv_large_nsv2_semi_fig2;
	TH1F* h_dEta_sv_large_nsv2_semi_fig2;
	TH1F* h_dPhi_good_tracks_its_sv_large_semi_fig2;
	TH1F* h_dPhi_poor_tracks_its_sv_large_semi_fig2;
	TH1F* h_dPhi_good_tracks_other_sv_large_semi_fig2;
	TH1F* h_dPhi_poor_tracks_other_sv_large_semi_fig2;
	TH2F* h_2D_dEta_dPhi_sv_large_nsv2_semi_fig2;
	TH2F* h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_semi_fig2;
	TH2F* h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_semi_fig2;
	TH2F* h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_semi_fig2;
	TH2F* h_2D_pt_dxy_err_absdeltaphi0_large_sv_semi_fig2;
	TH2F* h_2D_pt_dxy_err_absdeltaphi1_large_sv_semi_fig2;
	TH2F* h_2D_poor_dxy_err_good_dxy_err_large_sv_semi_fig2;
	TH2F* h_2D_poor_dz_err_good_dz_err_large_sv_semi_fig2;
	TH2F* h_2D_poor_dxy_sig_good_dxy_sig_large_sv_semi_fig2;
	TH2F* h_2D_poor_dz_sig_good_dz_sig_large_sv_semi_fig2;
	TH2F* h_2D_poor_pt_good_pt_large_sv_semi_fig2;
	TH2F* h_2D_poor_avg_pt_good_avg_pt_large_sv_semi_fig2;
	TH2F* h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_semi_fig2;
	TH2F* h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_semi_fig2;
	TH2F* h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_semi_fig2;
	TH2F* h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2;
	TH2F* h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2;
	TH1F* h_pt_shared_jets_large_nsv2_semi_fig2;
	TH1F* h_more_ratio_sum_major_sv_nsv2_semi_fig2;
	TH1F* h_more_ratio_avg_major_sv_nsv2_semi_fig2;
	TH1F* h_more_pT_avg_major_sv_nsv2_semi_fig2;
	TH1F* h_more_pT_sig_avg_major_sv_nsv2_semi_fig2;


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

	h_ratio_diff_pT_sum_sv_nsv2_no_shj = fs->make<TH1F>("h_ratio_diff_pT_sum_sv_nsv2_no_shj", "nsv = 2, no shared jets; #frac{sv0 shared sum pT - sv1 shared sum pT}{sv0 shared sum pT + sv1 shared sum pT}", 40, -2, 2);
	h_ratio_diff_pT_sum_sv_nsv2_large_no_shj = fs->make<TH1F>("h_ratio_diff_pT_sum_sv_nsv2_large_no_shj", "nsv = 2, absdeltaphi01 > 0.5, no shared jets; #frac{sv0 shared sum pT - sv1 shared sum pT}{sv0 shared sum pT + sv1 shared sum pT}", 40, -2, 2);
	h_ratio_diff_pT_sum_major_minor_sv_nsv2_no_shj = fs->make<TH1F>("h_ratio_diff_pT_sum_major_minor_sv_nsv2_no_shj", "nsv = 2, no shared jets; #frac{major shared sum pT - minor shared sum pT}{major shared sum pT + minor shared sum pT}", 40, -2, 2);
	h_ratio_diff_pT_sum_major_minor_sv_nsv2_large_no_shj = fs->make<TH1F>("h_ratio_diff_pT_sum_major_minor_sv_nsv2_large_no_shj", "nsv = 2, absdeltaphi01 > 0.5, no shared jets; #frac{major shared sum pT - minor shared sum pT}{major shared sum pT + minor shared sum pT}", 40, -2, 2);

	h_sv_njets_nsv1 = fs->make<TH1F>("h_sv_njets_large_nsv1", "nsv = 1; # of jets/SV;arb. units", 10, 0, 10);
	h_sv_njets_large_nsv2_no_shj = fs->make<TH1F>("h_sv_njets_large_nsv2_no_shj", "nsv = 2, absdPhi01 > 0.5, no shared jets; # of jets/SV;arb. units", 10, 0, 10);

	h_ratio_ntracks_large_nsv2_shared_jets = fs->make<TH1F>("h_ratio_ntracks_large_nsv2_shared_jets", "nsv = 2, absdPhi01 > 0.5, shared jets;ratios of shared tracks (>1);arb. units", 50, 0, 10);

	h_sv_njets_large_nsv2_shj = fs->make<TH1F>("h_sv_njets_large_nsv2_shj", "nsv = 2, absdPhi01 > 0.5, shared jets; # of jets/SV;arb. units", 10, 0, 10);
	h_sv_njets_no_less_sum_pt_shared_tracks_large_nsv2 = fs->make<TH1F>("h_sv_njets_no_less_sum_pt_shared_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5, no less_{sum p_{T}} shared tracks; # of jets/SV;arb. units", 10, 0, 10);
	h_sv_njets_no_less_ratio_sum_pt_shared_tracks_large_nsv2 = fs->make<TH1F>("h_sv_njets_no_less_ratio_sum_pt_shared_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5, no less_{ratio sum p_{T}} shared tracks; # of jets/SV;arb. units", 10, 0, 10);

	h_2D_sv_tracks_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5; # vtx's more tracks; # vtx's less tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_shared_jets_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_shared_jets_large_nsv2", "nsv = 2, absdPhi01 > 0.5, shared jets; # vtx's more tracks; # vtx's less tracks", 50, 0, 50, 50, 0, 50);
    h_2D_sv_tracks_no_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_shared_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5, no shared tracks; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_minor_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_minor_shared_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5, no minor tracks;  # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5, no less_{sum p_{T}} shared tracks; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj1 = fs->make<TH2F>("h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj1", "nsv = 2, absdPhi01 > 0.5, nshj=1, no less_{sum p_{T}} shared tracks; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj2", "nsv = 2, absdPhi01 > 0.5, nshj=2, no less_{sum p_{T}} shared tracks; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_more_nshj3 = fs->make<TH2F>("h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_more_nshj3", "nsv = 2, absdPhi01 > 0.5, nshj>=3, no less_{sum p_{T}} shared tracks; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);

	h_diff_pT_avg_dPhi_shj_sv_large_nsv2 = fs->make<TH1F>("h_diff_pT_avg_dPhi_shj_sv_large_nsv2", "nsv = 2, absdPhi01 > 0.5 ;delta(phi of a shared jet, phi of vtx w/ more_{avg p_{T}} );arb. units", 31, 0, 3.16);
	h_diff_ratio_pT_avg_dPhi_shj_sv_large_nsv2 = fs->make<TH1F>("h_diff_ratio_pT_avg_dPhi_shj_sv_large_nsv2", "nsv = 2, absdPhi01 > 0.5 ;delta(phi of a shared jet, phi of vtx w/ more_{ratio avg p_{T}} );arb. units", 31, 0, 3.16);
	h_diff_pT_sum_dPhi_shj_sv_large_nsv2 = fs->make<TH1F>("h_diff_pT_sum_dPhi_shj_sv_large_nsv2", "nsv = 2, absdPhi01 > 0.5 ;delta(phi of a shared jet, phi of vtx w/ more_{sum p_{T}} );arb. units", 31, 0, 3.16);
	h_diff_ratio_pT_sum_dPhi_shj_sv_large_nsv2 = fs->make<TH1F>("h_diff_ratio_pT_sum_dPhi_shj_sv_large_nsv2", "nsv = 2, absdPhi01 > 0.5 ;delta(phi of a shared jet, phi of vtx w/ more_{ratio sum p_{T}} );arb. units", 31, 0, 3.16);

	h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, no less_{sum p_{T}} shared tracks; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj1 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj1", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no less_{sum p_{T}} shared tracks; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj2", "poor nsv = 2, absdPhi01 > 0.5, nshj=2, no less_{sum p_{T}} shared tracks; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_more_nshj3 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_more_nshj3", "poor nsv = 2, absdPhi01 > 0.5, nshj>=3, no less_{sum p_{T}} shared tracks; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);

	h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5, no vtx's w/ less_{ratio sum pt} shared tracks ; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5, no vtx's w/ less_{ratio avg pt} shared tracks ; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_avg_pt_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_avg_pt_shared_tracks_large_nsv2", "nsv = 2, absdPhi01 > 0.5, no vtx's w/ less_{avg pt} shared tracks ; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, no vtx's w/ less_{ratio sum pt} shared tracks ; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, no vtx's w/ less_{ratio avg pt} shared tracks ; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, no vtx's w/ less_{avg pt} shared tracks ; # sv0 vtx's tracks; # sv1 vtx's tracks", 50, 0, 50, 50, 0, 50);

	h_2D_sv_tracks_shared_jets_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_shared_jets_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no shared tracks; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's shared tracks; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no vtx's w/ less_{ratio sum pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no vtx's w/ less_{ratio avg pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no vtx's w/ less_{avg pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no vtx's w/ less_{sum pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's w/ less_{ratio sum pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's w/ less_{ratio avg pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's w/ less_{avg pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2", "nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's w/ less_{sum pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	
	h_2D_poor_sv_tracks_shared_jets_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_shared_jets_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no shared tracks; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's shared tracks; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no vtx's w/ less_{ratio sum pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no vtx's w/ less_{ratio avg pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no vtx's w/ less_{avg pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no vtx's w/ less_{sum pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's w/ less_{ratio sum pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's w/ less_{ratio avg pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's w/ less_{avg pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);
	h_2D_poor_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2 = fs->make<TH2F>("h_2D_poor_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2", "poor nsv = 2, absdPhi01 > 0.5, nshj=1, no minor vtx's w/ less_{sum pt} shared tracks ; # major vtx's tracks; # minor vtx's tracks", 50, 0, 50, 50, 0, 50);

	//start to apply the >= n ratio cut

	h_less_dphi_more_shared_tracks_large_nsv2_all_fig2 = fs->make<TH1F>("h_less_dphi_more_shared_tracks_large_nsv2_all_fig2", "nsv = 2, absdeltaphi01 > 0.5, shared ratios > 1; min(dphi(SV,a shared jet)) most shared tracks ?", 2, 0, 2);
	h_more_ratio_sum_major_sv_nsv2_all_fig2 = fs->make<TH1F>("h_more_ratio_sum_major_sv_nsv2_all_fig2", "nsv = 2, absdeltaphi01 > 0.5, shared ratios > 1; more sum p_{T} ratio(shared/SV's tracks) is major SV ?", 2, 0, 2);
	h_more_ratio_avg_major_sv_nsv2_all_fig2 = fs->make<TH1F>("h_more_ratio_avg_major_sv_nsv2_all_fig2", "nsv = 2, absdeltaphi01 > 0.5, shared ratios > 1; more avg p_{T} ratio(shared/SV's tracks) is major SV ?", 2, 0, 2);
	h_more_pT_avg_major_sv_nsv2_all_fig2 = fs->make<TH1F>("h_more_pT_avg_major_sv_nsv2_all_fig2", "nsv = 2, absdeltaphi01 > 0.5, shared ratios > 1; more avg p_{T} shared tracks is major SV ?", 2, 0, 2);
	h_more_pT_sig_avg_major_sv_nsv2_all_fig2 = fs->make<TH1F>("h_more_pT_sig_avg_major_sv_nsv2_all_fig2", "nsv = 2, absdeltaphi01 > 0.5, shared ratios > 1; more avg p_{T} sig shared tracks is major SV ?", 2, 0, 2);
	
	h_diff_pT_avg_sv0_sv1_sv_nsv2 = fs->make<TH1F>("h_diff_pT_avg_sv0_sv1_sv_nsv2", "nsv = 2, absdeltaphi01 > 0.5; sv0 shared avg p_{T} - sv1 shared avg p_{T}", 200, -100, 100);
	h_diff_ratio_pT_avg_sv0_sv1_sv_nsv2 = fs->make<TH1F>("h_diff_ratio_pT_avg_sv0_sv1_sv_nsv2", "nsv = 2, absdeltaphi01 > 0.5; sv0 ratio - sv1 ratio (shared avg p_{T}/SV all avg p_{T})", 40, -2, 2);
	h_diff_pT_sum_sv0_sv1_sv_nsv2 = fs->make<TH1F>("h_diff_pT_sum_sv0_sv1_sv_nsv2", "nsv = 2, absdeltaphi01 > 0.5; sv0 shared sum p_{T} - sv1 shared sum p_{T}", 200, -1000, 1000);
	h_diff_ratio_pT_sum_sv0_sv1_sv_nsv2 = fs->make<TH1F>("h_diff_ratio_pT_sum_sv0_sv1_sv_nsv2", "nsv = 2, absdeltaphi01 > 0.5; sv0 ratio - sv1 ratio (shared sum p_{T}/SV all sum p_{T})", 40, -2, 2);
	h_ratio_diff_pT_sum_sv0_sv1_sv_nsv2 = fs->make<TH1F>("h_ratio_diff_pT_sum_sv0_sv1_sv_nsv2", "nsv = 2, absdeltaphi01 > 0.5; #frac{sv0 shared sum pT - sv1 shared sum pT}{sv0 shared sum pT + sv1 shared sum pT}", 40, -2, 2);
	h_ratio_diff_pT_avg_sv0_sv1_sv_nsv2 = fs->make<TH1F>("h_ratio_diff_pT_avg_sv0_sv1_sv_nsv2", "nsv = 2, absdeltaphi01 > 0.5; #frac{sv0 shared avg pT - sv1 shared avg pT}{sv0 shared avg pT + sv1 shared avg pT}", 40, -2, 2);


	h_dPhi_good_tracks_its_sv_large_fig2 = fs->make<TH1F>("h_dPhi_good_tracks_its_sv_large_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;delta(phi of good tracks, phi of major vertex);arb. units", 31, 0, 3.16);
	h_dPhi_poor_tracks_its_sv_large_fig2 = fs->make<TH1F>("h_dPhi_poor_tracks_its_sv_large_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;delta(phi of poor tracks, phi of minor vertex);arb. units", 31, 0, 3.16);
	h_dPhi_good_tracks_other_sv_large_fig2 = fs->make<TH1F>("h_dPhi_good_tracks_other_sv_large_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;delta(phi of good tracks, phi of minor vertex);arb. units", 31, 0, 3.16);
	h_dPhi_poor_tracks_other_sv_large_fig2 = fs->make<TH1F>("h_dPhi_poor_tracks_other_sv_large_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;delta(phi of poor tracks, phi of major vertex);arb. units", 31, 0, 3.16);	
	h_vertex_chi2dof_absdeltaphi0_large_nsv2_fig2 = fs->make<TH1F>("h_vertex_chi2dof_absdeltaphi0_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; sv opposite to the shared jet's tracks #chi^ { 2 } / dof ;arb. units", 100, 0, 10);
	h_absdeltaphi_large_jet_shared_tracks_nsv2_fig2 = fs->make<TH1F>("h_absdeltaphi_large_jet_shared_tracks_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; abs(delta(shared tracks, phi of the only-one shared jet));arb. units", 31, 0, 3.16);
	h_ratio_ntracks_large_nsv2_shared_jets_fig2 = fs->make<TH1F>("h_ratio_ntracks_large_nsv2_shared_jets_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;ratios of shared tracks (>=1);arb. units", 50, 0, 10);
	h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_fig2 = fs->make<TH2F>("h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;miss dist sig(1 poor trk, minor vertex);avg miss dist sig(>=4 good trk, major vertex)", 200, 0, 20, 200, 0, 20);
	h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_fig2 = fs->make<TH2F>("h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;miss dist sig(1 poor trk, major vertex);avg miss dist sig(>=4 good trk, minor vertex)", 100, 0, 100, 100, 0, 100);
	h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_fig2 = fs->make<TH2F>("h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;miss dist sig(1 poor trk, minor vertex);miss dist sig(1 poor trk, major vertex)", 200, 0, 20, 100, 0, 100);
	h_2D_poor_dxy_err_good_dxy_err_large_sv_fig2 = fs->make<TH2F>("h_2D_poor_dxy_err_good_dxy_err_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;poor track d_{xy} err (cm);good track d_{xy} err (cm)", 50, 0, 0.05, 50, 0, 0.05);
	h_2D_poor_dz_err_good_dz_err_large_sv_fig2 = fs->make<TH2F>("h_2D_poor_dz_err_good_dz_err_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;poor track d_{z} err (cm);good track d_{z} err (cm)", 50, 0, 0.05, 50, 0, 0.05);
	h_2D_poor_dxy_sig_good_dxy_sig_large_sv_fig2 = fs->make<TH2F>("h_2D_poor_dxy_sig_good_dxy_sig_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;poor track d_{xy} significance;good track d_{xy} significance", 200, 0, 20, 400, 0, 40);
	h_2D_poor_dz_sig_good_dz_sig_large_sv_fig2 = fs->make<TH2F>("h_2D_poor_dz_sig_good_dz_sig_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;poor track d_{z} significance;good track d_{z} significance", 100, 0, 100, 100, 0, 100);
	h_2D_poor_pt_good_pt_large_sv_fig2 = fs->make<TH2F>("h_2D_poor_pt_good_pt_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;poor track p_{T} (GeV);good track p_{T} (GeV)", 100, 0, 100, 100, 0, 100);
	h_2D_poor_avg_pt_good_avg_pt_large_sv_fig2 = fs->make<TH2F>("h_2D_poor_avg_pt_good_avg_pt_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;poor track avg p_{T} (GeV);good track avg p_{T} (GeV)", 100, 0, 100, 100, 0, 100);
	h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_fig2 = fs->make<TH2F>("h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;poor track avg p_{T} sig;good track avg p_{T} sig", 21, 0, 210, 21, 0, 210);
	h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_fig2 = fs->make<TH2F>("h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;minor sv's sum p_{T} ratio(shared/minor tracks);major sv's sum p_{T} ratio(shared/major tracks)", 200, 0, 2, 200, 0, 2);
	h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_fig2 = fs->make<TH2F>("h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;minor sv's avg p_{T} ratio(shared/minor tracks);major sv's avg p_{T} ratio(shared/major tracks)", 200, 0, 2, 200, 0, 2);
	h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_fig2 = fs->make<TH2F>("h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;major sv's non-shared jet avg p_{T} (GeV);shared jet p_{T} (GeV)", 30, 0, 3000, 20, 0, 3000);
	h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_fig2 = fs->make<TH2F>("h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;minor sv's non-shared jet avg p_{T} (GeV);shared jet p_{T} (GeV)", 30, 0, 3000, 20, 0, 3000);
	h_pt_poor_shared_tracks_large_nsv2_fig2 = fs->make<TH1F>("h_pt_poor_shared_tracks_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; poor track p_{T} (GeV);arb. units", 100, 0, 100);
	h_pt_good_shared_tracks_large_nsv2_fig2 = fs->make<TH1F>("h_pt_good_shared_tracks_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; good track p_{T} (GeV);arb. units", 100, 0, 100);
	h_dxy_err_poor_shared_tracks_large_nsv2_fig2 = fs->make<TH1F>("h_dxy_err_poor_shared_tracks_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; poor track d_{xy} (cm);arb. units", 50, 0, 0.05);
	h_dxy_err_good_shared_tracks_large_nsv2_fig2 = fs->make<TH1F>("h_dxy_err_good_shared_tracks_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; good track d_{xy} (cm);arb. units", 50, 0, 0.05);
	h_dxy_sig_poor_shared_tracks_large_nsv2_fig2 = fs->make<TH1F>("h_dxy_sig_poor_shared_tracks_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; poor track d_{xy} significance;arb. units", 400, 0, 40);
	h_dxy_sig_good_shared_tracks_large_nsv2_fig2 = fs->make<TH1F>("h_dxy_sig_good_shared_tracks_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; good track d_{xy} significance;arb. units", 400, 0, 40);
	h_dPhi_sv_large_nsv2_fig2 = fs->make<TH1F>("h_dPhi_sv_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;delta(phi of SV0, phi of SV1);arb. units", 31, 0, 3.16);
	h_dEta_sv_large_nsv2_fig2 = fs->make<TH1F>("h_dEta_sv_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;delta(eta of SV0, eta of SV1);arb. units", 600, 0, 6);
	h_2D_dEta_dPhi_sv_large_nsv2_fig2 = fs->make<TH2F>("h_2D_dEta_dPhi_sv_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2;delta(eta of SV0, eta of SV1);delta(phi of SV0, phi of SV1)", 600, 0, 6, 31, 0, 3.16);
	h_pt_shared_jets_large_nsv2_fig2 = fs->make<TH1F>("h_pt_shared_jets_large_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, fig2; shared jets p_{T} (GeV);arb. units", 30, 0, 3000);
	h_more_ratio_sum_major_sv_nsv2_fig2 = fs->make<TH1F>("h_more_ratio_sum_major_sv_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, Fig2; more sum p_{T} ratio(shared/SV's tracks) is major SV ?", 2, 0, 2);
	h_more_ratio_avg_major_sv_nsv2_fig2 = fs->make<TH1F>("h_more_ratio_avg_major_sv_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, Fig2; more avg p_{T} ratio(shared/SV's tracks) is major SV ?", 2, 0, 2);
	h_more_pT_avg_major_sv_nsv2_fig2 = fs->make<TH1F>("h_more_pT_avg_major_sv_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, Fig2; more avg p_{T} shared tracks is major SV ?", 2, 0, 2);
	h_more_pT_sig_avg_major_sv_nsv2_fig2 = fs->make<TH1F>("h_more_pT_sig_avg_major_sv_nsv2_fig2", "nsv = 2, absdeltaphi01 > 0.5, Fig2; more avg p_{T} sig shared tracks is major SV ?", 2, 0, 2);


	h_dPhi_good_tracks_its_sv_large_semi_fig2 = fs->make<TH1F>("h_dPhi_good_tracks_its_sv_large_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;delta(phi of good tracks, phi of major vertex);arb. units", 31, 0, 3.16);
	h_dPhi_poor_tracks_its_sv_large_semi_fig2 = fs->make<TH1F>("h_dPhi_poor_tracks_its_sv_large_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;delta(phi of poor tracks, phi of minor vertex);arb. units", 31, 0, 3.16);
	h_dPhi_good_tracks_other_sv_large_semi_fig2 = fs->make<TH1F>("h_dPhi_good_tracks_other_sv_large_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;delta(phi of good tracks, phi of minor vertex);arb. units", 31, 0, 3.16);
	h_dPhi_poor_tracks_other_sv_large_semi_fig2 = fs->make<TH1F>("h_dPhi_poor_tracks_other_sv_large_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;delta(phi of poor tracks, phi of major vertex);arb. units", 31, 0, 3.16);
	h_vertex_chi2dof_absdeltaphi0_large_nsv2_semi_fig2 = fs->make<TH1F>("h_vertex_chi2dof_absdeltaphi0_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; sv opposite to the shared jet's tracks #chi^ { 2 } / dof ;arb. units", 100, 0, 10);
	h_absdeltaphi_large_jet_shared_tracks_nsv2_semi_fig2 = fs->make<TH1F>("h_absdeltaphi_large_jet_shared_tracks_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; abs(delta(shared tracks, phi of the only-one shared jet));arb. units", 31, 0, 3.16);
	h_ratio_ntracks_large_nsv2_shared_jets_semi_fig2 = fs->make<TH1F>("h_ratio_ntracks_large_nsv2_shared_jets_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;ratios of shared tracks (>=1);arb. units", 50, 0, 10);
	h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;miss dist sig(1 poor trk, minor vertex);avg miss dist sig(>=4 good trk, major vertex)", 200, 0, 20, 200, 0, 20);
	h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;miss dist sig(1 poor trk, major vertex);avg miss dist sig(>=4 good trk, minor vertex)", 100, 0, 100, 100, 0, 100);
	h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;miss dist sig(1 poor trk, minor vertex);miss dist sig(1 poor trk, major vertex)", 200, 0, 20, 100, 0, 100);
	h_2D_pt_dxy_err_absdeltaphi0_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_pt_dxy_err_absdeltaphi0_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;max track p_{T} (GeV);max track d_{xy} err (cm)", 100, 0, 100, 50, 0, 0.05);
	h_2D_pt_dxy_err_absdeltaphi1_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_pt_dxy_err_absdeltaphi1_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;min track p_{T} (GeV);min track d_{xy} err (cm)", 100, 0, 100, 50, 0, 0.05);
	h_2D_poor_dxy_err_good_dxy_err_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_poor_dxy_err_good_dxy_err_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;poor track d_{xy} err (cm);good track d_{xy} err (cm)", 50, 0, 0.05, 50, 0, 0.05);
	h_2D_poor_dz_err_good_dz_err_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_poor_dz_err_good_dz_err_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;poor track d_{z} err (cm);good track d_{z} err (cm)", 50, 0, 0.05, 50, 0, 0.05);
	h_2D_poor_dxy_sig_good_dxy_sig_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_poor_dxy_sig_good_dxy_sig_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;poor track d_{xy} significance;good track d_{xy} significance", 200, 0, 20, 400, 0, 40);
	h_2D_poor_dz_sig_good_dz_sig_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_poor_dz_sig_good_dz_sig_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;poor track d_{z} significance;good track d_{z} significance", 100, 0, 100, 100, 0, 100);
	h_2D_poor_pt_good_pt_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_poor_pt_good_pt_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;poor track p_{T} (GeV);good track p_{T} (GeV)", 100, 0, 100, 100, 0, 100);
	h_2D_poor_avg_pt_good_avg_pt_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_poor_avg_pt_good_avg_pt_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;poor track avg p_{T} (GeV);good track avg p_{T} (GeV)", 100, 0, 100, 100, 0, 100);
	h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;poor track avg p_{T} sig;good track avg p_{T} sig", 21, 0, 210, 21, 0, 210);
	h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;minor sv's sum p_{T} ratio(shared/minor tracks);major sv's sum p_{T} ratio(shared/major tracks)", 200, 0, 2, 200, 0, 2);
	h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;minor sv's avg p_{T} ratio(shared/minor tracks);major sv's avg p_{T} ratio(shared/major tracks)", 200, 0, 2, 200, 0, 2);
	h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;major sv's non-shared jet avg p_{T} (GeV);shared jet p_{T} (GeV)", 30, 0, 3000, 20, 0, 3000);
	h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2 = fs->make<TH2F>("h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;minor sv's non-shared jet avg p_{T} (GeV);shared jet p_{T} (GeV)", 30, 0, 3000, 20, 0, 3000);
	h_pt_poor_shared_tracks_large_nsv2_semi_fig2 = fs->make<TH1F>("h_pt_poor_shared_tracks_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; poor track p_{T} (GeV);arb. units", 100, 0, 100);
	h_pt_good_shared_tracks_large_nsv2_semi_fig2 = fs->make<TH1F>("h_pt_good_shared_tracks_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; good track p_{T} (GeV);arb. units", 100, 0, 100);
	h_dxy_err_poor_shared_tracks_large_nsv2_semi_fig2 = fs->make<TH1F>("h_dxy_err_poor_shared_tracks_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; poor track d_{xy} (cm);arb. units", 50, 0, 0.05);
	h_dxy_err_good_shared_tracks_large_nsv2_semi_fig2 = fs->make<TH1F>("h_dxy_err_good_shared_tracks_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; good track d_{xy} (cm);arb. units", 50, 0, 0.05);
	h_dxy_sig_poor_shared_tracks_large_nsv2_semi_fig2 = fs->make<TH1F>("h_dxy_sig_poor_shared_tracks_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; poor track d_{xy} significance;arb. units", 400, 0, 40);
	h_dxy_sig_good_shared_tracks_large_nsv2_semi_fig2 = fs->make<TH1F>("h_dxy_sig_good_shared_tracks_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; good track d_{xy} significance;arb. units", 400, 0, 40);
	h_dPhi_sv_large_nsv2_semi_fig2 = fs->make<TH1F>("h_dPhi_sv_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;delta(phi of SV0, phi of SV1);arb. units", 31, 0, 3.16);
	h_dEta_sv_large_nsv2_semi_fig2 = fs->make<TH1F>("h_dEta_sv_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;delta(eta of SV0, eta of SV1);arb. units", 600, 0, 6);
	h_2D_dEta_dPhi_sv_large_nsv2_semi_fig2 = fs->make<TH2F>("h_2D_dEta_dPhi_sv_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2;delta(eta of SV0, eta of SV1);delta(phi of SV0, phi of SV1)", 600, 0, 6, 31, 0, 3.16);
	h_pt_shared_jets_large_nsv2_semi_fig2 = fs->make<TH1F>("h_pt_shared_jets_large_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-fig2; shared jets p_{T} (GeV);arb. units", 30, 0, 3000);
	h_more_ratio_sum_major_sv_nsv2_semi_fig2 = fs->make<TH1F>("h_more_ratio_sum_major_sv_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-Fig2; more sum p_{T} ratio(shared/SV's tracks) is major SV ?", 2, 0, 2);
	h_more_ratio_avg_major_sv_nsv2_semi_fig2 = fs->make<TH1F>("h_more_ratio_avg_major_sv_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-Fig2; more avg p_{T} ratio(shared/SV's tracks) is major SV ?", 2, 0, 2);
	h_more_pT_avg_major_sv_nsv2_semi_fig2 = fs->make<TH1F>("h_more_pT_avg_major_sv_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-Fig2; more avg p_{T} shared tracks is major SV ?", 2, 0, 2);
	h_more_pT_sig_avg_major_sv_nsv2_semi_fig2 = fs->make<TH1F>("h_more_pT_sig_avg_major_sv_nsv2_semi_fig2", "nsv = 2, absdeltaphi01 > 0.5, semi-Fig2; more avg p_{T} sig shared tracks is major SV ?", 2, 0, 2);

	
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

	if (nsv == 1) {
		int njets = std::set<double>(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end()).size();
		h_sv_njets_nsv1->Fill(njets, w);
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

		if ((nsv==2) && (fabs(reco::deltaPhi(phi0, phi1)) > 0.5)) {
			std::vector<int> sv0_track_which_idx(int(sv0.ntracks()));
			int idx0 = 0;
			std::generate(sv0_track_which_idx.begin(), sv0_track_which_idx.end(), [&] { return idx0++; });
			std::vector<int> sv1_track_which_idx(int(sv1.ntracks()));
			int idx1 = 0;
			std::generate(sv1_track_which_idx.begin(), sv1_track_which_idx.end(), [&] { return idx1++; });
			/*
			if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx.size() >= 5)) {
				if (sv0_track_which_idx.size() >= sv1_track_which_idx.size()) {
					h_2D_sv_tracks_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size());
				}
				else {
					h_2D_sv_tracks_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size());

				}
			}
			*/
			if (sv0_track_which_idx.size() >= sv1_track_which_idx.size()) {
				h_2D_sv_tracks_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
			}
			else {
				h_2D_sv_tracks_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);

			}

			
			int njets_sv0 = std::set<double>(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end()).size();
			int njets_sv1 = std::set<double>(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()).size();
			h_sv_njets_large_nsv2_no_shj->Fill(njets_sv0, w);
			h_sv_njets_large_nsv2_no_shj->Fill(njets_sv1, w);
			
			
		}

		bool shared_jet = std::find_first_of(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()) != sv_track_which_jet[0].end();
		if (shared_jet) {

			std::vector<int> sv0_track_which_idx(int(sv0.ntracks()));
			int idx0 = 0;
			std::generate(sv0_track_which_idx.begin(), sv0_track_which_idx.end(), [&] { return idx0++; });
			std::vector<int> sv1_track_which_idx(int(sv1.ntracks()));
			int idx1 = 0;
			std::generate(sv1_track_which_idx.begin(), sv1_track_which_idx.end(), [&] { return idx1++; });

			if ((nsv == 2) && (fabs(reco::deltaPhi(phi0, phi1)) > 0.5)) {
				
				if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx.size() >= 5)) {
					if (sv0_track_which_idx.size() >= sv1_track_which_idx.size()) {
						h_2D_sv_tracks_shared_jets_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
					}
					else {
						h_2D_sv_tracks_shared_jets_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);

					}
				}

				int njets_sv0 = std::set<double>(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end()).size();
				int njets_sv1 = std::set<double>(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()).size();
				h_sv_njets_large_nsv2_shj->Fill(njets_sv0, w);
				h_sv_njets_large_nsv2_shj->Fill(njets_sv1, w);

				

			}

			int nsharedjets = 1;
			std::vector<int> nsharedjet_jet_index;
			std::vector<std::vector<int> > sv_track_which_jet_copy;
			sv_track_which_jet_copy = sv_track_which_jet;
			std::vector<int> nsharedjet_tracks_sv0;
			std::vector<int> nsharedjet_tracks_sv1;
			std::vector<std::vector<int> >sv0_sharedjet_which_idx;
			std::vector<std::vector<int> >sv1_sharedjet_which_idx;
			std::vector<std::vector<int> >sv0_sharedjet_which_no_trk_idx;
			std::vector<std::vector<int> >sv1_sharedjet_which_no_trk_idx;
			std::vector<int>::iterator it = std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end());
			int idx = std::distance(sv_track_which_jet_copy[0].begin(), it);
			int jet_index = sv_track_which_jet_copy[0].at(idx);
			nsharedjet_jet_index.push_back(jet_index);
			sv_track_which_jet_copy[0].erase(std::remove(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), jet_index), sv_track_which_jet_copy[0].end());
			sv_track_which_jet_copy[1].erase(std::remove(sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end(), jet_index), sv_track_which_jet_copy[1].end());
			nsharedjet_tracks_sv0.push_back(std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index));

			std::vector<int> sv0_track_which_jet = sv_track_which_jet[0];
			std::vector<int> sv0_track_which_jet_idx = sv_track_which_idx[0];

			std::vector<int> sv0_track_which_idx_copy = sv0_track_which_idx;
			std::vector<int> sv0_track_which_idx_no_trk = sv0_track_which_idx;
			std::vector<int> sv0_track_which_temp_idx;
			std::multimap<int, size_t> sv0_m_nshj1;



			for (size_t k = 0; k < sv0_track_which_jet.size(); k++) if (sv0_track_which_jet[k] == jet_index) { sv0_m_nshj1.insert({ sv0_track_which_jet[k], k }); }

			for (auto it = sv0_m_nshj1.begin(); it != sv0_m_nshj1.end(); )
			{
				auto p = sv0_m_nshj1.equal_range(it->first);

				while (p.first != p.second)
				{
					sv0_track_which_temp_idx.push_back(sv0_track_which_jet_idx[p.first++->second]);
				}
				it = p.second;

			}

			sv0_sharedjet_which_idx.push_back(sv0_track_which_temp_idx);
			for (size_t k = 0; k < sv0_track_which_temp_idx.size(); k++) {
				int track_index = sv0_track_which_temp_idx[k];
				sv0_track_which_idx_copy.erase(std::remove(sv0_track_which_idx_copy.begin(), sv0_track_which_idx_copy.end(), track_index), sv0_track_which_idx_copy.end());
				sv0_track_which_idx_no_trk.erase(std::remove(sv0_track_which_idx_no_trk.begin(), sv0_track_which_idx_no_trk.end(), track_index), sv0_track_which_idx_no_trk.end());

			}
			sv0_sharedjet_which_no_trk_idx.push_back(sv0_track_which_idx_copy);
            sv0_track_which_temp_idx = {};

			nsharedjet_tracks_sv1.push_back(std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index));

			std::vector<int> sv1_track_which_jet = sv_track_which_jet[1];
			std::vector<int> sv1_track_which_jet_idx = sv_track_which_idx[1];

			std::vector<int> sv1_track_which_idx_copy = sv1_track_which_idx;
			std::vector<int> sv1_track_which_idx_no_trk = sv1_track_which_idx;
			std::vector<int> sv1_track_which_temp_idx;
			std::multimap<int, size_t> sv1_m_nshj1;

			for (size_t k = 0; k < sv1_track_which_jet.size(); k++) if (sv1_track_which_jet[k] == jet_index) { sv1_m_nshj1.insert({ sv1_track_which_jet[k], k }); }

			for (auto it = sv1_m_nshj1.begin(); it != sv1_m_nshj1.end(); )
			{
				auto p = sv1_m_nshj1.equal_range(it->first);

				while (p.first != p.second)
				{
					sv1_track_which_temp_idx.push_back(sv1_track_which_jet_idx[p.first++->second]);
				}
				it = p.second;

			}
			sv1_sharedjet_which_idx.push_back(sv1_track_which_temp_idx);  
			for (size_t k = 0; k < sv1_track_which_temp_idx.size(); k++) {
				int track_index = sv1_track_which_temp_idx[k];
				sv1_track_which_idx_copy.erase(std::remove(sv1_track_which_idx_copy.begin(), sv1_track_which_idx_copy.end(), track_index), sv1_track_which_idx_copy.end());
				sv1_track_which_idx_no_trk.erase(std::remove(sv1_track_which_idx_no_trk.begin(), sv1_track_which_idx_no_trk.end(), track_index), sv1_track_which_idx_no_trk.end());

			}
			sv1_sharedjet_which_no_trk_idx.push_back(sv1_track_which_idx_copy);
			sv1_track_which_temp_idx = {};
			
			while (std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end()) != sv_track_which_jet_copy[0].end()) {
				nsharedjets++;
				sv0_track_which_idx_copy = sv0_track_which_idx;
				sv1_track_which_idx_copy = sv1_track_which_idx;
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
						sv0_track_which_temp_idx.push_back(sv0_track_which_jet_idx[p.first++->second]);
					}
					it = p.second;

				}
				
				sv0_sharedjet_which_idx.push_back(sv0_track_which_temp_idx); 
				for (size_t k = 0; k < sv0_track_which_temp_idx.size(); k++) {
					int track_index = sv0_track_which_temp_idx[k];
					sv0_track_which_idx_copy.erase(std::remove(sv0_track_which_idx_copy.begin(), sv0_track_which_idx_copy.end(), track_index), sv0_track_which_idx_copy.end());
					sv0_track_which_idx_no_trk.erase(std::remove(sv0_track_which_idx_no_trk.begin(), sv0_track_which_idx_no_trk.end(), track_index), sv0_track_which_idx_no_trk.end());

				}
				sv0_sharedjet_which_no_trk_idx.push_back(sv0_track_which_idx_copy);
			    
				if (sv0_track_which_temp_idx.size() + sv0_track_which_idx_copy.size() != sv0_track_which_idx.size()) {
					std::cout << "sv0 needs to be fixed" << std::endl;
					std::cout << "sv0 tracks = " << sv0_track_which_idx.size() << ", shared ones = " << sv0_track_which_temp_idx.size() << ", not shared ones = " << sv0_track_which_idx_copy.size() << std::endl;

				}
				sv0_track_which_temp_idx = {};

				nsharedjet_tracks_sv1.push_back(std::count(sv1_track_which_jet.begin(), sv1_track_which_jet.end(), jet_index));
				std::multimap<int, size_t> sv1_m;
				for (size_t k = 0; k < sv1_track_which_jet.size(); k++) if (sv1_track_which_jet[k] == jet_index) { sv1_m.insert({ sv1_track_which_jet[k], k }); }

				for (auto it = sv1_m.begin(); it != sv1_m.end(); )
				{
					auto p = sv1_m.equal_range(it->first);

					while (p.first != p.second)
					{
						sv1_track_which_temp_idx.push_back(sv1_track_which_jet_idx[p.first++->second]);
					}
					it = p.second;

				}
				
				sv1_sharedjet_which_idx.push_back(sv1_track_which_temp_idx);  
				for (size_t k = 0; k < sv1_track_which_temp_idx.size(); k++) {
					int track_index = sv1_track_which_temp_idx[k];
					sv1_track_which_idx_copy.erase(std::remove(sv1_track_which_idx_copy.begin(), sv1_track_which_idx_copy.end(), track_index), sv1_track_which_idx_copy.end());
					sv1_track_which_idx_no_trk.erase(std::remove(sv1_track_which_idx_no_trk.begin(), sv1_track_which_idx_no_trk.end(), track_index), sv1_track_which_idx_no_trk.end());

				}
				sv1_sharedjet_which_no_trk_idx.push_back(sv1_track_which_idx_copy);

				if (sv1_track_which_temp_idx.size() + sv1_track_which_idx_copy.size() != sv1_track_which_idx.size()) {
					std::cout << "sv1 needs to be fixed" << std::endl;
					std::cout << "sv1 tracks = " << sv1_track_which_idx.size() << ", shared ones = " << sv1_track_which_temp_idx.size() << ", not shared ones = " << sv1_track_which_idx_copy.size() << std::endl;

				}

				sv1_track_which_temp_idx = {};
				
			
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

					std::vector<int> sv0_track_which_idx_no_shared_track = sv0_track_which_idx_no_trk;
					std::vector<int> sv1_track_which_idx_no_shared_track = sv1_track_which_idx_no_trk;
					std::vector<int> sv0_track_which_idx(int(sv0.ntracks()));
					int idx0 = 0;
					std::generate(sv0_track_which_idx.begin(), sv0_track_which_idx.end(), [&] { return idx0++; });
					std::vector<int> sv1_track_which_idx(int(sv1.ntracks()));
					int idx1 = 0;
					std::generate(sv1_track_which_idx.begin(), sv1_track_which_idx.end(), [&] { return idx1++; });

					if ((sv0_track_which_idx.size() - sv0_track_which_idx_no_shared_track.size()) >= (sv1_track_which_idx.size() - sv1_track_which_idx_no_shared_track.size())) {		 // sv0 is major vtx


						if ((sv0_track_which_idx_no_shared_track.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
							h_2D_sv_tracks_no_shared_tracks_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx_no_shared_track.size(),w);

						}

					}
					else {

						if ((sv1_track_which_idx_no_shared_track.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
							h_2D_sv_tracks_no_shared_tracks_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx_no_shared_track.size(),w);
						}

					}
					
					if ((sv0_track_which_idx.size()- sv0_track_which_idx_no_shared_track.size()) >= (sv1_track_which_idx.size() - sv1_track_which_idx_no_shared_track.size())) {		 // sv0 is major vtx
						
						
						if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
								h_2D_sv_tracks_no_minor_shared_tracks_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
							
						}

					}
					else {
						
						if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {	  
								h_2D_sv_tracks_no_minor_shared_tracks_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
						}

					}


					
					

					

					h_nsharedjets_large_nsv2_shared_jets->Fill(nsharedjets, w);
					h_svdist2d_large_absdeltaphi01_nsv2_shared_jets->Fill(svdist2d, w);
					h_svdist3d_large_absdeltaphi01_nsv2_shared_jets->Fill(svdist3d, w);

					std::vector<int> sv0_track_which_non_shared_jet = sv_track_which_jet_copy[0];
					std::vector<int> sv1_track_which_non_shared_jet = sv_track_which_jet_copy[1];

					double sum_sv0_non_shared_jet_pt = 0;
					double sum_sv1_non_shared_jet_pt = 0;
					for (size_t k = 0; k < sv0_track_which_non_shared_jet.size(); k++) {
						sum_sv0_non_shared_jet_pt = sum_sv0_non_shared_jet_pt + mevent->jet_pt[int(sv0_track_which_non_shared_jet[k])];
					}
					for (size_t k = 0; k < sv1_track_which_non_shared_jet.size(); k++) {
						sum_sv1_non_shared_jet_pt = sum_sv1_non_shared_jet_pt + mevent->jet_pt[int(sv1_track_which_non_shared_jet[k])];
					}
					double avg_sv0_non_shared_jet_pt = sum_sv0_non_shared_jet_pt/ sv0_track_which_non_shared_jet.size();
					double avg_sv1_non_shared_jet_pt = sum_sv1_non_shared_jet_pt / sv1_track_which_non_shared_jet.size();

					std::vector<int> sv0_sum_pt_track_which_idx = sv0_track_which_idx;
					std::vector<int> sv1_sum_pt_track_which_idx = sv1_track_which_idx;
					std::vector<int> sv0_ratio_sum_pt_track_which_idx = sv0_track_which_idx;
					std::vector<int> sv1_ratio_sum_pt_track_which_idx = sv1_track_which_idx;
					std::vector<int> sv0_avg_pt_track_which_idx = sv0_track_which_idx;
					std::vector<int> sv1_avg_pt_track_which_idx = sv1_track_which_idx;
					std::vector<int> sv0_ratio_avg_pt_track_which_idx = sv0_track_which_idx;
					std::vector<int> sv1_ratio_avg_pt_track_which_idx = sv1_track_which_idx;
					int njets_sum_pT_sv0 = std::set<double>(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end()).size();
					int njets_sum_pT_sv1 = std::set<double>(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()).size();
					int njets_ratio_sum_pT_sv0 = std::set<double>(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end()).size();
					int njets_ratio_sum_pT_sv1 = std::set<double>(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()).size();

                                        for (int i = 0; i < nsharedjets; i++) {				//start nsharedjet loop


                                                
						int jet_index = nsharedjet_jet_index[i];

						double dphi_large_sv0_sharedjet = double(fabs(reco::deltaPhi(phi0, mevent->jet_phi[jet_index])));
						double dphi_large_sv1_sharedjet = double(fabs(reco::deltaPhi(phi1, mevent->jet_phi[jet_index])));


						double ratio_ntracks_nsv2;
						if (nsharedjet_tracks_sv0[i] > nsharedjet_tracks_sv1[i]) { ratio_ntracks_nsv2 = nsharedjet_tracks_sv0[i] / nsharedjet_tracks_sv1[i]; }
						else { ratio_ntracks_nsv2 = nsharedjet_tracks_sv1[i] / nsharedjet_tracks_sv0[i]; }
						h_ratio_ntracks_large_nsv2_shared_jets->Fill(ratio_ntracks_nsv2);

						// all w/ sum pt 

						std::vector<int> sv1_diff; 
						std::vector<int> sv0_diff;

						double sum_pt_i_sv0 = 0;
						std::vector<int> sv0_i_sharedjet_which_idx = sv0_sharedjet_which_idx[i];
						for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {
							int idx = sv0_i_sharedjet_which_idx[j];
							sum_pt_i_sv0 = sum_pt_i_sv0 + sv0.track_pt(idx);
						}
						double sum_pt_i_sv1 = 0;
						std::vector<int> sv1_i_sharedjet_which_idx = sv1_sharedjet_which_idx[i];
						for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {
							int idx = sv1_i_sharedjet_which_idx[j];
							sum_pt_i_sv1 = sum_pt_i_sv1 + sv1.track_pt(idx);
						}

						double sum_pt_i_no_sv0 = 0;
						std::vector<int> sv0_i_sharedjet_which_no_idx = sv0_sharedjet_which_no_trk_idx[i];
						for (unsigned int j = 0; j < sv0_i_sharedjet_which_no_idx.size(); j++) {
							int idx = sv0_i_sharedjet_which_no_idx[j];
							sum_pt_i_no_sv0 = sum_pt_i_no_sv0 + sv0.track_pt(idx);
						}
						double sum_pt_i_no_sv1 = 0;
						std::vector<int> sv1_i_sharedjet_which_no_idx = sv1_sharedjet_which_no_trk_idx[i];
						for (unsigned int j = 0; j < sv1_i_sharedjet_which_no_idx.size(); j++) {
							int idx = sv1_i_sharedjet_which_no_idx[j];
							sum_pt_i_no_sv1 = sum_pt_i_no_sv1 + sv1.track_pt(idx);
						}

						
						h_diff_pT_sum_sv0_sv1_sv_nsv2->Fill(sum_pt_i_sv0- sum_pt_i_sv1,w);
						h_ratio_diff_pT_sum_sv0_sv1_sv_nsv2->Fill((sum_pt_i_sv0 - sum_pt_i_sv1)/(sum_pt_i_sv0 + sum_pt_i_sv1), w);
						

						if (sum_pt_i_sv0 >= sum_pt_i_sv1) {
							h_diff_pT_sum_dPhi_shj_sv_large_nsv2->Fill(dphi_large_sv0_sharedjet, w);
							std::set_difference(sv1_sum_pt_track_which_idx.begin(), sv1_sum_pt_track_which_idx.end(), sv1_i_sharedjet_which_idx.begin(), sv1_i_sharedjet_which_idx.end(),
								std::inserter(sv1_diff, sv1_diff.begin())); 
							
							sv1_sum_pt_track_which_idx = sv1_diff;
							njets_sum_pT_sv1 = njets_sum_pT_sv1 - 1;
						}
						else {
							h_diff_pT_sum_dPhi_shj_sv_large_nsv2->Fill(dphi_large_sv1_sharedjet, w);
							std::set_difference(sv0_sum_pt_track_which_idx.begin(), sv0_sum_pt_track_which_idx.end(), sv0_i_sharedjet_which_idx.begin(), sv0_i_sharedjet_which_idx.end(),
								std::inserter(sv0_diff, sv0_diff.begin()));
							
							sv0_sum_pt_track_which_idx = sv0_diff;
							njets_sum_pT_sv0 = njets_sum_pT_sv0 - 1;

						}

						//1D all w/ pT variables 
						std::vector<int> sv1_diff2;
						std::vector<int> sv0_diff2;

						h_diff_pT_avg_sv0_sv1_sv_nsv2->Fill((sum_pt_i_sv0 / nsharedjet_tracks_sv0[i]) - (sum_pt_i_sv1 / nsharedjet_tracks_sv1[i]), w);
						h_ratio_diff_pT_avg_sv0_sv1_sv_nsv2->Fill(((sum_pt_i_sv0 / nsharedjet_tracks_sv0[i]) - (sum_pt_i_sv1 / nsharedjet_tracks_sv1[i]))/((sum_pt_i_sv0 / nsharedjet_tracks_sv0[i]) + (sum_pt_i_sv1 / nsharedjet_tracks_sv1[i])), w);

						if ((sum_pt_i_sv0/ nsharedjet_tracks_sv0[i]) >= (sum_pt_i_sv1/ nsharedjet_tracks_sv1[i])) {
							h_diff_pT_avg_dPhi_shj_sv_large_nsv2->Fill(dphi_large_sv0_sharedjet, w);
							std::set_difference(sv1_avg_pt_track_which_idx.begin(), sv1_avg_pt_track_which_idx.end(), sv1_i_sharedjet_which_idx.begin(), sv1_i_sharedjet_which_idx.end(),
								std::inserter(sv1_diff2, sv1_diff2.begin()));

							sv1_avg_pt_track_which_idx = sv1_diff2;
						}
						else {
							h_diff_pT_avg_dPhi_shj_sv_large_nsv2->Fill(dphi_large_sv1_sharedjet, w);
							std::set_difference(sv0_avg_pt_track_which_idx.begin(), sv0_avg_pt_track_which_idx.end(), sv0_i_sharedjet_which_idx.begin(), sv0_i_sharedjet_which_idx.end(),
								std::inserter(sv0_diff2, sv0_diff2.begin()));

							sv0_avg_pt_track_which_idx = sv0_diff2;

						}

						std::vector<int> sv1_diff3;
						std::vector<int> sv0_diff3;

						h_diff_ratio_pT_sum_sv0_sv1_sv_nsv2->Fill((sum_pt_i_sv0 / (sum_pt_i_sv0 + sum_pt_i_no_sv0)) - (sum_pt_i_sv1 / (sum_pt_i_sv1 + sum_pt_i_no_sv1)), w);

						if ((sum_pt_i_sv0 / (sum_pt_i_sv0 + sum_pt_i_no_sv0)) >= (sum_pt_i_sv1 / (sum_pt_i_sv1 + sum_pt_i_no_sv1))) {
							h_diff_ratio_pT_sum_dPhi_shj_sv_large_nsv2->Fill(dphi_large_sv0_sharedjet, w);
							std::set_difference(sv1_ratio_sum_pt_track_which_idx.begin(), sv1_ratio_sum_pt_track_which_idx.end(), sv1_i_sharedjet_which_idx.begin(), sv1_i_sharedjet_which_idx.end(),
								std::inserter(sv1_diff3, sv1_diff3.begin()));

							sv1_ratio_sum_pt_track_which_idx = sv1_diff3;
							njets_ratio_sum_pT_sv1 = njets_ratio_sum_pT_sv1 - 1;
						}
						else {
							h_diff_ratio_pT_sum_dPhi_shj_sv_large_nsv2->Fill(dphi_large_sv1_sharedjet, w);
							std::set_difference(sv0_ratio_sum_pt_track_which_idx.begin(), sv0_ratio_sum_pt_track_which_idx.end(), sv0_i_sharedjet_which_idx.begin(), sv0_i_sharedjet_which_idx.end(),
								std::inserter(sv0_diff3, sv0_diff3.begin()));

							sv0_ratio_sum_pt_track_which_idx = sv0_diff3;
							njets_ratio_sum_pT_sv0 = njets_ratio_sum_pT_sv0 - 1;

						}

						std::vector<int> sv1_diff4;
						std::vector<int> sv0_diff4;

						h_diff_ratio_pT_avg_sv0_sv1_sv_nsv2->Fill(((sum_pt_i_sv0* sv0_track_which_idx.size()) / ((sum_pt_i_sv0 + sum_pt_i_no_sv0) * nsharedjet_tracks_sv0[i])) - ((sum_pt_i_sv1 * sv1_track_which_idx.size()) / ((sum_pt_i_sv1 + sum_pt_i_no_sv1) * nsharedjet_tracks_sv1[i])), w);

						if (((sum_pt_i_sv0* sv0_track_which_idx.size()) /((sum_pt_i_sv0 + sum_pt_i_no_sv0)* nsharedjet_tracks_sv0[i])) >= ((sum_pt_i_sv1* sv1_track_which_idx.size()) /((sum_pt_i_sv1 + sum_pt_i_no_sv1)* nsharedjet_tracks_sv1[i]))) {
							h_diff_ratio_pT_avg_dPhi_shj_sv_large_nsv2->Fill(dphi_large_sv0_sharedjet, w);
							std::set_difference(sv1_ratio_avg_pt_track_which_idx.begin(), sv1_ratio_avg_pt_track_which_idx.end(), sv1_i_sharedjet_which_idx.begin(), sv1_i_sharedjet_which_idx.end(),
								std::inserter(sv1_diff4, sv1_diff4.begin()));

							sv1_ratio_avg_pt_track_which_idx = sv1_diff4;
						}
						else {
							h_diff_ratio_pT_avg_dPhi_shj_sv_large_nsv2->Fill(dphi_large_sv1_sharedjet, w);
							std::set_difference(sv0_ratio_avg_pt_track_which_idx.begin(), sv0_ratio_avg_pt_track_which_idx.end(), sv0_i_sharedjet_which_idx.begin(), sv0_i_sharedjet_which_idx.end(),
								std::inserter(sv0_diff4, sv0_diff4.begin()));

							sv0_ratio_avg_pt_track_which_idx = sv0_diff4;

						}




						

						if (nsharedjets == 1) {

							if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx.size() >= 5)) {
								if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
									h_2D_sv_tracks_shared_jets_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
								}
								else {
									h_2D_sv_tracks_shared_jets_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);

								}
							}
							else {
								if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
									h_2D_poor_sv_tracks_shared_jets_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
								}
								else {
									h_2D_poor_sv_tracks_shared_jets_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);

								}

							}

							if ((sv0_track_which_idx_no_shared_track.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
								if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
									h_2D_sv_tracks_no_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx_no_shared_track.size(),w);
								}
								else {
									h_2D_sv_tracks_no_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx_no_shared_track.size(),w);

								}
							}
							else {
								if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
									h_2D_poor_sv_tracks_no_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx_no_shared_track.size(),w);
								}
								else {
									h_2D_poor_sv_tracks_no_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx_no_shared_track.size(),w);

								}

							}




							if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
								if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
									h_2D_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
								}
								else {
									h_2D_poor_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

								}
							}
							else {
								if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
									h_2D_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
								}
								else {
									h_2D_poor_sv_tracks_no_minor_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);


								}

							}




							std::vector<int> sv1_sharedjet_1_which_trk_idx = sv1_sharedjet_which_idx[0];
							std::vector<int> sv1_sharedjet_1_which_no_trk_idx = sv1_sharedjet_which_no_trk_idx[0];
							std::vector<int> sv0_sharedjet_1_which_trk_idx = sv0_sharedjet_which_idx[0];
							std::vector<int> sv0_sharedjet_1_which_no_trk_idx = sv0_sharedjet_which_no_trk_idx[0];

							double sum_pt_shared_sv1 = 0.0;
							double sum_pt_no_shared_sv1 = 0.0;
							double sum_pt_shared_sv0 = 0.0;
							double sum_pt_no_shared_sv0 = 0.0;
							double sum_pt_sig_shared_sv1 = 0.0;
							double sum_pt_sig_shared_sv0 = 0.0;
							for (unsigned int j = 0; j < sv1_sharedjet_1_which_trk_idx.size(); j++) {
								idx = sv1_sharedjet_1_which_trk_idx[j];
								sum_pt_shared_sv1 = sum_pt_shared_sv1 + sv1.track_pt(idx);
								sum_pt_sig_shared_sv1 = sum_pt_sig_shared_sv1 + (sv1.track_pt(idx) / sv1.track_pt_err[idx]);
							}
							for (unsigned int j = 0; j < sv1_sharedjet_1_which_no_trk_idx.size(); j++) {
								idx = sv1_sharedjet_1_which_no_trk_idx[j];
								sum_pt_no_shared_sv1 = sum_pt_no_shared_sv1 + sv1.track_pt(idx);

							}
							for (unsigned int j = 0; j < sv0_sharedjet_1_which_trk_idx.size(); j++) {
								idx = sv0_sharedjet_1_which_trk_idx[j];
								sum_pt_shared_sv0 = sum_pt_shared_sv0 + sv0.track_pt(idx);
								sum_pt_sig_shared_sv0 = sum_pt_sig_shared_sv0 + (sv0.track_pt(idx) / sv0.track_pt_err[idx]);
							}
							for (unsigned int j = 0; j < sv0_sharedjet_1_which_no_trk_idx.size(); j++) {
								idx = sv0_sharedjet_1_which_no_trk_idx[j];
								sum_pt_no_shared_sv0 = sum_pt_no_shared_sv0 + sv0.track_pt(idx);

							}

							if ((sum_pt_shared_sv0 / (sum_pt_shared_sv0 + sum_pt_no_shared_sv0)) >= (sum_pt_shared_sv1 / (sum_pt_shared_sv1 + sum_pt_no_shared_sv1))) {		 // sv0 is major vtx	sv1 got removed if equal 


								if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
									if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
										h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx.size(),w);

									}
								}
								else {
									if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
										h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx.size(),w);

									}

								}

							}
							else {

								if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
									if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {
										h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx.size(),w);

									}
								}
								else {
									if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {
										h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx.size(),w);

									}

								}
							}

							if (((sum_pt_shared_sv0 * sv0_track_which_idx.size()) / ((sum_pt_shared_sv0 + sum_pt_no_shared_sv0) * sv0_sharedjet_1_which_trk_idx.size())) >= ((sum_pt_shared_sv1 * sv1_track_which_idx.size()) / ((sum_pt_shared_sv1 + sum_pt_no_shared_sv1) * sv1_sharedjet_1_which_trk_idx.size()))) {		 // sv0 is major vtx	sv1 got removed if equal 


								if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
									if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
										h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx.size(),w);

									}

								}
								else {
									if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
										h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx.size(),w);

									}
								}

							}
							else {

								if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
									if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {
										h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx.size(),w);

									}
								}
								else {
									if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {
										h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx.size(),w);

									}
								}
							}

							if ((sum_pt_shared_sv0 / (sv0_sharedjet_1_which_trk_idx.size())) >= (sum_pt_shared_sv1 / (sv1_sharedjet_1_which_trk_idx.size()))) {		 // sv0 is major vtx	sv1 got removed if equal 


								if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
									if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
										h_2D_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx.size(),w);

									}

								}
								else {
									if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
										h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx.size(),w);

									}

								}

							}
							else {

								if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
									if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {
										h_2D_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx.size(),w);

									}
								}
								else {
									if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {
										h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx.size(),w);

									}
								}
							}

							if ((sum_pt_shared_sv0) >= (sum_pt_shared_sv1)) {		 // sv0 is major vtx	sv1 got removed if equal 


								if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
									if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
										h_2D_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx.size(),w);

									}

								}
								else {
									if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
										h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx_no_shared_track.size(), sv0_track_which_idx.size(),w);

									}
								}

							}
							else {

								if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
									if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {
										h_2D_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx.size(),w);

									}
								}
								else {
									if (nsharedjet_tracks_sv1[i] > nsharedjet_tracks_sv0[i]) {
										h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx_no_shared_track.size(), sv1_track_which_idx.size(),w);

									}
								}
							}

							if (nsharedjet_tracks_sv0[i] >= nsharedjet_tracks_sv1[i]) {
								if ((sum_pt_shared_sv0 / (sum_pt_shared_sv0 + sum_pt_no_shared_sv0)) >= (sum_pt_shared_sv1 / (sum_pt_shared_sv1 + sum_pt_no_shared_sv1))){
									if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

									}

								}
								else {
									if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);

									}
								}

								if (((sum_pt_shared_sv0 * sv0_track_which_idx.size()) / ((sum_pt_shared_sv0 + sum_pt_no_shared_sv0) * sv0_sharedjet_1_which_trk_idx.size())) >= ((sum_pt_shared_sv1 * sv1_track_which_idx.size()) / ((sum_pt_shared_sv1 + sum_pt_no_shared_sv1) * sv1_sharedjet_1_which_trk_idx.size()))) {
									if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

									}

								}
								else {
									if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);

									}
								}

								if ((sum_pt_shared_sv0 / (sv0_sharedjet_1_which_trk_idx.size())) >= (sum_pt_shared_sv1 / (sv1_sharedjet_1_which_trk_idx.size()))) {
									if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

									}

								}
								else {
									if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);

									}
								}

								if ((sum_pt_shared_sv0) >= (sum_pt_shared_sv1)) {
									if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx_no_shared_track.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx_no_shared_track.size(),w);

									}

								}
								else {
									if ((sv0_track_which_idx.size() >= 5) && (sv1_track_which_idx.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);

									}
								}
								


							}
							else {

								if ((sum_pt_shared_sv1 / (sum_pt_shared_sv1 + sum_pt_no_shared_sv1)) >= (sum_pt_shared_sv0 / (sum_pt_shared_sv0 + sum_pt_no_shared_sv0))) {
									if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);

									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);

									}

								}
								else {
									if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_ratio_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);

									}
								}

								if (((sum_pt_shared_sv1 * sv1_track_which_idx.size()) / ((sum_pt_shared_sv1 + sum_pt_no_shared_sv1) * sv1_sharedjet_1_which_trk_idx.size())) >= ((sum_pt_shared_sv0 * sv0_track_which_idx.size()) / ((sum_pt_shared_sv0 + sum_pt_no_shared_sv0) * sv0_sharedjet_1_which_trk_idx.size()))) {
									if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);

									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);

									}

								}
								else {
									if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv0_track_which_idx.size(), sv1_track_which_idx.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_ratio_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);

									}
								}

								if ((sum_pt_shared_sv1 / (sv1_sharedjet_1_which_trk_idx.size())) >= (sum_pt_shared_sv0 / (sv0_sharedjet_1_which_trk_idx.size()))) {
									if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);

									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);

									}

								}
								else {
									if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_avg_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);

									}
								}

								if ((sum_pt_shared_sv1) >= (sum_pt_shared_sv0)) {
									if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx_no_shared_track.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);

									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx_no_shared_track.size(),w);

									}

								}
								else {
									if ((sv1_track_which_idx.size() >= 5) && (sv0_track_which_idx.size() >= 5)) {
										h_2D_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);
									}
									else {
										h_2D_poor_sv_tracks_no_minor_less_sum_pt_shared_tracks_nshj1_large_nsv2->Fill(sv1_track_which_idx.size(), sv0_track_which_idx.size(),w);

									}
								}

							}

							

							


						}

						//if (ratio_ntracks_nsv2 > 1) {		// study sample with all sh jets

							
                                                       
							if ((dphi_large_sv0_sharedjet > dphi_large_sv1_sharedjet)) {		 // sv1 next to a shared jet and sv0 is away


								// fig2  && semi-fig2
									
									if (nsharedjet_tracks_sv0[i] < nsharedjet_tracks_sv1[i]) {
                                                                                
										h_less_dphi_more_shared_tracks_large_nsv2_all_fig2->Fill(int(1), w);
										h_pt_shared_jets_large_nsv2_fig2->Fill(mevent->jet_pt[jet_index], w);
										h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_fig2->Fill(avg_sv1_non_shared_jet_pt, mevent->jet_pt[jet_index]);
										h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_fig2->Fill(avg_sv0_non_shared_jet_pt, mevent->jet_pt[jet_index]);
										h_dPhi_sv_large_nsv2_fig2->Fill(double(fabs(reco::deltaPhi(phi0, phi1))), w);
										h_dEta_sv_large_nsv2_fig2->Fill(double(fabs(reco::deltaPhi(eta0, eta1))), w);
										h_2D_dEta_dPhi_sv_large_nsv2_fig2->Fill(double(fabs(reco::deltaPhi(eta0, eta1))), double(fabs(reco::deltaPhi(phi0, phi1))));


										double ratio_ntracks_nsv2_fig2 = nsharedjet_tracks_sv1[i] / nsharedjet_tracks_sv0[i];
										h_ratio_ntracks_large_nsv2_shared_jets_fig2->Fill(ratio_ntracks_nsv2_fig2);

										
										std::vector<int> sv1_sharedjet_i_which_trk_idx = sv1_sharedjet_which_idx[i];
										std::vector<int> sv1_sharedjet_i_which_no_trk_idx = sv1_sharedjet_which_no_trk_idx[i];
										std::vector<int> sv0_sharedjet_i_which_trk_idx = sv0_sharedjet_which_idx[i];
										std::vector<int> sv0_sharedjet_i_which_no_trk_idx = sv0_sharedjet_which_no_trk_idx[i];
										double sum_pt_good = 0;
										double sum_pt_sig_good = 0;
										for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {
											int track_idx = sv1_sharedjet_i_which_trk_idx[j];
											double absdelta_sv1_track = double(fabs(reco::deltaPhi(phi1, sv1.track_phi[track_idx]))); //phi1
											double absdelta_sv0_track = double(fabs(reco::deltaPhi(phi0, sv1.track_phi[track_idx]))); //phi0

											double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx])));
											h_absdeltaphi_large_jet_shared_tracks_nsv2_fig2->Fill(absdelta_jet_track, w);
											h_dPhi_good_tracks_its_sv_large_fig2->Fill(absdelta_sv1_track, w);
											h_dPhi_good_tracks_other_sv_large_fig2->Fill(absdelta_sv0_track, w);
											h_pt_good_shared_tracks_large_nsv2_fig2->Fill(sv1.track_pt(track_idx), w);
											h_dxy_err_good_shared_tracks_large_nsv2_fig2->Fill(sv1.track_dxy_err(track_idx), w);
											h_dxy_sig_good_shared_tracks_large_nsv2_fig2->Fill(sv1.track_dxy[track_idx] / sv1.track_dxy_err(track_idx), w);
											sum_pt_good = sum_pt_good + sv1.track_pt(track_idx);
											sum_pt_sig_good = sum_pt_sig_good + (sv1.track_pt(track_idx) / sv1.track_pt_err[track_idx]);

											for (int k = 0; k < nsharedjet_tracks_sv0[i]; k++) {
												int track_idx_poor = sv0_sharedjet_i_which_trk_idx[k];
												h_2D_poor_dxy_sig_good_dxy_sig_large_sv_fig2->Fill(sv0.track_dxy[track_idx_poor] / sv0.track_dxy_err(track_idx_poor), sv1.track_dxy[track_idx] / sv1.track_dxy_err(track_idx));
												h_2D_poor_dxy_err_good_dxy_err_large_sv_fig2->Fill(sv0.track_dxy_err(track_idx_poor), sv1.track_dxy_err(track_idx));
												h_2D_poor_dz_err_good_dz_err_large_sv_fig2->Fill(sv0.track_dz_err(track_idx_poor), sv1.track_dz_err(track_idx));
												h_2D_poor_dz_sig_good_dz_sig_large_sv_fig2->Fill(sv0.track_dz[track_idx_poor] / sv0.track_dz_err(track_idx_poor), sv1.track_dz[track_idx] / sv1.track_dz_err(track_idx));
												
												
												h_2D_poor_pt_good_pt_large_sv_fig2->Fill(sv0.track_pt(track_idx_poor), sv1.track_pt(track_idx));

											}

										}
										
										
										

										
										
										double sum_pt_poor = 0;
										double sum_pt_sig_poor = 0;
										for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {
											int track_idx = sv0_sharedjet_i_which_trk_idx[j];
											double absdelta_sv0_track = double(fabs(reco::deltaPhi(phi0, sv0.track_phi[track_idx]))); //phi0
											double absdelta_sv1_track = double(fabs(reco::deltaPhi(phi1, sv0.track_phi[track_idx]))); //phi0
											double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx])));                                                                     h_absdeltaphi_large_jet_shared_tracks_nsv2_fig2->Fill(absdelta_jet_track, w);
											h_pt_poor_shared_tracks_large_nsv2_fig2->Fill(sv0.track_pt(track_idx), w);
											h_dPhi_poor_tracks_its_sv_large_fig2->Fill(absdelta_sv0_track, w);
											h_dPhi_poor_tracks_other_sv_large_fig2->Fill(absdelta_sv1_track, w);
											h_dxy_err_poor_shared_tracks_large_nsv2_fig2->Fill(sv0.track_dxy_err(track_idx), w);
											h_dxy_sig_poor_shared_tracks_large_nsv2_fig2->Fill(sv0.track_dxy[track_idx] / sv0.track_dxy_err(track_idx), w);
											sum_pt_poor = sum_pt_poor + sv0.track_pt(track_idx);
											sum_pt_sig_poor = sum_pt_sig_poor + (sv0.track_pt(track_idx) / sv0.track_pt_err[track_idx]);

								
										}


										h_2D_poor_avg_pt_good_avg_pt_large_sv_fig2->Fill(sum_pt_poor / nsharedjet_tracks_sv0[i], sum_pt_good / nsharedjet_tracks_sv1[i]);
										h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_fig2->Fill(sum_pt_sig_poor / nsharedjet_tracks_sv0[i], sum_pt_sig_good / nsharedjet_tracks_sv1[i]);
										//h_diff_pT_avg_major_minor_sv_nsv2_all_fig2->Fill((sum_pt_good / nsharedjet_tracks_sv1[i]) - (sum_pt_poor / nsharedjet_tracks_sv0[i]), w);
										//h_diff_pT_sum_major_minor_sv_nsv2_all_fig2->Fill(sum_pt_good - sum_pt_poor, w);

										double sum_pt_non_shared_poor = 0;
										for (unsigned int j = 0; j < sv0_sharedjet_i_which_no_trk_idx.size(); j++) {
											idx = sv0_sharedjet_i_which_no_trk_idx[j];
											sum_pt_non_shared_poor = sum_pt_non_shared_poor + sv0.track_pt(idx);
											}
										double sum_pt_non_shared_good = 0;
										for (unsigned int j = 0; j < sv1_sharedjet_i_which_no_trk_idx.size(); j++) {
											idx = sv1_sharedjet_i_which_no_trk_idx[j];
											sum_pt_non_shared_good = sum_pt_non_shared_good + sv1.track_pt(idx);
										}

										
										h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_fig2->Fill(sum_pt_poor/ (sum_pt_poor+sum_pt_non_shared_poor), sum_pt_good/(sum_pt_good+sum_pt_non_shared_good));
										h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_fig2->Fill((sum_pt_poor* sv0_track_which_idx.size()) / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv0[i]), (sum_pt_good* sv1_track_which_idx.size()) / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv1[i]));
										//h_diff_ratio_pT_avg_major_minor_sv_nsv2_all_fig2->Fill(((sum_pt_good* sv1_track_which_idx.size()) / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv1[i]))-((sum_pt_poor * sv0_track_which_idx.size()) / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv0[i])), w);
										//h_diff_ratio_pT_sum_major_minor_sv_nsv2_all_fig2->Fill((sum_pt_good / (sum_pt_good + sum_pt_non_shared_good)) - (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor)), w);

										//h_ratio_diff_pT_sum_major_minor_sv_nsv2_all_fig2->Fill(((sum_pt_good) - (sum_pt_poor))/((sum_pt_good)+(sum_pt_poor)), w);
										double avg_pt_good = sum_pt_good / nsharedjet_tracks_sv1[i];
										double avg_pt_poor = sum_pt_poor / nsharedjet_tracks_sv0[i];
										//h_ratio_diff_pT_avg_major_minor_sv_nsv2_all_fig2->Fill(((avg_pt_good)-(avg_pt_poor)) / ((avg_pt_good)+(avg_pt_poor)), w);


										if (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) < 0 || sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) < 0 || sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) == 1 || sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) == 1) {
											std::cout << "found minor's ratio sum pT = " << sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) << "and major's ratio sum pT = " << sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) << std::endl;
											
											std::cout << "minor tracks = " << sv0_track_which_idx.size() << ", shared ones = " << nsharedjet_tracks_sv0[i] << ", not shared ones = " << sv0_sharedjet_i_which_no_trk_idx.size() << std::endl;
											for (unsigned int j = 0; j < sv0_sharedjet_i_which_trk_idx.size(); j++) {
												idx = sv0_sharedjet_i_which_trk_idx[j];
												std::cout << "minor's track pT = " << sv0.track_pt(idx) << std::endl;
											}
											
											std::cout << "major tracks = " << sv1_track_which_idx.size() << ", shared ones = " << nsharedjet_tracks_sv1[i] << ", not shared ones = " << sv1_sharedjet_i_which_no_trk_idx.size() << std::endl;
											for (unsigned int j = 0; j < sv1_sharedjet_i_which_trk_idx.size(); j++) {
												idx = sv1_sharedjet_i_which_trk_idx[j];
												std::cout << "major's track pT = " << sv1.track_pt(idx) << std::endl;
											}
											std::cout << ", sum pT minor tracks = " << sum_pt_poor + sum_pt_non_shared_poor << "sum pT major tracks = " << sum_pt_good + sum_pt_non_shared_good << std::endl;

											
										}

										if ((sum_pt_good / (sum_pt_good + sum_pt_non_shared_good)) >= (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor))) {
											h_more_ratio_sum_major_sv_nsv2_all_fig2->Fill(int(1), w);
											h_more_ratio_sum_major_sv_nsv2_fig2->Fill(int(1), w);
										}
										else {
											h_more_ratio_sum_major_sv_nsv2_all_fig2->Fill(int(0), w);
											h_more_ratio_sum_major_sv_nsv2_fig2->Fill(int(0), w);

										}

										if ((sum_pt_good * sv1_track_which_idx.size()/ ((sum_pt_good + sum_pt_non_shared_good)* nsharedjet_tracks_sv1[i])) >= (sum_pt_poor * sv0_track_which_idx.size() / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv0[i]))) {
											h_more_ratio_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
											h_more_ratio_avg_major_sv_nsv2_fig2->Fill(int(1), w);
										}
										else {
											h_more_ratio_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
											h_more_ratio_avg_major_sv_nsv2_fig2->Fill(int(0), w);

										}

										if ((sum_pt_good/ (nsharedjet_tracks_sv1[i])) >= (sum_pt_poor/ (nsharedjet_tracks_sv0[i]))) {
											h_more_pT_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
											h_more_pT_avg_major_sv_nsv2_fig2->Fill(int(1), w);
										}
										else {
											h_more_pT_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
											h_more_pT_avg_major_sv_nsv2_fig2->Fill(int(0), w);

										}

										if ((sum_pt_sig_good / (nsharedjet_tracks_sv1[i])) >= (sum_pt_sig_poor / (nsharedjet_tracks_sv0[i]))) {
											h_more_pT_sig_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
											h_more_pT_avg_major_sv_nsv2_fig2->Fill(int(1), w);
										}
										else {
											h_more_pT_sig_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
											h_more_pT_sig_avg_major_sv_nsv2_fig2->Fill(int(0), w);

										}


										
										if (nsharedjet_tracks_sv0[i] == 1) {


											int one_poor_sv0_track_idx = sv0_sharedjet_i_which_trk_idx[0];
											AlgebraicVector3 mom_one_poor(sv0.track_px[one_poor_sv0_track_idx], sv0.track_py[one_poor_sv0_track_idx], sv0.track_pz[one_poor_sv0_track_idx]);
											AlgebraicVector3 ref_one_poor(sv0.track_vx[one_poor_sv0_track_idx], sv0.track_vy[one_poor_sv0_track_idx], sv0.track_vz[one_poor_sv0_track_idx]);

											Measurement1D miss_dist_one_poor_other = miss_dist(sv1, ref_one_poor, mom_one_poor);
											Measurement1D miss_dist_one_poor_its = miss_dist(sv0, ref_one_poor, mom_one_poor);

											h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_fig2->Fill(miss_dist_one_poor_its.significance(), miss_dist_one_poor_other.significance());

											double sum_good_miss_dist_sig_other = 0;
											double sum_good_miss_dist_sig_its = 0;

											for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {

												int more_good_sv1_track_idx = sv1_sharedjet_i_which_trk_idx[j];
												AlgebraicVector3 mom_more_good(sv1.track_px[more_good_sv1_track_idx], sv1.track_py[more_good_sv1_track_idx], sv1.track_pz[more_good_sv1_track_idx]);
												AlgebraicVector3 ref_more_good(sv1.track_vx[more_good_sv1_track_idx], sv1.track_vy[more_good_sv1_track_idx], sv1.track_vz[more_good_sv1_track_idx]);

												Measurement1D miss_dist_more_good_other = miss_dist(sv0, ref_more_good, mom_more_good);
												Measurement1D miss_dist_more_good_its = miss_dist(sv1, ref_more_good, mom_more_good);

												sum_good_miss_dist_sig_other = sum_good_miss_dist_sig_other + miss_dist_more_good_other.significance();
												sum_good_miss_dist_sig_its = sum_good_miss_dist_sig_its + miss_dist_more_good_its.significance();



											}


											h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_fig2->Fill(miss_dist_one_poor_its.significance(), sum_good_miss_dist_sig_its / nsharedjet_tracks_sv1[i]);
											h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_fig2->Fill(miss_dist_one_poor_other.significance(), sum_good_miss_dist_sig_other / nsharedjet_tracks_sv1[i]);

										}

										

										
										
										h_vertex_chi2dof_absdeltaphi0_large_nsv2_fig2->Fill(sv0.chi2dof(), w);
									}
									else
									{
									h_less_dphi_more_shared_tracks_large_nsv2_all_fig2->Fill(int(0), w);		   // sv1 next to a shared jet and sv0 is away
									h_pt_shared_jets_large_nsv2_semi_fig2->Fill(mevent->jet_pt[jet_index],w);
									h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2->Fill(avg_sv0_non_shared_jet_pt, mevent->jet_pt[jet_index]);
									h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2->Fill(avg_sv1_non_shared_jet_pt, mevent->jet_pt[jet_index]);
									h_dPhi_sv_large_nsv2_semi_fig2->Fill(double(fabs(reco::deltaPhi(phi0, phi1))), w);
									h_dEta_sv_large_nsv2_semi_fig2->Fill(double(fabs(reco::deltaPhi(eta0, eta1))), w);
									h_2D_dEta_dPhi_sv_large_nsv2_semi_fig2->Fill(double(fabs(reco::deltaPhi(eta0, eta1))), double(fabs(reco::deltaPhi(phi0, phi1))));

									

									double ratio_ntracks_nsv2_semi_fig2 = nsharedjet_tracks_sv0[i] / nsharedjet_tracks_sv1[i];
									h_ratio_ntracks_large_nsv2_shared_jets_semi_fig2->Fill(ratio_ntracks_nsv2_semi_fig2);

									std::vector<double> absdeltaphi_max_jet_shared_tracks;
									std::vector<int> sv0_sharedjet_i_which_trk_idx = sv0_sharedjet_which_idx[i];
									std::vector<int> sv0_sharedjet_i_which_no_trk_idx = sv0_sharedjet_which_no_trk_idx[i];
									std::vector<int> sv1_sharedjet_i_which_trk_idx = sv1_sharedjet_which_idx[i];
									std::vector<int> sv1_sharedjet_i_which_no_trk_idx = sv1_sharedjet_which_no_trk_idx[i];
									double sum_pt_good = 0;
									double sum_pt_sig_good = 0;
									for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {
										int track_idx = sv0_sharedjet_i_which_trk_idx[j];
										double absdelta_sv0_track = double(fabs(reco::deltaPhi(phi0, sv0.track_phi[track_idx]))); //phi1
										double absdelta_sv1_track = double(fabs(reco::deltaPhi(phi1, sv0.track_phi[track_idx]))); //phi0

										double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx])));
										h_absdeltaphi_large_jet_shared_tracks_nsv2_semi_fig2->Fill(absdelta_jet_track, w);
										h_dPhi_good_tracks_its_sv_large_semi_fig2->Fill(absdelta_sv0_track, w);
										h_dPhi_good_tracks_other_sv_large_semi_fig2->Fill(absdelta_sv1_track, w);
										h_pt_good_shared_tracks_large_nsv2_semi_fig2->Fill(sv0.track_pt(track_idx), w);
										h_dxy_err_good_shared_tracks_large_nsv2_semi_fig2->Fill(sv0.track_dxy_err(track_idx), w);
										h_dxy_sig_good_shared_tracks_large_nsv2_semi_fig2->Fill(sv0.track_dxy[track_idx] / sv0.track_dxy_err(track_idx), w);
										absdeltaphi_max_jet_shared_tracks.push_back(absdelta_jet_track);
										sum_pt_good = sum_pt_good + sv0.track_pt(track_idx);
										sum_pt_sig_good = sum_pt_sig_good + (sv0.track_pt(track_idx)/ sv0.track_pt_err[track_idx]);

										for (int k = 0; k < nsharedjet_tracks_sv1[i]; k++) {
											int track_idx_poor = sv1_sharedjet_i_which_trk_idx[k];
											h_2D_poor_dxy_sig_good_dxy_sig_large_sv_semi_fig2->Fill(sv1.track_dxy[track_idx_poor] / sv1.track_dxy_err(track_idx_poor), sv0.track_dxy[track_idx] / sv0.track_dxy_err(track_idx));
											h_2D_poor_dxy_err_good_dxy_err_large_sv_semi_fig2->Fill(sv1.track_dxy_err(track_idx_poor), sv0.track_dxy_err(track_idx));
											h_2D_poor_dz_err_good_dz_err_large_sv_semi_fig2->Fill(sv1.track_dz_err(track_idx_poor), sv0.track_dz_err(track_idx));
											h_2D_poor_dz_sig_good_dz_sig_large_sv_semi_fig2->Fill(sv1.track_dz[track_idx_poor] / sv1.track_dz_err(track_idx_poor), sv0.track_dz[track_idx] / sv0.track_dz_err(track_idx));

											
											
											h_2D_poor_pt_good_pt_large_sv_semi_fig2->Fill(sv1.track_pt(track_idx_poor), sv0.track_pt(track_idx));

										}

									}
									
									double max_dphi_jet_track = *std::max_element(absdeltaphi_max_jet_shared_tracks.begin(), absdeltaphi_max_jet_shared_tracks.end());
									int max_dphi_jet_track_idx = std::max_element(absdeltaphi_max_jet_shared_tracks.begin(), absdeltaphi_max_jet_shared_tracks.end()) - absdeltaphi_max_jet_shared_tracks.begin();

									int max_sv0_track_idx = sv0_sharedjet_i_which_trk_idx[max_dphi_jet_track_idx];
									

									double sum_pt_poor = 0;
									double sum_pt_sig_poor = 0;
									
									std::vector<double> absdeltaphi_min_jet_shared_tracks;
									for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {
										int track_idx = sv1_sharedjet_i_which_trk_idx[j];
										double absdelta_sv1_track = double(fabs(reco::deltaPhi(phi1, sv1.track_phi[track_idx]))); //phi1
										double absdelta_sv0_track = double(fabs(reco::deltaPhi(phi0, sv1.track_phi[track_idx]))); //phi0
										double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx])));                                                                     h_absdeltaphi_large_jet_shared_tracks_nsv2_semi_fig2->Fill(absdelta_jet_track, w);
										h_pt_poor_shared_tracks_large_nsv2_semi_fig2->Fill(sv1.track_pt(track_idx), w);
										h_dPhi_poor_tracks_its_sv_large_semi_fig2->Fill(absdelta_sv1_track, w);
										h_dPhi_poor_tracks_other_sv_large_semi_fig2->Fill(absdelta_sv0_track, w);
										h_dxy_err_poor_shared_tracks_large_nsv2_semi_fig2->Fill(sv1.track_dxy_err(track_idx), w);
										h_dxy_sig_poor_shared_tracks_large_nsv2_semi_fig2->Fill(sv1.track_dxy[track_idx] / sv1.track_dxy_err(track_idx), w);

										absdeltaphi_min_jet_shared_tracks.push_back(absdelta_jet_track);
										sum_pt_poor = sum_pt_poor + sv1.track_pt(track_idx);
										sum_pt_sig_poor = sum_pt_sig_poor + (sv1.track_pt(track_idx) / sv1.track_pt_err[track_idx]);
									}

									h_2D_poor_avg_pt_good_avg_pt_large_sv_semi_fig2->Fill(sum_pt_poor / nsharedjet_tracks_sv1[i], sum_pt_good / nsharedjet_tracks_sv0[i]);
									h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_semi_fig2->Fill(sum_pt_sig_poor / nsharedjet_tracks_sv1[i], sum_pt_sig_good / nsharedjet_tracks_sv0[i]);
									//h_diff_pT_avg_major_minor_sv_nsv2_all_fig2->Fill((sum_pt_good / nsharedjet_tracks_sv0[i]) - (sum_pt_poor / nsharedjet_tracks_sv1[i]), w);
									//h_diff_pT_sum_major_minor_sv_nsv2_all_fig2->Fill(sum_pt_good - sum_pt_poor, w);

									double sum_pt_non_shared_poor = 0;
									for (unsigned int j = 0; j < sv1_sharedjet_i_which_no_trk_idx.size(); j++) {
										idx = sv1_sharedjet_i_which_no_trk_idx[j];
										sum_pt_non_shared_poor = sum_pt_non_shared_poor + sv1.track_pt(idx);
									}
									double sum_pt_non_shared_good = 0;
									for (unsigned int j = 0; j < sv0_sharedjet_i_which_no_trk_idx.size(); j++) {
										idx = sv0_sharedjet_i_which_no_trk_idx[j];
										sum_pt_non_shared_good = sum_pt_non_shared_good + sv0.track_pt(idx);
									}
									h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_semi_fig2->Fill(sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor), sum_pt_good / (sum_pt_good + sum_pt_non_shared_good));
									h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_semi_fig2->Fill((sum_pt_poor* sv1_track_which_idx.size()) / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv1[i]), (sum_pt_good* sv0_track_which_idx.size()) / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv0[i]));
									//h_diff_ratio_pT_avg_major_minor_sv_nsv2_all_fig2->Fill(((sum_pt_good* sv0_track_which_idx.size()) / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv0[i])) - ((sum_pt_poor * sv1_track_which_idx.size()) / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv1[i])), w);
									//h_diff_ratio_pT_sum_major_minor_sv_nsv2_all_fig2->Fill((sum_pt_good / (sum_pt_good + sum_pt_non_shared_good)) - (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor)), w);

									//h_ratio_diff_pT_sum_major_minor_sv_nsv2_all_fig2->Fill(((sum_pt_good)-(sum_pt_poor)) / ((sum_pt_good)+(sum_pt_poor)), w);
									double avg_pt_good = sum_pt_good / nsharedjet_tracks_sv0[i];
									double avg_pt_poor = sum_pt_poor / nsharedjet_tracks_sv1[i];
									//h_ratio_diff_pT_avg_major_minor_sv_nsv2_all_fig2->Fill(((avg_pt_good)-(avg_pt_poor)) / ((avg_pt_good)+(avg_pt_poor)), w);

									if (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) < 0 || sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) < 0 || sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) == 1 || sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) == 1) {
										std::cout << "found minor's ratio sum pT = " << sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) << "and major's ratio sum pT = " << sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) << std::endl;

										

										std::cout << "minor tracks = " << sv1_track_which_idx.size() << ", shared ones = " << nsharedjet_tracks_sv1[i] << ", not shared ones = " << sv1_sharedjet_i_which_no_trk_idx.size() << std::endl;
										for (unsigned int j = 0; j < sv1_sharedjet_i_which_trk_idx.size(); j++) {
											idx = sv1_sharedjet_i_which_trk_idx[j];
											std::cout << "minor's track pT = " << sv1.track_pt(idx) << std::endl;
										}

										std::cout << "major tracks = " << sv0_track_which_idx.size() << ", shared ones = " << nsharedjet_tracks_sv0[i] << ", not shared ones = " << sv0_sharedjet_i_which_no_trk_idx.size() << std::endl;
										for (unsigned int j = 0; j < sv0_sharedjet_i_which_trk_idx.size(); j++) {
											idx = sv0_sharedjet_i_which_trk_idx[j];
											std::cout << "major's track pT = " << sv0.track_pt(idx) << std::endl;
										}
										std::cout << ", sum pT minor tracks = " << sum_pt_poor + sum_pt_non_shared_poor << "sum pT major tracks = " << sum_pt_good + sum_pt_non_shared_good << std::endl;

										
									}

									if ((sum_pt_good / (sum_pt_good + sum_pt_non_shared_good)) >= (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor))) {
										h_more_ratio_sum_major_sv_nsv2_all_fig2->Fill(int(1), w);
										h_more_ratio_sum_major_sv_nsv2_semi_fig2->Fill(int(1), w);
									}
									else {
										h_more_ratio_sum_major_sv_nsv2_all_fig2->Fill(int(0), w);
										h_more_ratio_sum_major_sv_nsv2_semi_fig2->Fill(int(0), w);

									}

									if ((sum_pt_good * sv0_track_which_idx.size() / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv0[i])) >= (sum_pt_poor * sv1_track_which_idx.size() / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv1[i]))) {
										h_more_ratio_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
										h_more_ratio_avg_major_sv_nsv2_semi_fig2->Fill(int(1), w);
									}
									else {
										h_more_ratio_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
										h_more_ratio_avg_major_sv_nsv2_semi_fig2->Fill(int(0), w);

									}

									if ((sum_pt_good / (nsharedjet_tracks_sv0[i])) >= (sum_pt_poor / (nsharedjet_tracks_sv1[i]))) {
										h_more_pT_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
										h_more_pT_avg_major_sv_nsv2_semi_fig2->Fill(int(1), w);
									}
									else {
										h_more_pT_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
										h_more_pT_avg_major_sv_nsv2_semi_fig2->Fill(int(0), w);

									}

									if ((sum_pt_sig_good / (nsharedjet_tracks_sv0[i])) >= (sum_pt_sig_poor / (nsharedjet_tracks_sv1[i]))) {
										h_more_pT_sig_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
										h_more_pT_sig_avg_major_sv_nsv2_semi_fig2->Fill(int(1), w);
									}
									else {
										h_more_pT_sig_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
										h_more_pT_sig_avg_major_sv_nsv2_semi_fig2->Fill(int(0), w);

									}

									
									
									double min_dphi_jet_track = *std::min_element(absdeltaphi_min_jet_shared_tracks.begin(), absdeltaphi_min_jet_shared_tracks.end());
									int min_dphi_jet_track_idx = std::min_element(absdeltaphi_min_jet_shared_tracks.begin(), absdeltaphi_min_jet_shared_tracks.end()) - absdeltaphi_min_jet_shared_tracks.begin();

									int min_sv1_track_idx = sv1_sharedjet_i_which_trk_idx[min_dphi_jet_track_idx];
									




									if (nsharedjet_tracks_sv1[i] == 1) {


										int one_poor_sv1_track_idx = sv1_sharedjet_i_which_trk_idx[0];
										AlgebraicVector3 mom_one_poor(sv1.track_px[one_poor_sv1_track_idx], sv1.track_py[one_poor_sv1_track_idx], sv1.track_pz[one_poor_sv1_track_idx]);
										AlgebraicVector3 ref_one_poor(sv1.track_vx[one_poor_sv1_track_idx], sv1.track_vy[one_poor_sv1_track_idx], sv1.track_vz[one_poor_sv1_track_idx]);

										Measurement1D miss_dist_one_poor_other = miss_dist(sv0, ref_one_poor, mom_one_poor);
										Measurement1D miss_dist_one_poor_its = miss_dist(sv1, ref_one_poor, mom_one_poor);

										h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_semi_fig2->Fill(miss_dist_one_poor_its.significance(), miss_dist_one_poor_other.significance());

										double sum_good_miss_dist_sig_other = 0;
										double sum_good_miss_dist_sig_its = 0;

										for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {

											int more_good_sv0_track_idx = sv0_sharedjet_i_which_trk_idx[j];
											AlgebraicVector3 mom_more_good(sv0.track_px[more_good_sv0_track_idx], sv0.track_py[more_good_sv0_track_idx], sv0.track_pz[more_good_sv0_track_idx]);
											AlgebraicVector3 ref_more_good(sv0.track_vx[more_good_sv0_track_idx], sv0.track_vy[more_good_sv0_track_idx], sv0.track_vz[more_good_sv0_track_idx]);

											Measurement1D miss_dist_more_good_other = miss_dist(sv1, ref_more_good, mom_more_good);
											Measurement1D miss_dist_more_good_its = miss_dist(sv0, ref_more_good, mom_more_good);

											sum_good_miss_dist_sig_other = sum_good_miss_dist_sig_other + miss_dist_more_good_other.significance();
											sum_good_miss_dist_sig_its = sum_good_miss_dist_sig_its + miss_dist_more_good_its.significance();



										}


										h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_semi_fig2->Fill(miss_dist_one_poor_its.significance(), sum_good_miss_dist_sig_its / nsharedjet_tracks_sv0[i]);
										h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_semi_fig2->Fill(miss_dist_one_poor_other.significance(), sum_good_miss_dist_sig_other / nsharedjet_tracks_sv0[i]);

									}


									

									h_2D_pt_dxy_err_absdeltaphi0_large_sv_semi_fig2->Fill(sv0.track_pt(max_sv0_track_idx), sv0.track_dxy_err(max_sv0_track_idx));
									h_2D_pt_dxy_err_absdeltaphi1_large_sv_semi_fig2->Fill(sv1.track_pt(min_sv1_track_idx), sv1.track_dxy_err(min_sv1_track_idx));

									h_vertex_chi2dof_absdeltaphi0_large_nsv2_semi_fig2->Fill(sv0.chi2dof(), w);
									}



							}
							
							if (dphi_large_sv1_sharedjet >= dphi_large_sv0_sharedjet) {	  // sv0 next to a shared jet and sv1 is away

								// fig2	&& semi-fig2

									if (nsharedjet_tracks_sv1[i] <= nsharedjet_tracks_sv0[i]) {
										h_less_dphi_more_shared_tracks_large_nsv2_all_fig2->Fill(int(1), w);
										h_pt_shared_jets_large_nsv2_fig2->Fill(mevent->jet_pt[jet_index], w);
										h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_fig2->Fill(avg_sv0_non_shared_jet_pt, mevent->jet_pt[jet_index]);
										h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_fig2->Fill(avg_sv1_non_shared_jet_pt, mevent->jet_pt[jet_index]);
										h_dPhi_sv_large_nsv2_fig2->Fill(double(fabs(reco::deltaPhi(phi0, phi1))), w);
										h_dEta_sv_large_nsv2_fig2->Fill(double(fabs(reco::deltaPhi(eta0, eta1))), w);
										h_2D_dEta_dPhi_sv_large_nsv2_fig2->Fill(double(fabs(reco::deltaPhi(eta0, eta1))), double(fabs(reco::deltaPhi(phi0, phi1))));
										

										double ratio_ntracks_nsv2_fig2 = nsharedjet_tracks_sv0[i] / nsharedjet_tracks_sv1[i];
										h_ratio_ntracks_large_nsv2_shared_jets_fig2->Fill(ratio_ntracks_nsv2_fig2);

										
										std::vector<int> sv0_sharedjet_i_which_trk_idx = sv0_sharedjet_which_idx[i];
										std::vector<int> sv0_sharedjet_i_which_no_trk_idx = sv0_sharedjet_which_no_trk_idx[i];
										std::vector<int> sv1_sharedjet_i_which_trk_idx = sv1_sharedjet_which_idx[i];
										std::vector<int> sv1_sharedjet_i_which_no_trk_idx = sv1_sharedjet_which_no_trk_idx[i];
										double sum_pt_good = 0;
										double sum_pt_sig_good = 0;
										for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {
											int track_idx = sv0_sharedjet_i_which_trk_idx[j];
											double absdelta_sv0_track = double(fabs(reco::deltaPhi(phi0, sv0.track_phi[track_idx]))); //phi0
											double absdelta_sv1_track = double(fabs(reco::deltaPhi(phi1, sv0.track_phi[track_idx]))); //phi1
											double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx])));
											h_absdeltaphi_large_jet_shared_tracks_nsv2_fig2->Fill(absdelta_jet_track, w);
											h_pt_good_shared_tracks_large_nsv2_fig2->Fill(sv0.track_pt(track_idx), w);
											h_dPhi_good_tracks_its_sv_large_fig2->Fill(absdelta_sv0_track, w);
											h_dPhi_good_tracks_other_sv_large_fig2->Fill(absdelta_sv1_track, w);
											h_dxy_err_good_shared_tracks_large_nsv2_fig2->Fill(sv0.track_dxy_err(track_idx), w);
											h_dxy_sig_good_shared_tracks_large_nsv2_fig2->Fill(sv0.track_dxy[track_idx] / sv0.track_dxy_err(track_idx), w);
											sum_pt_good = sum_pt_good + sv0.track_pt(track_idx);
											sum_pt_sig_good = sum_pt_sig_good + (sv0.track_pt(track_idx) / sv0.track_pt_err[track_idx]);
                     						for (int k = 0; k < nsharedjet_tracks_sv1[i]; k++) {
												int track_idx_poor = sv1_sharedjet_i_which_trk_idx[k];
												h_2D_poor_dxy_sig_good_dxy_sig_large_sv_fig2->Fill(sv1.track_dxy[track_idx_poor] / sv1.track_dxy_err(track_idx_poor), sv0.track_dxy[track_idx] / sv0.track_dxy_err(track_idx));
												h_2D_poor_dxy_err_good_dxy_err_large_sv_fig2->Fill(sv1.track_dxy_err(track_idx_poor), sv0.track_dxy_err(track_idx));
												h_2D_poor_dz_err_good_dz_err_large_sv_fig2->Fill(sv1.track_dz_err(track_idx_poor), sv0.track_dz_err(track_idx));
												h_2D_poor_dz_sig_good_dz_sig_large_sv_fig2->Fill(sv1.track_dz[track_idx_poor] / sv1.track_dz_err(track_idx_poor), sv0.track_dz[track_idx] / sv0.track_dz_err(track_idx));

												
												
												h_2D_poor_pt_good_pt_large_sv_fig2->Fill(sv1.track_pt(track_idx_poor), sv0.track_pt(track_idx));

											}
										}
										
										

										
										double sum_pt_poor = 0;
										double sum_pt_sig_poor = 0;
										for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {
											int track_idx = sv1_sharedjet_i_which_trk_idx[j];
											double absdelta_sv1_track = double(fabs(reco::deltaPhi(phi1, sv1.track_phi[track_idx]))); //phi1
											double absdelta_sv0_track = double(fabs(reco::deltaPhi(phi0, sv1.track_phi[track_idx]))); //phi1
											double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx])));
											h_absdeltaphi_large_jet_shared_tracks_nsv2_fig2->Fill(absdelta_jet_track, w);
											h_pt_poor_shared_tracks_large_nsv2_fig2->Fill(sv1.track_pt(track_idx), w);
											h_dPhi_poor_tracks_its_sv_large_fig2->Fill(absdelta_sv1_track, w);
											h_dPhi_poor_tracks_other_sv_large_fig2->Fill(absdelta_sv0_track, w);
											h_dxy_err_poor_shared_tracks_large_nsv2_fig2->Fill(sv1.track_dxy_err(track_idx), w);
											h_dxy_sig_poor_shared_tracks_large_nsv2_fig2->Fill(sv1.track_dxy[track_idx] / sv1.track_dxy_err(track_idx), w);
											sum_pt_poor = sum_pt_poor + sv1.track_pt(track_idx);
											sum_pt_sig_poor = sum_pt_sig_poor + (sv1.track_pt(track_idx) / sv1.track_pt_err[track_idx]);

										
											}
										h_2D_poor_avg_pt_good_avg_pt_large_sv_fig2->Fill(sum_pt_poor / nsharedjet_tracks_sv1[i], sum_pt_good / nsharedjet_tracks_sv0[i]);
										h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_fig2->Fill(sum_pt_sig_poor / nsharedjet_tracks_sv1[i], sum_pt_sig_good / nsharedjet_tracks_sv0[i]);
										//h_diff_pT_avg_major_minor_sv_nsv2_all_fig2->Fill((sum_pt_good / nsharedjet_tracks_sv0[i]) - (sum_pt_poor / nsharedjet_tracks_sv1[i]), w);
										//h_diff_pT_sum_major_minor_sv_nsv2_all_fig2->Fill(sum_pt_good - sum_pt_poor, w);

										double sum_pt_non_shared_poor = 0;
										for (unsigned int j = 0; j < sv1_sharedjet_i_which_no_trk_idx.size(); j++) {
											idx = sv1_sharedjet_i_which_no_trk_idx[j];
											sum_pt_non_shared_poor = sum_pt_non_shared_poor + sv1.track_pt(idx);
										}
										double sum_pt_non_shared_good = 0;
										for (unsigned int j = 0; j < sv0_sharedjet_i_which_no_trk_idx.size(); j++) {
											idx = sv0_sharedjet_i_which_no_trk_idx[j];
											sum_pt_non_shared_good = sum_pt_non_shared_good + sv0.track_pt(idx);
										}
										
										h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_fig2->Fill(sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor), sum_pt_good / (sum_pt_good + sum_pt_non_shared_good));
										h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_fig2->Fill((sum_pt_poor * sv1_track_which_idx.size()) / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv1[i]), (sum_pt_good * sv0_track_which_idx.size()) / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv0[i]));
										//h_diff_ratio_pT_avg_major_minor_sv_nsv2_all_fig2->Fill(((sum_pt_good* sv0_track_which_idx.size()) / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv0[i])) - ((sum_pt_poor * sv1_track_which_idx.size()) / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv1[i])), w);
										//h_diff_ratio_pT_sum_major_minor_sv_nsv2_all_fig2->Fill((sum_pt_good / (sum_pt_good + sum_pt_non_shared_good)) - (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor)), w);

										//h_ratio_diff_pT_sum_major_minor_sv_nsv2_all_fig2->Fill(((sum_pt_good)-(sum_pt_poor)) / ((sum_pt_good)+(sum_pt_poor)), w);
										double avg_pt_good = sum_pt_good / nsharedjet_tracks_sv0[i];
										double avg_pt_poor = sum_pt_poor / nsharedjet_tracks_sv1[i];
										//h_ratio_diff_pT_avg_major_minor_sv_nsv2_all_fig2->Fill(((avg_pt_good)-(avg_pt_poor)) / ((avg_pt_good)+(avg_pt_poor)), w);

										if (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) < 0 || sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) < 0 || sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) == 1 || sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) == 1) {
											std::cout << "found minor's ratio sum pT = " << sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) << "and major's ratio sum pT = " << sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) << std::endl;



											std::cout << "minor tracks = " << sv1_track_which_idx.size() << ", shared ones = " << nsharedjet_tracks_sv1[i] << ", not shared ones = " << sv1_sharedjet_i_which_no_trk_idx.size() << std::endl;
											for (unsigned int j = 0; j < sv1_sharedjet_i_which_trk_idx.size(); j++) {
												idx = sv1_sharedjet_i_which_trk_idx[j];
												std::cout << "minor's track pT = " << sv1.track_pt(idx) << std::endl;
											}

											std::cout << "major tracks = " << sv0_track_which_idx.size() << ", shared ones = " << nsharedjet_tracks_sv0[i] << ", not shared ones = " << sv0_sharedjet_i_which_no_trk_idx.size() << std::endl;
											for (unsigned int j = 0; j < sv0_sharedjet_i_which_trk_idx.size(); j++) {
												idx = sv0_sharedjet_i_which_trk_idx[j];
												std::cout << "major's track pT = " << sv0.track_pt(idx) << std::endl;
											}
											std::cout << ", sum pT minor tracks = " << sum_pt_poor + sum_pt_non_shared_poor << "sum pT major tracks = " << sum_pt_good + sum_pt_non_shared_good << std::endl;


										}

										if ((sum_pt_good / (sum_pt_good + sum_pt_non_shared_good)) >= (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor))) {
											h_more_ratio_sum_major_sv_nsv2_all_fig2->Fill(int(1), w);
											h_more_ratio_sum_major_sv_nsv2_fig2->Fill(int(1), w);
										}
										else {
											h_more_ratio_sum_major_sv_nsv2_all_fig2->Fill(int(0), w);
											h_more_ratio_sum_major_sv_nsv2_fig2->Fill(int(0), w);

										}

										if ((sum_pt_good * sv0_track_which_idx.size() / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv0[i])) >= (sum_pt_poor * sv1_track_which_idx.size() / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv1[i]))) {
											h_more_ratio_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
											h_more_ratio_avg_major_sv_nsv2_fig2->Fill(int(1), w);
										}
										else {
											h_more_ratio_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
											h_more_ratio_avg_major_sv_nsv2_fig2->Fill(int(0), w);

										}

										if ((sum_pt_good / (nsharedjet_tracks_sv0[i])) >= (sum_pt_poor / (nsharedjet_tracks_sv1[i]))) {
											h_more_pT_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
											h_more_pT_avg_major_sv_nsv2_fig2->Fill(int(1), w);
										}
										else {
											h_more_pT_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
											h_more_pT_avg_major_sv_nsv2_fig2->Fill(int(0), w);

										}

										if ((sum_pt_sig_good / (nsharedjet_tracks_sv0[i])) >= (sum_pt_sig_poor / (nsharedjet_tracks_sv1[i]))) {
											h_more_pT_sig_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
											h_more_pT_sig_avg_major_sv_nsv2_fig2->Fill(int(1), w);
										}
										else {
											h_more_pT_sig_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
											h_more_pT_sig_avg_major_sv_nsv2_fig2->Fill(int(0), w);

										}
										
										if (nsharedjet_tracks_sv1[i] == 1) {


											int one_poor_sv1_track_idx = sv1_sharedjet_i_which_trk_idx[0];
											AlgebraicVector3 mom_one_poor(sv1.track_px[one_poor_sv1_track_idx], sv1.track_py[one_poor_sv1_track_idx], sv1.track_pz[one_poor_sv1_track_idx]);
											AlgebraicVector3 ref_one_poor(sv1.track_vx[one_poor_sv1_track_idx], sv1.track_vy[one_poor_sv1_track_idx], sv1.track_vz[one_poor_sv1_track_idx]);

											Measurement1D miss_dist_one_poor_other = miss_dist(sv0, ref_one_poor, mom_one_poor);
											Measurement1D miss_dist_one_poor_its = miss_dist(sv1, ref_one_poor, mom_one_poor);

											h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_fig2->Fill(miss_dist_one_poor_its.significance(), miss_dist_one_poor_other.significance());

											double sum_good_miss_dist_sig_other = 0;
											double sum_good_miss_dist_sig_its = 0;

											for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {

												int more_good_sv0_track_idx = sv0_sharedjet_i_which_trk_idx[j];
												AlgebraicVector3 mom_more_good(sv0.track_px[more_good_sv0_track_idx], sv0.track_py[more_good_sv0_track_idx], sv0.track_pz[more_good_sv0_track_idx]);
												AlgebraicVector3 ref_more_good(sv0.track_vx[more_good_sv0_track_idx], sv0.track_vy[more_good_sv0_track_idx], sv0.track_vz[more_good_sv0_track_idx]);

												Measurement1D miss_dist_more_good_other = miss_dist(sv1, ref_more_good, mom_more_good);
												Measurement1D miss_dist_more_good_its = miss_dist(sv0, ref_more_good, mom_more_good);

												sum_good_miss_dist_sig_other = sum_good_miss_dist_sig_other + miss_dist_more_good_other.significance();
												sum_good_miss_dist_sig_its = sum_good_miss_dist_sig_its + miss_dist_more_good_its.significance();



											}


											h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_fig2->Fill(miss_dist_one_poor_its.significance(), sum_good_miss_dist_sig_its / nsharedjet_tracks_sv0[i]);
											h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_fig2->Fill(miss_dist_one_poor_other.significance(), sum_good_miss_dist_sig_other / nsharedjet_tracks_sv0[i]);

										}


													
										
										h_vertex_chi2dof_absdeltaphi0_large_nsv2_fig2->Fill(sv1.chi2dof(), w);
									}
									else
									{
									h_less_dphi_more_shared_tracks_large_nsv2_all_fig2->Fill(int(0), w);		   // sv0 next to a shared jet and sv1 is away
									h_pt_shared_jets_large_nsv2_semi_fig2->Fill(mevent->jet_pt[jet_index],w);
									h_2D_major_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2->Fill(avg_sv1_non_shared_jet_pt, mevent->jet_pt[jet_index]);
									h_2D_minor_sv_nonshj_avg_pt_shj_pt_large_sv_semi_fig2->Fill(avg_sv0_non_shared_jet_pt, mevent->jet_pt[jet_index]);
									h_dPhi_sv_large_nsv2_semi_fig2->Fill(double(fabs(reco::deltaPhi(phi0, phi1))), w);
									h_dEta_sv_large_nsv2_semi_fig2->Fill(double(fabs(reco::deltaPhi(eta0, eta1))), w);
									h_2D_dEta_dPhi_sv_large_nsv2_semi_fig2->Fill(double(fabs(reco::deltaPhi(eta0, eta1))), double(fabs(reco::deltaPhi(phi0, phi1))));

									

									double ratio_ntracks_nsv2_semi_fig2 = nsharedjet_tracks_sv1[i] / nsharedjet_tracks_sv0[i];
									h_ratio_ntracks_large_nsv2_shared_jets_semi_fig2->Fill(ratio_ntracks_nsv2_semi_fig2);

									std::vector<double> absdeltaphi_max_jet_shared_tracks;
									std::vector<int> sv1_sharedjet_i_which_trk_idx = sv1_sharedjet_which_idx[i];
									std::vector<int> sv1_sharedjet_i_which_no_trk_idx = sv1_sharedjet_which_no_trk_idx[i];
									std::vector<int> sv0_sharedjet_i_which_trk_idx = sv0_sharedjet_which_idx[i];
									std::vector<int> sv0_sharedjet_i_which_no_trk_idx = sv0_sharedjet_which_no_trk_idx[i];
									double sum_pt_good = 0;
									double sum_pt_sig_good = 0;
									for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {
										int track_idx = sv1_sharedjet_i_which_trk_idx[j];
										double absdelta_sv1_track = double(fabs(reco::deltaPhi(phi1, sv1.track_phi[track_idx]))); //phi1
										double absdelta_sv0_track = double(fabs(reco::deltaPhi(phi0, sv1.track_phi[track_idx]))); //phi0

										double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx])));
										h_absdeltaphi_large_jet_shared_tracks_nsv2_semi_fig2->Fill(absdelta_jet_track, w);
										h_dPhi_good_tracks_its_sv_large_semi_fig2->Fill(absdelta_sv1_track, w);
										h_dPhi_good_tracks_other_sv_large_semi_fig2->Fill(absdelta_sv0_track, w);
										h_pt_good_shared_tracks_large_nsv2_semi_fig2->Fill(sv1.track_pt(track_idx), w);
										h_dxy_err_good_shared_tracks_large_nsv2_semi_fig2->Fill(sv1.track_dxy_err(track_idx), w);
										h_dxy_sig_good_shared_tracks_large_nsv2_semi_fig2->Fill(sv1.track_dxy[track_idx] / sv1.track_dxy_err(track_idx), w);
										absdeltaphi_max_jet_shared_tracks.push_back(absdelta_jet_track);
										sum_pt_good = sum_pt_good + sv1.track_pt(track_idx);
										sum_pt_sig_good = sum_pt_sig_good + (sv1.track_pt(track_idx) / sv1.track_pt_err[track_idx]);
										for (int k = 0; k < nsharedjet_tracks_sv0[i]; k++) {
											int track_idx_poor = sv0_sharedjet_i_which_trk_idx[k];
											h_2D_poor_dxy_sig_good_dxy_sig_large_sv_semi_fig2->Fill(sv0.track_dxy[track_idx_poor] / sv0.track_dxy_err(track_idx_poor), sv1.track_dxy[track_idx] / sv1.track_dxy_err(track_idx));
											h_2D_poor_dxy_err_good_dxy_err_large_sv_semi_fig2->Fill(sv0.track_dxy_err(track_idx_poor), sv1.track_dxy_err(track_idx));
											h_2D_poor_dz_err_good_dz_err_large_sv_semi_fig2->Fill(sv0.track_dz_err(track_idx_poor), sv1.track_dz_err(track_idx));
											h_2D_poor_dz_sig_good_dz_sig_large_sv_semi_fig2->Fill(sv0.track_dz[track_idx_poor] / sv0.track_dz_err(track_idx_poor), sv1.track_dz[track_idx] / sv1.track_dz_err(track_idx));

											
											
											h_2D_poor_pt_good_pt_large_sv_semi_fig2->Fill(sv0.track_pt(track_idx_poor), sv1.track_pt(track_idx));

										}

									}
									
									double max_dphi_jet_track = *std::max_element(absdeltaphi_max_jet_shared_tracks.begin(), absdeltaphi_max_jet_shared_tracks.end());
									int max_dphi_jet_track_idx = std::max_element(absdeltaphi_max_jet_shared_tracks.begin(), absdeltaphi_max_jet_shared_tracks.end()) - absdeltaphi_max_jet_shared_tracks.begin();

									int max_sv1_track_idx = sv1_sharedjet_i_which_trk_idx[max_dphi_jet_track_idx];
									

									double sum_pt_poor = 0;
									double sum_pt_sig_poor = 0;
								
									std::vector<double> absdeltaphi_min_jet_shared_tracks;
									for (int j = 0; j < nsharedjet_tracks_sv0[i]; j++) {
										int track_idx = sv0_sharedjet_i_which_trk_idx[j];
										double absdelta_sv0_track = double(fabs(reco::deltaPhi(phi0, sv0.track_phi[track_idx]))); //phi0
										double absdelta_sv1_track = double(fabs(reco::deltaPhi(phi1, sv0.track_phi[track_idx]))); //phi1
										double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx])));                                                                     h_absdeltaphi_large_jet_shared_tracks_nsv2_semi_fig2->Fill(absdelta_jet_track, w);
										h_pt_poor_shared_tracks_large_nsv2_semi_fig2->Fill(sv0.track_pt(track_idx), w);
										h_dPhi_poor_tracks_its_sv_large_semi_fig2->Fill(absdelta_sv0_track, w);
										h_dPhi_poor_tracks_other_sv_large_semi_fig2->Fill(absdelta_sv1_track, w);
										h_dxy_err_poor_shared_tracks_large_nsv2_semi_fig2->Fill(sv0.track_dxy_err(track_idx), w);
										h_dxy_sig_poor_shared_tracks_large_nsv2_semi_fig2->Fill(sv0.track_dxy[track_idx] / sv0.track_dxy_err(track_idx), w);

										absdeltaphi_min_jet_shared_tracks.push_back(absdelta_jet_track);
										sum_pt_poor = sum_pt_poor + sv0.track_pt(track_idx);
										sum_pt_sig_poor = sum_pt_sig_poor + (sv0.track_pt(track_idx) / sv0.track_pt_err[track_idx]);
									}
									h_2D_poor_avg_pt_good_avg_pt_large_sv_semi_fig2->Fill(sum_pt_poor / nsharedjet_tracks_sv0[i], sum_pt_good / nsharedjet_tracks_sv1[i]);
									h_2D_poor_avg_pt_sig_good_avg_pt_sig_large_sv_semi_fig2->Fill(sum_pt_sig_poor / nsharedjet_tracks_sv0[i], sum_pt_sig_good / nsharedjet_tracks_sv1[i]);
									//h_diff_pT_avg_major_minor_sv_nsv2_all_fig2->Fill((sum_pt_good / nsharedjet_tracks_sv1[i]) - (sum_pt_poor / nsharedjet_tracks_sv0[i]), w);
									//h_diff_pT_sum_major_minor_sv_nsv2_all_fig2->Fill(sum_pt_good- sum_pt_poor, w);

									double sum_pt_non_shared_poor = 0;
									for (unsigned int j = 0; j < sv0_sharedjet_i_which_no_trk_idx.size(); j++) {
										idx = sv0_sharedjet_i_which_no_trk_idx[j];
										sum_pt_non_shared_poor = sum_pt_non_shared_poor + sv0.track_pt(idx);
									}
									double sum_pt_non_shared_good = 0;
									for (unsigned int j = 0; j < sv1_sharedjet_i_which_no_trk_idx.size(); j++) {
										idx = sv1_sharedjet_i_which_no_trk_idx[j];
										sum_pt_non_shared_good = sum_pt_non_shared_good + sv1.track_pt(idx);
									}
									
									h_2D_minor_sv_major_sv_ratio_sum_pt_large_sv_semi_fig2->Fill(sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor), sum_pt_good / (sum_pt_good + sum_pt_non_shared_good));
									h_2D_minor_sv_major_sv_ratio_avg_pt_large_sv_semi_fig2->Fill((sum_pt_poor * sv0_track_which_idx.size()) / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv0[i]), (sum_pt_good * sv1_track_which_idx.size()) / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv1[i]));
									//h_diff_ratio_pT_avg_major_minor_sv_nsv2_all_fig2->Fill(((sum_pt_good* sv1_track_which_idx.size()) / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv1[i])) - ((sum_pt_poor * sv0_track_which_idx.size()) / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv0[i])), w);
									//h_diff_ratio_pT_sum_major_minor_sv_nsv2_all_fig2->Fill((sum_pt_good/(sum_pt_good + sum_pt_non_shared_good) ) - (sum_pt_poor /(sum_pt_poor + sum_pt_non_shared_poor)), w);

									//h_ratio_diff_pT_sum_major_minor_sv_nsv2_all_fig2->Fill(((sum_pt_good)-(sum_pt_poor)) / ((sum_pt_good)+(sum_pt_poor)), w);
									double avg_pt_good = sum_pt_good / nsharedjet_tracks_sv1[i];
									double avg_pt_poor = sum_pt_poor / nsharedjet_tracks_sv0[i];
									//h_ratio_diff_pT_avg_major_minor_sv_nsv2_all_fig2->Fill(((avg_pt_good)-(avg_pt_poor)) / ((avg_pt_good)+(avg_pt_poor)), w);
									
									if (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) < 0 || sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) < 0 || sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) == 1 || sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) == 1) {
										std::cout << "found minor's ratio sum pT = " << sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor) << "and major's ratio sum pT = " << sum_pt_good / (sum_pt_good + sum_pt_non_shared_good) << std::endl;

										std::cout << "minor tracks = " << sv0_track_which_idx.size() << ", shared ones = " << nsharedjet_tracks_sv0[i] << ", not shared ones = " << sv0_sharedjet_i_which_no_trk_idx.size() << std::endl;
										for (unsigned int j = 0; j < sv0_sharedjet_i_which_trk_idx.size(); j++) {
											idx = sv0_sharedjet_i_which_trk_idx[j];
											std::cout << "minor's track pT = " << sv0.track_pt(idx) << std::endl;
										}

										std::cout << "major tracks = " << sv1_track_which_idx.size() << ", shared ones = " << nsharedjet_tracks_sv1[i] << ", not shared ones = " << sv1_sharedjet_i_which_no_trk_idx.size() << std::endl;
										for (unsigned int j = 0; j < sv1_sharedjet_i_which_trk_idx.size(); j++) {
											idx = sv1_sharedjet_i_which_trk_idx[j];
											std::cout << "major's track pT = " << sv1.track_pt(idx) << std::endl;
										}
										std::cout << ", sum pT minor tracks = " << sum_pt_poor + sum_pt_non_shared_poor << "sum pT major tracks = " << sum_pt_good + sum_pt_non_shared_good << std::endl;


									}

									if ((sum_pt_good / (sum_pt_good + sum_pt_non_shared_good)) >= (sum_pt_poor / (sum_pt_poor + sum_pt_non_shared_poor))) {
										h_more_ratio_sum_major_sv_nsv2_all_fig2->Fill(int(1), w);
										h_more_ratio_sum_major_sv_nsv2_semi_fig2->Fill(int(1), w);
									}
									else {
										h_more_ratio_sum_major_sv_nsv2_all_fig2->Fill(int(0), w);
										h_more_ratio_sum_major_sv_nsv2_semi_fig2->Fill(int(0), w);

									}

									if ((sum_pt_good * sv1_track_which_idx.size() / ((sum_pt_good + sum_pt_non_shared_good) * nsharedjet_tracks_sv1[i])) >= (sum_pt_poor * sv0_track_which_idx.size() / ((sum_pt_poor + sum_pt_non_shared_poor) * nsharedjet_tracks_sv0[i]))) {
										h_more_ratio_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
										h_more_ratio_avg_major_sv_nsv2_semi_fig2->Fill(int(1), w);
									}
									else {
										h_more_ratio_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
										h_more_ratio_avg_major_sv_nsv2_semi_fig2->Fill(int(0), w);

									}

									if ((sum_pt_good / (nsharedjet_tracks_sv1[i])) >= (sum_pt_poor / (nsharedjet_tracks_sv0[i]))) {
										h_more_pT_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
										h_more_pT_avg_major_sv_nsv2_semi_fig2->Fill(int(1), w);
									}
									else {
										h_more_pT_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
										h_more_pT_avg_major_sv_nsv2_semi_fig2->Fill(int(0), w);

									}

									if ((sum_pt_sig_good / (nsharedjet_tracks_sv1[i])) >= (sum_pt_sig_poor / (nsharedjet_tracks_sv0[i]))) {
										h_more_pT_sig_avg_major_sv_nsv2_all_fig2->Fill(int(1), w);
										h_more_pT_sig_avg_major_sv_nsv2_semi_fig2->Fill(int(1), w);
									}
									else {
h_more_pT_sig_avg_major_sv_nsv2_all_fig2->Fill(int(0), w);
h_more_pT_sig_avg_major_sv_nsv2_semi_fig2->Fill(int(0), w);

									}

									double min_dphi_jet_track = *std::min_element(absdeltaphi_min_jet_shared_tracks.begin(), absdeltaphi_min_jet_shared_tracks.end());
									int min_dphi_jet_track_idx = std::min_element(absdeltaphi_min_jet_shared_tracks.begin(), absdeltaphi_min_jet_shared_tracks.end()) - absdeltaphi_min_jet_shared_tracks.begin();

									int min_sv0_track_idx = sv0_sharedjet_i_which_trk_idx[min_dphi_jet_track_idx];




									if (nsharedjet_tracks_sv0[i] == 1) {


										int one_poor_sv0_track_idx = sv0_sharedjet_i_which_trk_idx[0];
										AlgebraicVector3 mom_one_poor(sv0.track_px[one_poor_sv0_track_idx], sv0.track_py[one_poor_sv0_track_idx], sv0.track_pz[one_poor_sv0_track_idx]);
										AlgebraicVector3 ref_one_poor(sv0.track_vx[one_poor_sv0_track_idx], sv0.track_vy[one_poor_sv0_track_idx], sv0.track_vz[one_poor_sv0_track_idx]);

										Measurement1D miss_dist_one_poor_other = miss_dist(sv1, ref_one_poor, mom_one_poor);
										Measurement1D miss_dist_one_poor_its = miss_dist(sv0, ref_one_poor, mom_one_poor);

										h_2D_1poor_its_sv_1poor_other_sv_miss_dist_sig_large_sv_semi_fig2->Fill(miss_dist_one_poor_its.significance(), miss_dist_one_poor_other.significance());

										double sum_good_miss_dist_sig_other = 0;
										double sum_good_miss_dist_sig_its = 0;

										for (int j = 0; j < nsharedjet_tracks_sv1[i]; j++) {

											int more_good_sv1_track_idx = sv1_sharedjet_i_which_trk_idx[j];
											AlgebraicVector3 mom_more_good(sv1.track_px[more_good_sv1_track_idx], sv1.track_py[more_good_sv1_track_idx], sv1.track_pz[more_good_sv1_track_idx]);
											AlgebraicVector3 ref_more_good(sv1.track_vx[more_good_sv1_track_idx], sv1.track_vy[more_good_sv1_track_idx], sv1.track_vz[more_good_sv1_track_idx]);

											Measurement1D miss_dist_more_good_other = miss_dist(sv0, ref_more_good, mom_more_good);
											Measurement1D miss_dist_more_good_its = miss_dist(sv1, ref_more_good, mom_more_good);

											sum_good_miss_dist_sig_other = sum_good_miss_dist_sig_other + miss_dist_more_good_other.significance();
											sum_good_miss_dist_sig_its = sum_good_miss_dist_sig_its + miss_dist_more_good_its.significance();



										}


										h_2D_1poor_4good_miss_dist_sig_its_sv_large_sv_semi_fig2->Fill(miss_dist_one_poor_its.significance(), sum_good_miss_dist_sig_its / nsharedjet_tracks_sv1[i]);
										h_2D_1poor_4good_miss_dist_sig_other_sv_large_sv_semi_fig2->Fill(miss_dist_one_poor_other.significance(), sum_good_miss_dist_sig_other / nsharedjet_tracks_sv1[i]);

									}



									h_2D_pt_dxy_err_absdeltaphi0_large_sv_semi_fig2->Fill(sv1.track_pt(max_sv1_track_idx), sv1.track_dxy_err(max_sv1_track_idx));
									h_2D_pt_dxy_err_absdeltaphi1_large_sv_semi_fig2->Fill(sv0.track_pt(min_sv0_track_idx), sv0.track_dxy_err(min_sv0_track_idx));

									h_vertex_chi2dof_absdeltaphi0_large_nsv2_semi_fig2->Fill(sv1.chi2dof(), w);
									}







							}

						//}



					   }	//end nsharedjets loop



					   
					   if ((sv0_sum_pt_track_which_idx.size() >= 5) && (sv1_sum_pt_track_which_idx.size() >= 5)) {
						   
							   h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2->Fill(sv0_sum_pt_track_which_idx.size(), sv1_sum_pt_track_which_idx.size(),w);
							   h_sv_njets_no_less_sum_pt_shared_tracks_large_nsv2->Fill(njets_sum_pT_sv0, w);
							   h_sv_njets_no_less_sum_pt_shared_tracks_large_nsv2->Fill(njets_sum_pT_sv1, w);

					   }
					   else {
							   h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2->Fill(sv0_sum_pt_track_which_idx.size(), sv1_sum_pt_track_which_idx.size(),w);
					   }

					   if (nsharedjets == 1) {

						   if ((sv0_sum_pt_track_which_idx.size() >= 5) && (sv1_sum_pt_track_which_idx.size() >= 5)) {

							   h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj1->Fill(sv0_sum_pt_track_which_idx.size(), sv1_sum_pt_track_which_idx.size(),w);

						   }
						   else {
							   h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj1->Fill(sv0_sum_pt_track_which_idx.size(), sv1_sum_pt_track_which_idx.size(),w);
						   }

					   }

					   if (nsharedjets == 2) {

						   if ((sv0_sum_pt_track_which_idx.size() >= 5) && (sv1_sum_pt_track_which_idx.size() >= 5)) {

							   h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj2->Fill(sv0_sum_pt_track_which_idx.size(), sv1_sum_pt_track_which_idx.size(),w);

						   }
						   else {
							   h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_nshj2->Fill(sv0_sum_pt_track_which_idx.size(), sv1_sum_pt_track_which_idx.size(),w);
						   }

					   }

					   if (nsharedjets >= 3) {

						   if ((sv0_sum_pt_track_which_idx.size() >= 5) && (sv1_sum_pt_track_which_idx.size() >= 5)) {

							   h_2D_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_more_nshj3->Fill(sv0_sum_pt_track_which_idx.size(), sv1_sum_pt_track_which_idx.size(),w);

						   }
						   else {
							   h_2D_poor_sv_tracks_no_less_sum_pt_shared_tracks_large_nsv2_more_nshj3->Fill(sv0_sum_pt_track_which_idx.size(), sv1_sum_pt_track_which_idx.size(),w);
						   }

					   }

					   if ((sv0_ratio_sum_pt_track_which_idx.size() >= 5) && (sv1_ratio_sum_pt_track_which_idx.size() >= 5)) {

						   h_2D_sv_tracks_no_less_ratio_sum_pt_shared_tracks_large_nsv2->Fill(sv0_ratio_sum_pt_track_which_idx.size(), sv1_ratio_sum_pt_track_which_idx.size(), w);
						   h_sv_njets_no_less_ratio_sum_pt_shared_tracks_large_nsv2->Fill(njets_ratio_sum_pT_sv0, w);
						   h_sv_njets_no_less_ratio_sum_pt_shared_tracks_large_nsv2->Fill(njets_ratio_sum_pT_sv1, w);
					   }
					   else {
						   h_2D_poor_sv_tracks_no_less_ratio_sum_pt_shared_tracks_large_nsv2->Fill(sv0_ratio_sum_pt_track_which_idx.size(), sv1_ratio_sum_pt_track_which_idx.size(), w);
					   }

					   if ((sv0_avg_pt_track_which_idx.size() >= 5) && (sv1_avg_pt_track_which_idx.size() >= 5)) {

						   h_2D_sv_tracks_no_less_avg_pt_shared_tracks_large_nsv2->Fill(sv0_avg_pt_track_which_idx.size(), sv1_avg_pt_track_which_idx.size(), w);

					   }
					   else {
						   h_2D_poor_sv_tracks_no_less_avg_pt_shared_tracks_large_nsv2->Fill(sv0_avg_pt_track_which_idx.size(), sv1_avg_pt_track_which_idx.size(), w);
					   }

					   if ((sv0_ratio_avg_pt_track_which_idx.size() >= 5) && (sv1_ratio_avg_pt_track_which_idx.size() >= 5)) {

						   h_2D_sv_tracks_no_less_ratio_avg_pt_shared_tracks_large_nsv2->Fill(sv0_ratio_avg_pt_track_which_idx.size(), sv1_ratio_avg_pt_track_which_idx.size(), w);

					   }
					   else {
						   h_2D_poor_sv_tracks_no_less_ratio_avg_pt_shared_tracks_large_nsv2->Fill(sv0_ratio_avg_pt_track_which_idx.size(), sv1_ratio_avg_pt_track_which_idx.size(), w);
					   }

					   

					   
					   

					}   //end no split vertex 

				}	//end nsv=2






		   }
		else {	   // no shared jets

		




		h_lspdist2d_nsv2_no_shared_jets->Fill(mevent->lspdist2d(), w);
		h_lspdist3d_nsv2_no_shared_jets->Fill(mevent->lspdist3d(), w);
		h_absdeltaphi01_genlsp_nsv2_no_shared_jets->Fill(std::abs(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1])), w);
		if (nsv == 2) {
			std::vector<int> sv0_track_which_idx(int(sv0.ntracks()));
			int idx0 = 0;
			std::generate(sv0_track_which_idx.begin(), sv0_track_which_idx.end(), [&] { return idx0++; });
			std::vector<int> sv1_track_which_idx(int(sv1.ntracks()));
			int idx1 = 0;
			std::generate(sv1_track_which_idx.begin(), sv1_track_which_idx.end(), [&] { return idx1++; });

			double sum_pt_sv0 = 0.0;
			double sum_pt_sv1 = 0.0;
			for (unsigned int j = 0; j < sv0_track_which_idx.size(); j++) {
				int idx = sv0_track_which_idx[j];
				sum_pt_sv0 = sum_pt_sv0 + sv0.track_pt(idx);
			}
			for (unsigned int j = 0; j < sv1_track_which_idx.size(); j++) {
				int idx = sv1_track_which_idx[j];
				sum_pt_sv1 = sum_pt_sv1 + sv1.track_pt(idx);
			}

			//let sv0 be like a major vtx 
			h_ratio_diff_pT_sum_sv_nsv2_no_shj->Fill((sum_pt_sv0 - sum_pt_sv1)/(sum_pt_sv0 + sum_pt_sv1),w);
			
			if (sv0_track_which_idx.size() >= sv1_track_which_idx.size()) {
				h_ratio_diff_pT_sum_major_minor_sv_nsv2_no_shj->Fill((sum_pt_sv0 - sum_pt_sv1) / (sum_pt_sv0 + sum_pt_sv1), w);
			}
			else {
				h_ratio_diff_pT_sum_major_minor_sv_nsv2_no_shj->Fill((sum_pt_sv1 - sum_pt_sv0) / (sum_pt_sv1 + sum_pt_sv0), w);

			}
			
			if (fabs(reco::deltaPhi(phi0, phi1)) > 0.5){
				h_ratio_diff_pT_sum_sv_nsv2_large_no_shj->Fill((sum_pt_sv0 - sum_pt_sv1) / (sum_pt_sv0 + sum_pt_sv1), w);
				if (sv0_track_which_idx.size() >= sv1_track_which_idx.size()) {
					h_ratio_diff_pT_sum_major_minor_sv_nsv2_large_no_shj->Fill((sum_pt_sv0 - sum_pt_sv1) / (sum_pt_sv0 + sum_pt_sv1), w);
				}
				else {
					h_ratio_diff_pT_sum_major_minor_sv_nsv2_large_no_shj->Fill((sum_pt_sv1 - sum_pt_sv0) / (sum_pt_sv1 + sum_pt_sv0), w);

				}

			}	

			   h_nsharedjets_nsv2_shared_jets->Fill((int)0, w);
			   if ((reco::deltaPhi(phi0, phi1)) > 0.5) { h_nsharedjets_large_nsv2_shared_jets->Fill((int)0, w); }
		   }

		  
		}

		
	}
 //  }
	
	
  
}

DEFINE_FWK_MODULE(MFVVertexHistos);

