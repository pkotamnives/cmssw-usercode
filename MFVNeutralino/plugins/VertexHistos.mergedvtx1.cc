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
	
	
	TH1F* h_lspdist2d;
	TH1F* h_lspdist3d;
	TH1F* h_dphi_genlsp;
	TH2F* h_2D_sv0lsp0_sv0lsp1;
	TH2F* h_2D_sv1lsp0_sv1lsp1;
	TH2F* h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range0;
	TH2F* h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range1;
	TH2F* h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range2;
	TH2F* h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range0;
	TH2F* h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range1;
	TH2F* h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range2;

	TH2F* h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range0;
	TH2F* h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range1;
	TH2F* h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range2;
	TH2F* h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range3;
	TH2F* h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range0;
	TH2F* h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range1;
	TH2F* h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range2;
	TH2F* h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range3;

	TH1F* h_dphi_genlsp_nsv2;
	TH1F* h_dphi_sv0sv1_nsv2;
	TH2F* h_2D_less_tracks_dphi01_nsv2;
	TH2F* h_2D_more_tracks_less_tracks_small_nsv2;
	TH1F* h_svdist2d_split_sv_pair_nsv2;
	TH1F* h_svdist3d_split_sv_pair_nsv2;
	TH1F* h_svdistz_split_sv_pair_nsv2;
	


	TH1F* h_dphi_genlsp_nsv3;
	TH1F* h_dphi_sv0sv1_nsv3;
	TH1F* h_dphi_sv0sv2_nsv3;
	TH1F* h_dphi_sv1sv2_nsv3;
	TH1F* h_min_dphi_nsv3;
	TH1F* h_max_dphi_nsv3;
	TH1F* h_min_svdist2d_nsv3;
	TH1F* h_max_svdist2d_nsv3;
	TH1F* h_svdist2d_split_sv_pair_nsv3;
	TH1F* h_svdist3d_split_sv_pair_nsv3;
	TH1F* h_svdistz_split_sv_pair_nsv3;

	TH2F* h_2D_max_dphi_min_dphi_nsv3;
	TH2F* h_2D_max_svdist2d_min_svdist2d_nsv3;
	TH2F* h_2D_min_dphi_its_svdist2d_nsv3;
	TH2F* h_2D_less_tracks_min_dphi_nsv3;
	TH2F* h_2D_less_tracks_min_svdist2d_nsv3;
	TH2F* h_2D_less_tracks_dphi_nsv3;
	
	TH2F* h_2D_more_tracks_less_tracks_small_sv0sv1_nsv3;
	TH2F* h_2D_more_tracks_less_tracks_small_sv0sv2_nsv3;
	TH2F* h_2D_more_tracks_less_tracks_small_sv1sv2_nsv3;

	TH2F* h_2D_less_tracks_dphi_nsv4;
	TH1F* h_svdist2d_split_sv_pair_nsv4;
	TH1F* h_svdist3d_split_sv_pair_nsv4;
	TH1F* h_svdistz_split_sv_pair_nsv4;

	TH2F* h_2D_less_tracks_dphi_nsv5;
	
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

	h_lspdist2d = fs->make<TH1F>("h_lspdist2d", "nsv >= 2;dist2d(gen vtx #0, #1) (cm)", 200, 0, 2);
	h_lspdist3d = fs->make<TH1F>("h_lspdist3d", " nsv >= 2;dist3d(gen vtx #0, #1) (cm)", 200, 0, 2);
	h_dphi_genlsp = fs->make<TH1F>("h_dphi_genlsp", "nsv >= 2;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);
	h_2D_sv0lsp0_sv0lsp1 = fs->make<TH2F>("h_2D_sv0lsp0_sv0lsp1", "nsv >= 2; dist2d(gen vtx #0, sv #0) (cm); dist2d(gen vtx #1, sv #0) (cm)", 1000, 0, 1.0, 1000, 0, 1.0);
	h_2D_sv1lsp0_sv1lsp1 = fs->make<TH2F>("h_2D_sv1lsp0_sv1lsp1", "nsv >= 2; dist2d(gen vtx #0, sv #1) (cm); dist2d(gen vtx #1, sv #1) (cm)", 1000, 0, 1.0, 1000, 0, 1.0);
	h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range0 = fs->make<TH2F>("h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range0", "nsv >= 2, split-vertex pair w/ svdist2d <= 0.01 cm; dist2d(gen vtx #0, sv w/ less tracks) (cm); dist2d(gen vtx #1, sv w/ less tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range1 = fs->make<TH2F>("h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range1", "nsv >= 2, split-vertex pair w/ 0.01 < svdist2d < 0.02 cm; dist2d(gen vtx #0, sv w/ less tracks) (cm); dist2d(gen vtx #1, sv w/ less tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range2 = fs->make<TH2F>("h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range2", "nsv >= 2, split-vertex pair w/  0.02 cm <= svdist2d; dist2d(gen vtx #0, sv w/ less tracks) (cm); dist2d(gen vtx #1, sv w/ less tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range0 = fs->make<TH2F>("h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range0", "nsv >= 2, split-vertex pair w/ svdist2d <= 0.01 cm; dist2d(gen vtx #0, sv w/ more tracks) (cm); dist2d(gen vtx #1, sv w/ more tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range1 = fs->make<TH2F>("h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range1", "nsv >= 2, split-vertex pair w/ 0.01 < svdist2d < 0.02 cm; dist2d(gen vtx #0, sv w/ more tracks) (cm); dist2d(gen vtx #1, sv w/ more tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range2 = fs->make<TH2F>("h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range2", "nsv >= 2, split-vertex pair w/  0.02 cm <= svdist2d; dist2d(gen vtx #0, sv w/ more tracks) (cm); dist2d(gen vtx #1, sv w/ more tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	
	h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range0 = fs->make<TH2F>("h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range0", "nsv >= 2, split-vertex pair w/ svdist2d <= 0.01 cm; dist3d(gen vtx #0, sv w/ less tracks) (cm); dist3d(gen vtx #1, sv w/ less tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range1 = fs->make<TH2F>("h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range1", "nsv >= 2, split-vertex pair w/ 0.01 < svdist2d < 0.02 cm; dist3d(gen vtx #0, sv w/ less tracks) (cm); dist3d(gen vtx #1, sv w/ less tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range2 = fs->make<TH2F>("h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range2", "nsv >= 2, split-vertex pair w/  0.02 cm <= svdist2d; dist3d(gen vtx #0, sv w/ less tracks) (cm); dist3d(gen vtx #1, sv w/ less tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range3 = fs->make<TH2F>("h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range3", "nsv >= 2, split-vertex pair w/  0.1 cm <= svdist2d; dist3d(gen vtx #0, sv w/ less tracks) (cm); dist3d(gen vtx #1, sv w/ less tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);

	h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range0 = fs->make<TH2F>("h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range0", "nsv >= 2, split-vertex pair w/ svdist2d <= 0.01 cm; dist3d(gen vtx #0, sv w/ more tracks) (cm); dist3d(gen vtx #1, sv w/ more tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range1 = fs->make<TH2F>("h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range1", "nsv >= 2, split-vertex pair w/ 0.01 < svdist2d < 0.02 cm; dist3d(gen vtx #0, sv w/ more tracks) (cm); dist3d(gen vtx #1, sv w/ more tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range2 = fs->make<TH2F>("h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range2", "nsv >= 2, split-vertex pair w/  0.02 cm <= svdist2d; dist3d(gen vtx #0, sv w/ more tracks) (cm); dist3d(gen vtx #1, sv w/ more tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);
	h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range3 = fs->make<TH2F>("h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range3", "nsv >= 2, split-vertex pair w/  0.1 cm <= svdist2d; dist3d(gen vtx #0, sv w/ less tracks) (cm); dist3d(gen vtx #1, sv w/ less tracks) (cm)", 100, 0, 1.0, 100, 0, 1.0);


	//nsv=2
	h_dphi_genlsp_nsv2 = fs->make<TH1F>("h_dphi_genlsp_nsv2", "nsv = 2;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);
	h_dphi_sv0sv1_nsv2 = fs->make<TH1F>("h_dphi_sv0sv1_nsv2", "nsv = 2;abs(delta(phi of sv #0, phi of sv #1));arb. units", 316, 0, 3.16);
	h_2D_less_tracks_dphi01_nsv2 = fs->make<TH2F>("h_2D_less_tracks_dphi01_nsv2", "nsv = 2; # of less tracks(sv1's); abs(delta(phi of sv #0, phi of sv #1))", 50, 0, 50, 316, 0, 3.16);
	h_2D_more_tracks_less_tracks_small_nsv2 = fs->make<TH2F>("h_2D_more_tracks_less_tracks_small_nsv2", "nsv = 2, absdPhi01 <= 0.5; # of more tracks(sv0's); # of less tracks(sv1's)", 50, 0, 50, 50, 0, 50);
	h_svdist2d_split_sv_pair_nsv2 = fs->make<TH1F>("h_svdist2d_split_sv_pair_nsv2", "nsv = 2, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdist2d of sv pair (cm);arb. units", 100, 0, 1);
	h_svdist3d_split_sv_pair_nsv2 = fs->make<TH1F>("h_svdist3d_split_sv_pair_nsv2", "nsv = 2, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdist3d of sv pair (cm);arb. units", 100, 0, 1);
	h_svdistz_split_sv_pair_nsv2 = fs->make<TH1F>("h_svdistz_split_sv_pair_nsv2", "nsv = 2, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdistz of sv pair (cm);arb. units", 100, 0, 1);


	h_dphi_genlsp_nsv3 = fs->make<TH1F>("h_dphi_genlsp_nsv3", "nsv = 3;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);
	h_dphi_sv0sv1_nsv3 = fs->make<TH1F>("h_dphi_sv0sv1_nsv3", "nsv = 3;abs(delta(phi of sv #0, phi of sv #1));arb. units", 316, 0, 3.16);
	h_dphi_sv0sv2_nsv3 = fs->make<TH1F>("h_dphi_sv0sv2_nsv3", "nsv = 3;abs(delta(phi of sv #0, phi of sv #2));arb. units", 316, 0, 3.16);
	h_dphi_sv1sv2_nsv3 = fs->make<TH1F>("h_dphi_sv1sv2_nsv3", "nsv = 3;abs(delta(phi of sv #1, phi of sv #2));arb. units", 316, 0, 3.16);
	h_min_dphi_nsv3 = fs->make<TH1F>("h_min_dphi_nsv3", "nsv = 3;min(abs(dphi));arb. units", 316, 0, 3.16);
	h_max_dphi_nsv3 = fs->make<TH1F>("h_max_dphi_nsv3", "nsv = 3;max(abs(dphi));arb. units", 316, 0, 3.16);
	h_min_svdist2d_nsv3 = fs->make<TH1F>("h_min_svdist2d_nsv3", "nsv = 3;min(svdist2d) (cm);arb. units", 100, 0, 0.1);
	h_max_svdist2d_nsv3 = fs->make<TH1F>("h_max_svdist2d_nsv3", "nsv = 3;max(svdist2d) (cm);arb. units", 1000, 0, 1.0);
	h_svdist2d_split_sv_pair_nsv3 = fs->make<TH1F>("h_svdist2d_split_sv_pair_nsv3", "nsv = 3, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdist2d of sv pair (cm);arb. units", 100, 0, 1);
	h_svdist3d_split_sv_pair_nsv3 = fs->make<TH1F>("h_svdist3d_split_sv_pair_nsv3", "nsv = 3, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdist3d of sv pair (cm);arb. units", 100, 0, 1);
	h_svdistz_split_sv_pair_nsv3 = fs->make<TH1F>("h_svdistz_split_sv_pair_nsv3", "nsv = 3, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdistz of sv pair (cm);arb. units", 100, 0, 1);

	
	h_2D_max_dphi_min_dphi_nsv3 = fs->make<TH2F>("h_2D_max_dphi_min_dphi_nsv3", "nsv = 3; max dphi; min dphi", 316, 0, 3.16, 316, 0, 3.16);
	h_2D_max_svdist2d_min_svdist2d_nsv3 = fs->make<TH2F>("h_2D_max_svdist2d_min_svdist2d_nsv3", "nsv = 3; max svdist2d (cm); min svdist2d (cm)", 1000, 0, 1.0, 100, 0, 0.1);
	
	h_2D_min_dphi_its_svdist2d_nsv3 = fs->make<TH2F>("h_2D_min_dphi_its_svdist2d_nsv3", "nsv = 3; min dphi; its svdist2d (cm)", 316, 0, 3.16, 100, 0, 0.1);
	h_2D_less_tracks_min_dphi_nsv3 = fs->make<TH2F>("h_2D_less_tracks_min_dphi_nsv3", "nsv = 3; # of less tracks(sv's w/ min dphi); min dphi", 50, 0, 50, 316, 0, 3.16);
	h_2D_less_tracks_min_svdist2d_nsv3 = fs->make<TH2F>("h_2D_less_tracks_min_svdist2d_nsv3", "nsv = 3; # of less tracks(sv's w/ min svdist2d); min svdist2d (cm)", 50, 0, 50, 100, 0, 1);
	h_2D_less_tracks_dphi_nsv3 = fs->make<TH2F>("h_2D_less_tracks_dphi_nsv3", "nsv = 3; # of less tracks(sv's); dphi of any sv pairs", 50, 0, 50, 316, 0, 3.16);

	h_2D_more_tracks_less_tracks_small_sv0sv1_nsv3 = fs->make<TH2F>("h_2D_more_tracks_less_tracks_small_sv0sv1_nsv3", "nsv = 3, absdPhi01 <= 0.5; # of more tracks(sv0's); # of less tracks(sv1's)", 50, 0, 50, 50, 0, 50);
	h_2D_more_tracks_less_tracks_small_sv0sv2_nsv3 = fs->make<TH2F>("h_2D_more_tracks_less_tracks_small_sv0sv2_nsv3", "nsv = 3, absdPhi02 <= 0.5; # of more tracks(sv0's); # of less tracks(sv2's)", 50, 0, 50, 50, 0, 50);
	h_2D_more_tracks_less_tracks_small_sv1sv2_nsv3 = fs->make<TH2F>("h_2D_more_tracks_less_tracks_small_sv1sv2_nsv3", "nsv = 3, absdPhi12 <= 0.5; # of more tracks(sv1's); # of less tracks(sv2's)", 50, 0, 50, 50, 0, 50);


	h_2D_less_tracks_dphi_nsv4 = fs->make<TH2F>("h_2D_less_tracks_dphi_nsv4", "nsv = 4; # of less tracks(sv's); dphi of any sv pairs", 50, 0, 50, 316, 0, 3.16);
	h_svdist2d_split_sv_pair_nsv4 = fs->make<TH1F>("h_svdist2d_split_sv_pair_nsv4", "nsv = 4, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdist2d of sv pair (cm);arb. units", 100, 0, 1);
	h_svdist3d_split_sv_pair_nsv4 = fs->make<TH1F>("h_svdist3d_split_sv_pair_nsv4", "nsv = 4, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdist3d of sv pair (cm);arb. units", 100, 0, 1);
	h_svdistz_split_sv_pair_nsv4 = fs->make<TH1F>("h_svdistz_split_sv_pair_nsv4", "nsv = 4, absdPhiSVs <= 0.5, # of less tracks < 10 ;svdistz of sv pair (cm);arb. units", 100, 0, 1);

	h_2D_less_tracks_dphi_nsv5 = fs->make<TH2F>("h_2D_less_tracks_dphi_nsv5", "nsv = 5; # of less tracks(sv's); dphi of any sv pairs", 50, 0, 50, 316, 0, 3.16);
	
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


		double lsp0_z = mevent->gen_lsp_decay[2];
		double lsp0_x = mevent->gen_lsp_decay[0];
		double lsp0_y = mevent->gen_lsp_decay[1];
		

		double lsp1_z = mevent->gen_lsp_decay[5];
		double lsp1_x = mevent->gen_lsp_decay[3];
		double lsp1_y = mevent->gen_lsp_decay[4];
		

		std::vector<double> dphi_vec;
		std::vector<double> sv_x_vec;
		std::vector<double> sv_y_vec;
		std::vector<double> phi_vec;
		std::vector<double> svdist2d_vec;
		std::vector<int> less_tracks_vec;

		h_lspdist3d->Fill(mevent->lspdist3d(), w);
		h_lspdist2d->Fill(mevent->lspdist2d(), w);
		h_dphi_genlsp->Fill(std::abs(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1])), w);

		bool shared_jet = std::find_first_of(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()) != sv_track_which_jet[0].end();

		double sv0lsp0 = mag(sv0.x - lsp0_x, sv0.y - lsp0_y);
		double sv0lsp1 = mag(sv0.x - lsp1_x, sv0.y - lsp1_y);

		double sv1lsp0 = mag(sv1.x - lsp0_x, sv1.y - lsp0_y);
		double sv1lsp1 = mag(sv1.x - lsp1_x, sv1.y - lsp1_y);

		h_2D_sv0lsp0_sv0lsp1->Fill(sv0lsp0, sv0lsp1, w);
		h_2D_sv1lsp0_sv1lsp1->Fill(sv1lsp0, sv1lsp1, w);

		if (nsv == 2) {

			h_dphi_genlsp_nsv2->Fill(std::abs(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1])), w);
			double dphi01 = std::abs(reco::deltaPhi(phi0,phi1));
			h_dphi_sv0sv1_nsv2->Fill(dphi01, w);
			h_2D_less_tracks_dphi01_nsv2->Fill(int(sv1.ntracks()),dphi01, w);
			if (dphi01 <= 0.5) {
				h_2D_more_tracks_less_tracks_small_nsv2->Fill(int(sv0.ntracks()), int(sv1.ntracks()), w);
				if (int(sv1.ntracks()) < 10) {
					h_svdist2d_split_sv_pair_nsv2->Fill(double(mag(sv0.x - sv1.x, sv0.y - sv1.y)),w);
					h_svdist3d_split_sv_pair_nsv2->Fill(double(mag(sv0.x - sv1.x, sv0.y - sv1.y, sv0.z - sv1.z)), w);
					h_svdistz_split_sv_pair_nsv2->Fill(double(mag(sv0.z - sv1.z)), w);
					if (double(mag(sv0.x - sv1.x, sv0.y - sv1.y)) <= 0.01) {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

					}
					else if (0.01 < double(mag(sv0.x - sv1.x, sv0.y - sv1.y)) && double(mag(sv0.x - sv1.x, sv0.y - sv1.y)) < 0.02) {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

					}

					else {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

						if (0.1 < double(mag(sv0.x - sv1.x, sv0.y - sv1.y))) {
							h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);
							h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

						}

					}
				}
			}
		}

		if (nsv == 3) {
			const MFVVertexAux& sv2 = auxes->at(2);
			double phi2 = atan2(sv2.y - bsy, sv2.x - bsx);
			double eta2 = atan2(sv2.y - bsy, sv2.z - bsz);

			double dphi01 = std::abs(reco::deltaPhi(phi0, phi1));
			double dphi02 = std::abs(reco::deltaPhi(phi0, phi2));
			double dphi12 = std::abs(reco::deltaPhi(phi1, phi2));
			dphi_vec.push_back(dphi01);
			dphi_vec.push_back(dphi02);
			dphi_vec.push_back(dphi12);
			double max_dphi_nsv3 = *std::max_element(dphi_vec.begin(), dphi_vec.end());
			double min_dphi_nsv3 = *std::min_element(dphi_vec.begin(), dphi_vec.end());
			int min_dphi_nsv3_idx = std::min_element(dphi_vec.begin(), dphi_vec.end()) - dphi_vec.begin();

			double svdist2d_01 = svdist2d;
			double svdist2d_02 = mag(sv0.x - sv2.x, sv0.y - sv2.y);
			double svdist2d_12 = mag(sv1.x - sv2.x, sv1.y - sv2.y);

			svdist2d_vec.push_back(svdist2d_01);
			svdist2d_vec.push_back(svdist2d_02);
			svdist2d_vec.push_back(svdist2d_12);
			double max_svdist2d_nsv3 = *std::max_element(svdist2d_vec.begin(), svdist2d_vec.end());
			double min_svdist2d_nsv3 = *std::min_element(svdist2d_vec.begin(), svdist2d_vec.end());
			int min_svdist2d_nsv3_idx = std::min_element(svdist2d_vec.begin(), svdist2d_vec.end()) - svdist2d_vec.begin();

			less_tracks_vec.push_back(int(sv1.ntracks()));
			less_tracks_vec.push_back(int(sv2.ntracks()));
			less_tracks_vec.push_back(int(sv2.ntracks()));


			h_dphi_genlsp_nsv3->Fill(std::abs(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1])), w);
			h_dphi_sv0sv1_nsv3->Fill(dphi01, w);
			h_2D_less_tracks_dphi_nsv3->Fill(int(sv1.ntracks()), dphi01, w);
			h_dphi_sv0sv2_nsv3->Fill(dphi02, w);
			h_2D_less_tracks_dphi_nsv3->Fill(int(sv2.ntracks()), dphi02, w);
			h_dphi_sv1sv2_nsv3->Fill(dphi12, w);
			h_2D_less_tracks_dphi_nsv3->Fill(int(sv2.ntracks()), dphi12, w);
			h_min_dphi_nsv3->Fill(min_dphi_nsv3, w);
			h_max_dphi_nsv3->Fill(max_dphi_nsv3, w);
			h_min_svdist2d_nsv3->Fill(min_svdist2d_nsv3, w);
			h_max_svdist2d_nsv3->Fill(max_svdist2d_nsv3, w);

			h_2D_max_dphi_min_dphi_nsv3->Fill(max_dphi_nsv3, min_dphi_nsv3, w);
			h_2D_max_svdist2d_min_svdist2d_nsv3->Fill(max_svdist2d_nsv3, min_svdist2d_nsv3, w);
			h_2D_min_dphi_its_svdist2d_nsv3->Fill(min_dphi_nsv3, svdist2d_vec[min_dphi_nsv3_idx], w);
			
			h_2D_less_tracks_min_dphi_nsv3->Fill(less_tracks_vec[min_dphi_nsv3_idx],min_dphi_nsv3, w);
			h_2D_less_tracks_min_svdist2d_nsv3->Fill(less_tracks_vec[min_svdist2d_nsv3_idx], min_svdist2d_nsv3, w);
			
			if (dphi01 <= 0.5) {
				h_2D_more_tracks_less_tracks_small_sv0sv1_nsv3->Fill(int(sv0.ntracks()), int(sv1.ntracks()), w);
				if (int(sv1.ntracks()) < 10) {
					h_svdist2d_split_sv_pair_nsv3->Fill(double(mag(sv0.x - sv1.x, sv0.y - sv1.y)), w);
					h_svdist3d_split_sv_pair_nsv3->Fill(double(mag(sv0.x - sv1.x, sv0.y - sv1.y, sv0.z - sv1.z)), w);
					h_svdistz_split_sv_pair_nsv3->Fill(double(mag(sv0.z - sv1.z)), w);
					if (double(mag(sv0.x - sv1.x, sv0.y - sv1.y)) <= 0.01) {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

					}
					else if (0.01 < double(mag(sv0.x - sv1.x, sv0.y - sv1.y)) && double(mag(sv0.x - sv1.x, sv0.y - sv1.y)) < 0.02) {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

					}

					else {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

						if (0.1 < double(mag(sv0.x - sv1.x, sv0.y - sv1.y))) {
							h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);
							h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

						}

					}
				}
			}

			if (dphi02 <= 0.5) {
				h_2D_more_tracks_less_tracks_small_sv0sv2_nsv3->Fill(int(sv0.ntracks()), int(sv2.ntracks()), w);
				if (int(sv2.ntracks()) < 10) {
					h_svdist2d_split_sv_pair_nsv3->Fill(double(mag(sv0.x - sv2.x, sv0.y - sv2.y)), w);
					h_svdist3d_split_sv_pair_nsv3->Fill(double(mag(sv0.x - sv2.x, sv0.y - sv2.y, sv0.z - sv2.z)), w);
					h_svdistz_split_sv_pair_nsv3->Fill(double(mag(sv0.z - sv2.z)), w);
					if (double(mag(sv0.x - sv2.x, sv0.y - sv2.y)) <= 0.01) {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y, lsp0_z - sv2.z)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y, lsp1_z - sv2.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

					}
					else if (0.01 < double(mag(sv0.x - sv2.x, sv0.y - sv2.y)) && double(mag(sv0.x - sv2.x, sv0.y - sv2.y)) < 0.02) {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y, lsp0_z - sv2.z)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y, lsp1_z - sv2.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

					}

					else {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y, lsp0_z - sv2.z)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y, lsp1_z - sv2.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

						if (0.1 < double(mag(sv0.x - sv2.x, sv0.y - sv2.y))) {
							h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y, lsp0_z - sv2.z)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y, lsp1_z - sv2.z)), w);
							h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv0.x, lsp0_y - sv0.y, lsp0_z - sv0.z)), double(mag(lsp1_x - sv0.x, lsp1_y - sv0.y, lsp1_z - sv0.z)), w);

						}

					}
				}
			}

			if (dphi12 <= 0.5) {
				h_2D_more_tracks_less_tracks_small_sv1sv2_nsv3->Fill(int(sv1.ntracks()), int(sv2.ntracks()), w);
				if (int(sv2.ntracks()) < 10) {
					h_svdist2d_split_sv_pair_nsv3->Fill(double(mag(sv1.x - sv2.x, sv1.y - sv2.y)), w);
					h_svdist3d_split_sv_pair_nsv3->Fill(double(mag(sv1.x - sv2.x, sv1.y - sv2.y, sv1.z - sv2.z)), w);
					h_svdistz_split_sv_pair_nsv3->Fill(double(mag(sv1.z - sv2.z)), w);
					if (double(mag(sv1.x - sv2.x, sv1.y - sv2.y)) <= 0.01) {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y, lsp0_z - sv2.z)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y, lsp1_z - sv2.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);

					}
					else if (0.01 < double(mag(sv1.x - sv2.x, sv1.y - sv2.y)) && double(mag(sv1.x - sv2.x, sv1.y - sv2.y)) < 0.02) {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y, lsp0_z - sv2.z)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y, lsp1_z - sv2.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);

					}

					else {
						h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y)), w);
						h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y)), w);
						h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y, lsp0_z - sv2.z)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y, lsp1_z - sv2.z)), w);
						h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);

						if (0.1 < double(mag(sv1.x - sv2.x, sv1.y - sv2.y))) {
							h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv2.x, lsp0_y - sv2.y, lsp0_z - sv2.z)), double(mag(lsp1_x - sv2.x, lsp1_y - sv2.y, lsp1_z - sv2.z)), w);
							h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv1.x, lsp0_y - sv1.y, lsp0_z - sv1.z)), double(mag(lsp1_x - sv1.x, lsp1_y - sv1.y, lsp1_z - sv1.z)), w);

						}

					}
				}
			}




		}

		if (nsv == 4) {
			const MFVVertexAux& sv2 = auxes->at(2);
			double phi2 = atan2(sv2.y - bsy, sv2.x - bsx);
			const MFVVertexAux& sv3 = auxes->at(3);
			double phi3 = atan2(sv3.y - bsy, sv3.x - bsx);
			phi_vec.push_back(phi0);
			phi_vec.push_back(phi1);
			phi_vec.push_back(phi2);
			phi_vec.push_back(phi3);
			less_tracks_vec.push_back(int(sv1.ntracks()));
			less_tracks_vec.push_back(int(sv2.ntracks()));
			less_tracks_vec.push_back(int(sv3.ntracks()));
			sv_x_vec.push_back(double(sv0.x));
			sv_x_vec.push_back(double(sv1.x));
			sv_x_vec.push_back(double(sv2.x));
			sv_x_vec.push_back(double(sv3.x));
			sv_y_vec.push_back(double(sv0.y));
			sv_y_vec.push_back(double(sv1.y));
			sv_y_vec.push_back(double(sv2.y));
			sv_y_vec.push_back(double(sv3.y));

			
			for (int i = 0; i < nsv; ++i) {
				for (int j = 0; j < nsv; ++j) {
					if (i < j) {
						double dphi = std::abs(reco::deltaPhi(phi_vec[i], phi_vec[j]));
						h_2D_less_tracks_dphi_nsv4->Fill(less_tracks_vec[j], dphi,w);
						if (less_tracks_vec[j] < 10) {
							h_svdist2d_split_sv_pair_nsv4->Fill(double(mag(sv_x_vec[i] - sv_x_vec[j], sv_y_vec[i] - sv_y_vec[j])), w);
							h_svdist3d_split_sv_pair_nsv4->Fill(double(mag(sv_x_vec[i] - sv_x_vec[j], sv_y_vec[i] - sv_y_vec[j], sv_z_vec[i] - sv_z_vec[j])), w);
							h_svdistz_split_sv_pair_nsv4->Fill(double(mag(sv_z_vec[i] - sv_z_vec[j])), w);
							if (double(mag(sv_x_vec[i] - sv_x_vec[j], sv_y_vec[i] - sv_y_vec[j])) <= 0.01) {
								h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv_x_vec[j], lsp0_y - sv_y_vec[j])), double(mag(lsp1_x - sv_x_vec[j], lsp1_y - sv_y_vec[j])), w);
								h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv_x_vec[i], lsp0_y - sv_y_vec[i])), double(mag(lsp1_x - sv_x_vec[i], lsp1_y - sv_y_vec[i])), w);
								h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv_x_vec[j], lsp0_y - sv_y_vec[j], lsp0_z - sv_z_vec[j])), double(mag(lsp1_x - sv_x_vec[j], lsp1_y - sv_y_vec[j], lsp1_z - sv_z_vec[j])), w);
								h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range0->Fill(double(mag(lsp0_x - sv_x_vec[i], lsp0_y - sv_y_vec[i], lsp0_z - sv_z_vec[i])), double(mag(lsp1_x - sv_x_vec[i], lsp1_y - sv_y_vec[i], lsp1_z - sv_z_vec[i])), w);

							}
							else if (0.01 < double(mag(sv_x_vec[i] - sv_x_vec[j], sv_y_vec[i] - sv_y_vec[j])) && double(mag(sv_x_vec[i] - sv_x_vec[j], sv_y_vec[i] - sv_y_vec[j])) < 0.02) {
								h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv_x_vec[j], lsp0_y - sv_y_vec[j])), double(mag(lsp1_x - sv_x_vec[j], lsp1_y - sv_y_vec[j])), w);
								h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv_x_vec[i], lsp0_y - sv_y_vec[i])), double(mag(lsp1_x - sv_x_vec[i], lsp1_y - sv_y_vec[i])), w);
								h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv_x_vec[j], lsp0_y - sv_y_vec[j], lsp0_z - sv_z_vec[j])), double(mag(lsp1_x - sv_x_vec[j], lsp1_y - sv_y_vec[j], lsp1_z - sv_z_vec[j])), w);
								h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range1->Fill(double(mag(lsp0_x - sv_x_vec[i], lsp0_y - sv_y_vec[i], lsp0_z - sv_z_vec[i])), double(mag(lsp1_x - sv_x_vec[i], lsp1_y - sv_y_vec[i], lsp1_z - sv_z_vec[i])), w);

							}

							else {
								h_2D_less_svlsp0_less_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv_x_vec[j], lsp0_y - sv_y_vec[j])), double(mag(lsp1_x - sv_x_vec[j], lsp1_y - sv_y_vec[j])), w);
								h_2D_more_svlsp0_more_svlsp1_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv_x_vec[i], lsp0_y - sv_y_vec[i])), double(mag(lsp1_x - sv_x_vec[i], lsp1_y - sv_y_vec[i])), w);
								h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv_x_vec[j], lsp0_y - sv_y_vec[j], lsp0_z - sv_z_vec[j])), double(mag(lsp1_x - sv_x_vec[j], lsp1_y - sv_y_vec[j], lsp1_z - sv_z_vec[j])), w);
								h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range2->Fill(double(mag(lsp0_x - sv_x_vec[i], lsp0_y - sv_y_vec[i], lsp0_z - sv_z_vec[i])), double(mag(lsp1_x - sv_x_vec[i], lsp1_y - sv_y_vec[i], lsp1_z - sv_z_vec[i])), w);

								if (0.1 < double(mag(sv_x_vec[i] - sv_x_vec[j], sv_y_vec[i] - sv_y_vec[j]))) {
									h_2D_less_svlsp0_less_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv_x_vec[j], lsp0_y - sv_y_vec[j], lsp0_z - sv_z_vec[j])), double(mag(lsp1_x - sv_x_vec[j], lsp1_y - sv_y_vec[j], lsp1_z - sv_z_vec[j])), w);
									h_2D_more_svlsp0_more_svlsp1_dist3d_split_sv_pair_range3->Fill(double(mag(lsp0_x - sv_x_vec[i], lsp0_y - sv_y_vec[i], lsp0_z - sv_z_vec[i])), double(mag(lsp1_x - sv_x_vec[i], lsp1_y - sv_y_vec[i], lsp1_z - sv_z_vec[i])), w);

								}

							}
						}
					}
				}
			}

		}

		if (nsv == 5) {
			const MFVVertexAux& sv2 = auxes->at(2);
			double phi2 = atan2(sv2.y - bsy, sv2.x - bsx);
			const MFVVertexAux & sv3 = auxes->at(3);
			double phi3 = atan2(sv3.y - bsy, sv3.x - bsx);
			const MFVVertexAux& sv4 = auxes->at(4);
			double phi4 = atan2(sv4.y - bsy, sv4.x - bsx);
			phi_vec.push_back(phi0);
			phi_vec.push_back(phi1);
			phi_vec.push_back(phi2);
			phi_vec.push_back(phi3);
			phi_vec.push_back(phi4);
			less_tracks_vec.push_back(int(sv1.ntracks()));
			less_tracks_vec.push_back(int(sv2.ntracks()));
			less_tracks_vec.push_back(int(sv3.ntracks()));
			less_tracks_vec.push_back(int(sv4.ntracks()));
			


			for (int i = 0; i < nsv; ++i) {
				for (int j = 0; j < nsv; ++j) {
					if (i < j) {
						double dphi = std::abs(reco::deltaPhi(phi_vec[i], phi_vec[j]));
						h_2D_less_tracks_dphi_nsv5->Fill(less_tracks_vec[j], dphi, w);
					}
				}
			}

		}
		
	}
 //  }
	
	
  
}

DEFINE_FWK_MODULE(MFVVertexHistos);

