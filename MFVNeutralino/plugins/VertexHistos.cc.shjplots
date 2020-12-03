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

  void fill(TH1F** hs,          const int, const double val,                    const double weight) const { hs[sv_all]->Fill(val, weight); }
  void fill(TH2F** hs,          const int, const double val, const double val2, const double weight) const { hs[sv_all]->Fill(val, val2, weight); }

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

  Measurement1D miss_dist_2D(const reco::Vertex& sv, const AlgebraicVector3& ref, const AlgebraicVector3& mom) {
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

  TH1F* h_w;
  TH1F* h_nsv;
  TH1F* h_svdist2d;
  TH1F* h_svdist3d;
  TH2F* h_sv0pvdz_v_sv1pvdz;
  TH2F* h_sv0pvdzsig_v_sv1pvdzsig;
  TH1F* h_absdeltaphi01;
  TH2F* h_pvmosttracksshared;
  TH1F* h_fractrackssharedwpv01;
  TH1F* h_fractrackssharedwpvs01;
  TH1F* h_sv_shared_jets;
  TH1F* h_svdist2d_shared_jets;
  TH1F* h_svdist2d_no_shared_jets;
  TH1F* h_absdeltaphi01_shared_jets;
  TH1F* h_absdeltaphi01_no_shared_jets;
  TH1F* h_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_absdeltaphi01_nsv2_no_shared_jets;
  TH1F* h_max_absdeltaphi0_sv_nshj1_shared_jets;
  TH1F* h_max_absdeltaphi1_small_sv_nshj1_shared_jets;
  TH1F* h_max_absdeltaphi0_small_sv_nshj1_shared_jets;
  TH1F* h_max_absdeltaphi0_small_sv_shared_jets;
  TH1F* h_max_absdeltaphi1_small_sv_shared_jets;
  TH1F* h_max_absdeltaphi1_sv_nshj1_shared_jets;
  TH1F* h_max_absdeltaphi0_sv_shared_jets;
  TH1F* h_max_absdeltaphi1_sv_shared_jets;
  TH1F* h_max_absdeltaphi0_large_sv_nshj1_shared_tracks;
  TH1F* h_max_absdeltaphi1_large_sv_nshj1_shared_tracks;
  TH1F* h_max_pt_absdeltaphi0_large_sv_nshj1_shared_tracks;
  TH1F* h_max_pt_absdeltaphi1_large_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_absdeltaphi0_large_other_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_absdeltaphi1_large_other_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_absdeltaphi0_large_its_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_absdeltaphi1_large_its_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_err_absdeltaphi0_large_other_sv_nshj1_shared_tracks;                                                                                                              TH1F* h_miss_dist_err_absdeltaphi1_large_other_sv_nshj1_shared_tracks;                                                                                                              TH1F* h_miss_dist_err_absdeltaphi0_large_its_sv_nshj1_shared_tracks;                                                                                                                TH1F* h_miss_dist_err_absdeltaphi1_large_its_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_2D_absdeltaphi0_large_other_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_2D_absdeltaphi1_large_other_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_2D_absdeltaphi0_large_its_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_2D_absdeltaphi1_large_its_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_sig_absdeltaphi0_large_other_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_sig_absdeltaphi1_large_other_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_sig_absdeltaphi0_large_its_sv_nshj1_shared_tracks;
  TH1F* h_miss_dist_sig_absdeltaphi1_large_its_sv_nshj1_shared_tracks;
  TH1F* h_vertex_chi2dof_absdeltaphi0_large_sv2_nshj1;

  TH1F* h_lspdist2d_nsv2_shared_jets;
  TH1F* h_lspdist3d_nsv2_shared_jets;
  TH1F* h_absdeltaphi01_genlsp_nsv2_shared_jets;
  TH1F* h_lspdist2d_nsv2_no_shared_jets;
  TH1F* h_lspdist3d_nsv2_no_shared_jets;
  TH1F* h_absdeltaphi01_genlsp_nsv2_no_shared_jets;

  TH1F* h_nsharedjets_shared_jets;
  TH1F* h_ratio_ntracks_shared_jets;
  TH2F* h_ntracks_sv0sv1_nsharedjet1;
  TH2F* h_ntracks_sv0sv1_nsharedjet2;
  TH2F* h_ntracks_sv0sv1_nsharedjet3;
  TH2F* h_ntracks_sv0sv1_nsharedjet4;
  TH2F* h_ntracks_sv0sv1_nsharedjet5;
  TH1F* h_ntracks_sv0_ntrack1_sv1;
  TH1F* h_ntracks_sv0_ntrack2_sv1;
  TH1F* h_ntracks_sv0_ntrack3_sv1;
  TH1F* h_ntracks_sv0_ntrack4_sv1;
  TH1F* h_ntracks_sv0_ntrack5_sv1;
  TH1F* h_ntracks_nshj1_sv0_ntrack1_sv1;
  TH1F* h_ntracks_nshj1_sv0_ntrack2_sv1;
  TH1F* h_ntracks_nshj1_sv0_ntrack3_sv1;
  TH1F* h_ntracks_nshj1_sv0_ntrack4_sv1;
  TH1F* h_ntracks_nshj1_sv0_ntrack5_sv1;

  TH1F* h_nsharedjets_nsv2_shared_jets;
  TH1F* h_ratio_ntracks_nsv2_shared_jets;
  TH2F* h_ntracks_sv0sv1_nsv2_nsharedjet1;
  TH2F* h_ntracks_sv0sv1_nsv2_nsharedjet2;
  TH2F* h_ntracks_sv0sv1_nsv2_nsharedjet3;
  TH2F* h_ntracks_sv0sv1_nsv2_nsharedjet4;
  TH2F* h_ntracks_sv0sv1_nsv2_nsharedjet5;
  TH1F* h_ntracks_sv0_ntrack1_sv1_nsv2;
  TH1F* h_ntracks_sv0_ntrack2_sv1_nsv2;
  TH1F* h_ntracks_sv0_ntrack3_sv1_nsv2;
  TH1F* h_ntracks_sv0_ntrack4_sv1_nsv2;
  TH1F* h_ntracks_sv0_ntrack5_sv1_nsv2;
  TH1F* h_ntracks_nshj1_sv0_ntrack1_sv1_nsv2;
  TH1F* h_ntracks_nshj1_sv0_ntrack2_sv1_nsv2;
  TH1F* h_ntracks_nshj1_sv0_ntrack3_sv1_nsv2;
  TH1F* h_ntracks_nshj1_sv0_ntrack4_sv1_nsv2;
  TH1F* h_ntracks_nshj1_sv0_ntrack5_sv1_nsv2;

  TH1F* h_nsharedjets_small_nsv2_shared_jets;
  TH1F* h_ratio_ntracks_small_nsv2_nshj1_shared_jets;
  TH2F* h_ntracks_sv0sv1_small_nsv2_nsharedjet1;
  TH2F* h_ntracks_sv0sv1_small_nsv2_nsharedjet2;
  TH2F* h_ntracks_sv0sv1_small_nsv2_nsharedjet3;
  TH2F* h_ntracks_sv0sv1_small_nsv2_nsharedjet4;
  TH2F* h_ntracks_sv0sv1_small_nsv2_nsharedjet5;
  TH1F* h_ntracks_small_sv0_ntrack1_sv1_nsv2;
  TH1F* h_ntracks_small_sv0_ntrack2_sv1_nsv2;
  TH1F* h_ntracks_small_sv0_ntrack3_sv1_nsv2;
  TH1F* h_ntracks_small_sv0_ntrack4_sv1_nsv2;
  TH1F* h_ntracks_small_sv0_ntrack5_sv1_nsv2;
  TH1F* h_ntracks_small_nshj1_sv0_ntrack1_sv1_nsv2;
  TH1F* h_ntracks_small_nshj1_sv0_ntrack2_sv1_nsv2;
  TH1F* h_ntracks_small_nshj1_sv0_ntrack3_sv1_nsv2;
  TH1F* h_ntracks_small_nshj1_sv0_ntrack4_sv1_nsv2;
  TH1F* h_ntracks_small_nshj1_sv0_ntrack5_sv1_nsv2;

  TH1F* h_nsharedjets_large_nsv2_shared_jets;
  TH1F* h_ratio_ntracks_large_nsv2_nshj1_shared_jets;
  TH1F* h_max_absdeltaphi1_large_sv_nshj1_shared_jets;                                                                                                                                TH1F* h_max_absdeltaphi0_large_sv_nshj1_shared_jets;
  TH1F* h_absdeltaphi_large_jet_shared_tracks_nshj1_nsv2;
  TH1F* h_svdist2d_small_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_svdist3d_small_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_svdist2d_large_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_svdist3d_large_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_svdist2d_both_absdeltaphi01_nsv2_no_shared_jets;
  TH1F* h_svdist3d_both_absdeltaphi01_nsv2_no_shared_jets;
 };

const char* MFVVertexHistos::sv_index_names[MFVVertexHistos::sv_num_indices] = { "all" };

MFVVertexHistos::MFVVertexHistos(const edm::ParameterSet& cfg)
  : mevent_token(consumes<MFVEvent>(cfg.getParameter<edm::InputTag>("mevent_src"))),
    weight_token(consumes<double>(cfg.getParameter<edm::InputTag>("weight_src"))),
    vertex_token(consumes<MFVVertexAuxCollection>(cfg.getParameter<edm::InputTag>("vertex_src"))),
    max_ntrackplots(cfg.getParameter<int>("max_ntrackplots")),
    do_scatterplots(cfg.getParameter<bool>("do_scatterplots"))
{
  edm::Service<TFileService> fs;

  h_w = fs->make<TH1F>("h_w", ";event weight;events/0.1", 100, 0, 10);
  h_nsv = fs->make<TH1F>("h_nsv", ";# of secondary vertices;arb. units", 15, 0, 15);
 
  h_svdist2d = fs->make<TH1F>("h_svdist2d", ";dist2d(sv #0, #1) (cm);arb. units", 500, 0, 1);
  h_svdist3d = fs->make<TH1F>("h_svdist3d", ";dist3d(sv #0, #1) (cm);arb. units", 500, 0, 1);
  h_sv0pvdz_v_sv1pvdz = fs->make<TH2F>("h_sv0pvdz_v_sv1pvdz", ";sv #1 dz to PV (cm);sv #0 dz to PV (cm)", 100, 0, 0.5, 100, 0, 0.5);
  h_sv0pvdzsig_v_sv1pvdzsig = fs->make<TH2F>("h_sv0pvdzsig_v_sv1pvdzsig", ";N#sigma(sv #1 dz to PV);sv N#sigma(#0 dz to PV)", 100, 0, 50, 100, 0, 50);
  h_absdeltaphi01 = fs->make<TH1F>("h_absdeltaphi01", ";abs(delta(phi of sv #0, phi of sv #1));arb. units", 315, 0, 3.15);
  h_fractrackssharedwpv01 = fs->make<TH1F>("h_fractrackssharedwpv01", ";fraction of sv #0 and sv #1 tracks shared with the PV;arb. units", 41, 0, 1.025);
  h_fractrackssharedwpvs01 = fs->make<TH1F>("h_fractrackssharedwpvs01", ";fraction of sv #0 and sv #1 tracks shared with any PV;arb. units", 41, 0, 1.025);
  h_pvmosttracksshared = fs->make<TH2F>("h_pvmosttracksshared", ";index of pv most-shared to sv #0; index of pv most-shared to sv #1", 71, -1, 70, 71, -1, 70);
  h_sv_shared_jets  = fs->make<TH1F>("h_sv_shared_jets", ";SV tracks share jet?", 2, 0, 2);
  h_svdist2d_shared_jets = fs->make<TH1F>("h_svdist2d_shared_jets", ";dist2d(sv #0, #1) (cm);arb. units", 500, 0, 1);
  h_svdist2d_no_shared_jets = fs->make<TH1F>("h_svdist2d_no_shared_jets", ";dist2d(sv #0, #1) (cm);arb. units", 500, 0, 1);
  h_absdeltaphi01_shared_jets = fs->make<TH1F>("h_absdeltaphi01_shared_jets", "nsv >= 2;abs(delta(phi of sv #0, phi of sv #1));arb. units", 316, 0, 3.16);
  h_absdeltaphi01_no_shared_jets = fs->make<TH1F>("h_absdeltaphi01_no_shared_jets", "nsv >= 2;abs(delta(phi of sv #0, phi of sv #1));arb. units", 316, 0, 3.16);
  h_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_absdeltaphi01_nsv2_shared_jets", "nsv = 2;abs(delta(phi of sv #0, phi of sv #1));arb. units", 316, 0, 3.16);                   h_absdeltaphi01_nsv2_no_shared_jets = fs->make<TH1F>("h_absdeltaphi01_nsv2_no_shared_jets", "nsv = 2;abs(delta(phi of sv #0, phi of sv #1));arb. units", 316, 0, 3.16);
  h_max_absdeltaphi0_small_sv_nshj1_shared_jets = fs->make<TH1F>("h_max_absdeltaphi0_small_sv_nshj1_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5, nsharedjets = 1;max( dphi(sv0, a shared jet), dphi(sv1, a shared jet) );arb. units", 316, 0, 3.16);
  h_max_absdeltaphi1_small_sv_nshj1_shared_jets = fs->make<TH1F>("h_max_absdeltaphi1_small_sv_nshj1_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5, nsharedjets = 1;abs(delta(phi of max-dphi shared jet, phi of another sv ));arb. units", 316, 0, 3.16);
  h_max_absdeltaphi0_small_sv_shared_jets = fs->make<TH1F>("h_max_absdeltaphi0_small_sv_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5;max( dphi(sv0, shared jets), dphi(sv1, shared jets) );arb. units", 316, 0, 3.16);
  h_max_absdeltaphi1_small_sv_shared_jets = fs->make<TH1F>("h_max_absdeltaphi1_small_sv_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5;abs(delta(phi of max-dphi shared jet, phi of another sv ));arb. units", 316, 0, 3.16);
  h_max_absdeltaphi0_sv_nshj1_shared_jets = fs->make<TH1F>("h_max_absdeltaphi0_sv_nshj1_shared_jets", "nsv = 2, nsharedjets = 1;max( dphi(sv0, a shared jet), dphi(sv1, a shared jet) );arb. units", 316, 0, 3.16);
  h_max_absdeltaphi1_sv_nshj1_shared_jets = fs->make<TH1F>("h_max_absdeltaphi1_sv_nshj1_shared_jets", "nsv = 2, nsharedjets = 1;abs(delta(phi of max-dphi shared jet, phi of another sv ));arb. units", 316, 0, 3.16);
  h_max_absdeltaphi0_sv_shared_jets = fs->make<TH1F>("h_max_absdeltaphi0_sv_shared_jets", "nsv = 2;max( dphi(sv0, shared jets), dphi(sv1, shared jets) );arb. units", 316, 0, 3.16);
  h_max_absdeltaphi1_sv_shared_jets = fs->make<TH1F>("h_max_absdeltaphi1_sv_shared_jets", "nsv = 2;abs(delta(phi of max-dphi shared jet, phi of another sv ));arb. units", 316, 0, 3.16);
  h_max_absdeltaphi0_large_sv_nshj1_shared_jets = fs->make<TH1F>("h_max_absdeltaphi0_large_sv_nshj1_shared_jets", "nsv = 2, absdeltaphi01 > 0.5, nsharedjets = 1;max( dphi(sv0, a shared jet), dphi(sv1, a shared jet) );arb. units", 316, 0, 3.16);                                                                                                                     h_max_absdeltaphi1_large_sv_nshj1_shared_jets = fs->make<TH1F>("h_max_absdeltaphi1_large_sv_nshj1_shared_jets", "nsv = 2, absdeltaphi01 > 0.5, nsharedjets = 1; min( dphi(sv0, a shared jet), dphi(sv1, a shared jet) );arb. units", 316, 0, 3.16);                                                                                                                
  h_max_absdeltaphi0_large_sv_nshj1_shared_tracks = fs->make<TH1F>("h_max_absdeltaphi0_large_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;delta(phi of max track, phi of its vertex);arb. units", 316, 0, 3.16);
  h_max_absdeltaphi1_large_sv_nshj1_shared_tracks = fs->make<TH1F>("h_max_absdeltaphi1_large_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;delta(phi of min track, phi of its vertex);arb. units", 316, 0, 3.16);
  h_max_pt_absdeltaphi0_large_sv_nshj1_shared_tracks = fs->make<TH1F>("h_max_pt_absdeltaphi0_large_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5; max track p_{T} (GeV);arb. units", 200, 0, 200);
  h_max_pt_absdeltaphi1_large_sv_nshj1_shared_tracks = fs->make<TH1F>("h_max_pt_absdeltaphi1_large_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5; min track p_{T} (GeV);arb. units", 200, 0, 200);
  h_miss_dist_absdeltaphi0_large_other_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_absdeltaphi0_large_other_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance(max track, another vertex) (cm);arb. units", 100, 0, 0.5);
  h_miss_dist_absdeltaphi1_large_other_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_absdeltaphi1_large_other_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance(min track, another vertex) (cm);arb. units", 100, 0, 0.5);
  h_miss_dist_absdeltaphi0_large_its_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_absdeltaphi0_large_its_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance(max track, its vertex) (cm);arb. units", 100, 0, 0.5);
  h_miss_dist_absdeltaphi1_large_its_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_absdeltaphi1_large_its_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance(min track, its vertex) (cm);arb. units", 100, 0, 0.5);
  h_miss_dist_err_absdeltaphi0_large_other_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_err_absdeltaphi0_large_other_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance err(max track, another vertex) (cm);arb. units", 100, 0, 0.1);                                                                                         h_miss_dist_err_absdeltaphi1_large_other_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_err_absdeltaphi1_large_other_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance err(min track, another vertex) (cm);arb. units", 100, 0, 0.1);                                                                                         h_miss_dist_err_absdeltaphi0_large_its_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_err_absdeltaphi0_large_its_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance err(max track, its vertex) (cm);arb. units", 100, 0, 0.1);                                                                                                 h_miss_dist_err_absdeltaphi1_large_its_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_err_absdeltaphi1_large_its_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance err(min track, its vertex) (cm);arb. units", 100, 0, 0.1); 
  h_miss_dist_2D_absdeltaphi0_large_other_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_2D_absdeltaphi0_large_other_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance 2D(max track, another vertex) (cm);arb. units", 100, 0, 0.5);
  h_miss_dist_2D_absdeltaphi1_large_other_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_2D_absdeltaphi1_large_other_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance 2D(min track, another vertex) (cm);arb. units", 100, 0, 0.5);
  h_miss_dist_2D_absdeltaphi0_large_its_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_2D_absdeltaphi0_large_its_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance 2D(max track, its vertex) (cm);arb. units", 100, 0, 0.5);
  h_miss_dist_2D_absdeltaphi1_large_its_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_2D_absdeltaphi1_large_its_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance 2D(min track, its vertex) (cm);arb. units", 100, 0, 0.5);
  h_miss_dist_sig_absdeltaphi0_large_other_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_sig_absdeltaphi0_large_other_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance sig(max track, another vertex);arb. units", 100, 0, 10);
  h_miss_dist_sig_absdeltaphi1_large_other_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_sig_absdeltaphi1_large_other_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance sig(min track, another vertex);arb. units", 100, 0, 10);
  h_miss_dist_sig_absdeltaphi0_large_its_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_sig_absdeltaphi0_large_its_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance sig(max track, its vertex);arb. units", 100, 0, 10);
  h_miss_dist_sig_absdeltaphi1_large_its_sv_nshj1_shared_tracks = fs->make<TH1F>("h_miss_dist_sig_absdeltaphi1_large_its_sv_nshj1_shared_tracks", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5;miss distance sig(min track, its vertex);arb. units", 100, 0, 10);
  h_vertex_chi2dof_absdeltaphi0_large_sv2_nshj1 = fs->make<TH1F>("h_vertex_chi2dof_absdeltaphi0_large_sv2_nshj1", "nsv = 2, nsharedjets = 1, absdeltaphi01 > 0.5; sv opposite to the shared jet's tracks #chi^ { 2 } / dof ;arb. units", 100, 0, 10);

  h_lspdist2d_nsv2_shared_jets = fs->make<TH1F>("h_lspdist2d_nsv2_shared_jets", "nsv = 2;dist2d(gen vtx #0, #1) (cm)", 200, 0, 2);
  h_lspdist3d_nsv2_shared_jets = fs->make<TH1F>("h_lspdist3d_nsv2_shared_jets", " nsv = 2;dist3d(gen vtx #0, #1) (cm)", 200, 0, 2);
  h_absdeltaphi01_genlsp_nsv2_shared_jets = fs->make<TH1F>("h_absdeltaphi01_genlsp_nsv2_shared_jets", "nsv = 2;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);
  h_lspdist2d_nsv2_no_shared_jets = fs->make<TH1F>("h_lspdist2d_nsv2_no_shared_jets", "nsv = 2;dist2d(gen vtx #0, #1) (cm)", 200, 0, 2);
  h_lspdist3d_nsv2_no_shared_jets = fs->make<TH1F>("h_lspdist3d_nsv2_no_shared_jets", " nsv = 2;dist3d(gen vtx #0, #1) (cm)", 200, 0, 2);
  h_absdeltaphi01_genlsp_nsv2_no_shared_jets = fs->make<TH1F>("h_absdeltaphi01_genlsp_nsv2_no_shared_jets", "nsv = 2;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);

  h_nsharedjets_shared_jets = fs->make<TH1F>("h_nsharedjets_shared_jets", "nsv >= 2;# of shared jets;arb. units", 10, 0, 10);
  h_ratio_ntracks_shared_jets = fs->make<TH1F>("h_ratio_ntracks_shared_jets", "nsv >= 2;ratios of shared tracks (>=1);arb. units", 50, 0, 10);
  h_ntracks_sv0sv1_nsharedjet1 = fs->make<TH2F>("h_ntracks_sv0sv1_nsharedjet1", "nsv >= 2, nsharedjets = 1;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_nsharedjet2 = fs->make<TH2F>("h_ntracks_sv0sv1_nsharedjet2", "nsv >= 2, nsharedjets = 2;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_nsharedjet3 = fs->make<TH2F>("h_ntracks_sv0sv1_nsharedjet3", "nsv >= 2, nsharedjets = 3;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_nsharedjet4 = fs->make<TH2F>("h_ntracks_sv0sv1_nsharedjet4", "nsv >= 2, nsharedjets = 4;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_nsharedjet5 = fs->make<TH2F>("h_ntracks_sv0sv1_nsharedjet5", "nsv >= 2, nsharedjets = 5;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0_ntrack1_sv1 = fs->make<TH1F>("h_ntracks_sv0_ntrack1_sv1", "nsv >= 2, sv #1 with ntracks = 1 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);   
  h_ntracks_sv0_ntrack2_sv1 = fs->make<TH1F>("h_ntracks_sv0_ntrack2_sv1", "nsv >= 2, sv #1 with ntracks = 2 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_sv0_ntrack3_sv1 = fs->make<TH1F>("h_ntracks_sv0_ntrack3_sv1", "nsv >= 2, sv #1 with ntracks = 3 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_sv0_ntrack4_sv1 = fs->make<TH1F>("h_ntracks_sv0_ntrack4_sv1", "nsv >= 2, sv #1 with ntracks = 4 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_sv0_ntrack5_sv1 = fs->make<TH1F>("h_ntracks_sv0_ntrack5_sv1", "nsv >= 2, sv #1 with ntracks = 5 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack1_sv1 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack1_sv1", "nsv >= 2, sv #1 with ntracks = 1 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack2_sv1 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack2_sv1", "nsv >= 2, sv #1 with ntracks = 2 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack3_sv1 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack3_sv1", "nsv >= 2, sv #1 with ntracks = 3 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack4_sv1 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack4_sv1", "nsv >= 2, sv #1 with ntracks = 4 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack5_sv1 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack5_sv1", "nsv >= 2, sv #1 with ntracks = 5 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
 
  h_nsharedjets_nsv2_shared_jets = fs->make<TH1F>("h_nsharedjets_nsv2_shared_jets", "nsv = 2;# of shared jets;arb. units", 10, 0, 10);
  h_ratio_ntracks_nsv2_shared_jets = fs->make<TH1F>("h_ratio_ntracks_nsv2_shared_jets", "nsv = 2;ratios of shared tracks (>=1);arb. units", 50, 0, 10);
  h_ntracks_sv0sv1_nsv2_nsharedjet1 = fs->make<TH2F>("h_ntracks_sv0sv1_nsv2_nsharedjet1", "nsv = 2, nsharedjets = 1;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_nsv2_nsharedjet2 = fs->make<TH2F>("h_ntracks_sv0sv1_nsv2_nsharedjet2", "nsv = 2, nsharedjets = 2;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_nsv2_nsharedjet3 = fs->make<TH2F>("h_ntracks_sv0sv1_nsv2_nsharedjet3", "nsv = 2, nsharedjets = 3;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_nsv2_nsharedjet4 = fs->make<TH2F>("h_ntracks_sv0sv1_nsv2_nsharedjet4", "nsv = 2, nsharedjets = 4;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_nsv2_nsharedjet5 = fs->make<TH2F>("h_ntracks_sv0sv1_nsv2_nsharedjet5", "nsv = 2, nsharedjets = 5;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0_ntrack1_sv1_nsv2 = fs->make<TH1F>("h_ntracks_sv0_ntrack1_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 1 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_sv0_ntrack2_sv1_nsv2 = fs->make<TH1F>("h_ntracks_sv0_ntrack2_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 2 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_sv0_ntrack3_sv1_nsv2 = fs->make<TH1F>("h_ntracks_sv0_ntrack3_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 3 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_sv0_ntrack4_sv1_nsv2 = fs->make<TH1F>("h_ntracks_sv0_ntrack4_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 4 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_sv0_ntrack5_sv1_nsv2 = fs->make<TH1F>("h_ntracks_sv0_ntrack5_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 5 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack1_sv1_nsv2 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack1_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 1 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack2_sv1_nsv2 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack2_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 2 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack3_sv1_nsv2 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack3_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 3 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack4_sv1_nsv2 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack4_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 4 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_nshj1_sv0_ntrack5_sv1_nsv2 = fs->make<TH1F>("h_ntracks_nshj1_sv0_ntrack5_sv1_nsv2", "nsv = 2, sv #1 with ntracks = 5 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);

  h_nsharedjets_large_nsv2_shared_jets = fs->make<TH1F>("h_nsharedjets_large_nsv2_shared_jets", "nsv = 2, absdeltaphi01 > 0.5;# of shared jets;arb. units", 10, 0, 10); 
  h_ratio_ntracks_large_nsv2_nshj1_shared_jets = fs->make<TH1F>("h_ratio_ntracks_large_nsv2_nshj1_shared_jets", "nsv = 2, absdeltaphi01 > 0.5, nsharedjets = 1;ratios of shared tracks (>=1);arb. units", 50, 0, 10); 
  h_nsharedjets_small_nsv2_shared_jets = fs->make<TH1F>("h_nsharedjets_small_nsv2_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5;# of shared jets;arb. units", 10, 0, 10);
  h_ratio_ntracks_small_nsv2_nshj1_shared_jets = fs->make<TH1F>("h_ratio_ntracks_small_nsv2_nshj1_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5, nsharedjets = 1;ratios of shared tracks (>=1);arb. units", 50, 0, 10);
  h_ntracks_sv0sv1_small_nsv2_nsharedjet1 = fs->make<TH2F>("h_ntracks_sv0sv1_small_nsv2_nsharedjet1", "nsv = 2, absdeltaphi01 <= 0.5, nsharedjets = 1;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_small_nsv2_nsharedjet2 = fs->make<TH2F>("h_ntracks_sv0sv1_small_nsv2_nsharedjet2", "nsv = 2, absdeltaphi01 <= 0.5, nsharedjets = 2;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_small_nsv2_nsharedjet3 = fs->make<TH2F>("h_ntracks_sv0sv1_small_nsv2_nsharedjet3", "nsv = 2, absdeltaphi01 <= 0.5, nsharedjets = 3;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_small_nsv2_nsharedjet4 = fs->make<TH2F>("h_ntracks_sv0sv1_small_nsv2_nsharedjet4", "nsv = 2, absdeltaphi01 <= 0.5, nsharedjets = 4;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_sv0sv1_small_nsv2_nsharedjet5 = fs->make<TH2F>("h_ntracks_sv0sv1_small_nsv2_nsharedjet5", "nsv = 2, absdeltaphi01 <= 0.5, nsharedjets = 5;sv #0 # of tracks from a shared jet;sv #1 # of tracks from a shared jet", 10, 0, 10, 10, 0, 10);
  h_ntracks_small_sv0_ntrack1_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_sv0_ntrack1_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 1 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_sv0_ntrack2_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_sv0_ntrack2_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 2 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_sv0_ntrack3_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_sv0_ntrack3_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 3 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_sv0_ntrack4_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_sv0_ntrack4_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 4 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_sv0_ntrack5_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_sv0_ntrack5_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 5 from a shared jet ;sv #0 # of tracks from a shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_nshj1_sv0_ntrack1_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_nshj1_sv0_ntrack1_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 1 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_nshj1_sv0_ntrack2_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_nshj1_sv0_ntrack2_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 2 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_nshj1_sv0_ntrack3_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_nshj1_sv0_ntrack3_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 3 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_nshj1_sv0_ntrack4_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_nshj1_sv0_ntrack4_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 4 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_ntracks_small_nshj1_sv0_ntrack5_sv1_nsv2 = fs->make<TH1F>("h_ntracks_small_nshj1_sv0_ntrack5_sv1_nsv2", "nsv = 2, absdeltaphi01 <= 0.5, sv #1 with ntracks = 5 from an only-one shared jet ;sv #0 # of tracks from an only-one shared jet;arb. units", 10, 0, 10);
  h_absdeltaphi_large_jet_shared_tracks_nshj1_nsv2 = fs->make<TH1F>("h_absdeltaphi_large_jet_shared_tracks_nshj1_nsv2", "nsv = 2, sharedjets = 1, absdeltaphi01 > 0.5; abs(delta(shared tracks, phi of the only-one shared jet));arb. units", 316, 0, 3.16); 

  h_svdist2d_small_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist2d_small_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5 ;dist2d(sv #0, #1) (cm);arb. units", 500, 0, 1); 
  h_svdist3d_small_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist3d_small_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5 ;dist3d(sv #0, #1) (cm);arb. units", 500, 0, 1); 
  h_svdist2d_large_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist2d_large_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 > 0.5 ;dist2d(sv #0, #1) (cm);arb. units", 500, 0, 1);
  h_svdist3d_large_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist3d_large_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 > 0.5 ;dist3d(sv #0, #1) (cm);arb. units", 500, 0, 1);
  h_svdist2d_both_absdeltaphi01_nsv2_no_shared_jets = fs->make<TH1F>("h_svdist2d_both_absdeltaphi01_nsv2_no_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5 ;dist2d(sv #0, #1) (cm);arb. units", 500, 0, 1);    
  h_svdist3d_both_absdeltaphi01_nsv2_no_shared_jets = fs->make<TH1F>("h_svdist3d_both_absdeltaphi01_nsv2_no_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5 ;dist3d(sv #0, #1) (cm);arb. units", 500, 0, 1);   
}

void MFVVertexHistos::analyze(const edm::Event& event, const edm::EventSetup&) {

  unsigned run = event.id().run();
  unsigned lumi = event.luminosityBlock();
  unsigned long long evt = event.id().event();
  if ( 42000 <= evt && evt <= 42115 ){
  edm::Handle<MFVEvent> mevent;
  event.getByToken(mevent_token, mevent);
  edm::Handle<double> weight;
  event.getByToken(weight_token, weight); 
  const double w = *weight;
  h_w->Fill(w);

  const double bsx = mevent->bsx;
  const double bsy = mevent->bsy;
  const double bsz = mevent->bsz;
  const math::XYZPoint bs(bsx, bsy, bsz);
  const math::XYZPoint pv(mevent->pvx, mevent->pvy, mevent->pvz);

  edm::Handle<MFVVertexAuxCollection> auxes;
  event.getByToken(vertex_token, auxes);
  
  const int nsv = int(auxes->size());
  h_nsv->Fill(nsv, w);

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
			  //absdelta = double(fabs(reco::deltaPhi(phin, mevent->jet_phi[jet_index])));
			  //absdeltaphi_sv_jets.push_back(absdelta);
          	  }
	  }
	  sv_track_which_jet.push_back(track_which_jet);
	  sv_track_which_idx.push_back(track_which_idx);
	  /*
	  if (absdeltaphi_sv_jets.size() > 0 ){
              h_max_absdeltaphi_sv_jets->Fill(*max_element(absdeltaphi_sv_jets.begin(), absdeltaphi_sv_jets.end()),w);
              
	  }
	  */
	  
   }
/*
  std::vector<int> sv1_test_track_which_idx = sv_track_which_idx[1];
  std::vector<int> sv0_test_track_which_idx = sv_track_which_idx[0];
  for (unsigned i = 0; i < sv0_test_track_which_idx.size(); ++i){

       for (unsigned j = 0; j < sv1_test_track_which_idx.size(); ++j){
              if (sv0_test_track_which_idx[i] == sv1_test_track_which_idx[j]){
                   std::cout << "found at" << sv0_test_track_which_idx[i] << std::endl;
              }
           
}

}
*/




  if (nsv >= 2) {
    const MFVVertexAux& sv0 = auxes->at(0);
    const MFVVertexAux& sv1 = auxes->at(1);
    double svdist2d = mag(sv0.x - sv1.x, sv0.y - sv1.y);
    double svdist3d = mag(sv0.x - sv1.x, sv0.y - sv1.y, sv0.z - sv1.z);
    h_svdist2d->Fill(svdist2d, w);
    h_svdist3d->Fill(svdist3d, w);
    h_sv0pvdz_v_sv1pvdz->Fill(sv0.pvdz(), sv1.pvdz(), w);
    h_sv0pvdzsig_v_sv1pvdzsig->Fill(sv0.pvdzsig(), sv1.pvdzsig(), w);
    double phi0 = atan2(sv0.y - bsy, sv0.x - bsx);
    double phi1 = atan2(sv1.y - bsy, sv1.x - bsx);
    h_absdeltaphi01->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);
    h_fractrackssharedwpv01 ->Fill(float(sv0.ntrackssharedwpv () + sv1.ntrackssharedwpv ())/(sv0.ntracks() + sv1.ntracks()), w);
    h_fractrackssharedwpvs01->Fill(float(sv0.ntrackssharedwpvs() + sv1.ntrackssharedwpvs())/(sv0.ntracks() + sv1.ntracks()), w);
    h_pvmosttracksshared->Fill(sv0.ntrackssharedwpvs() ? sv0.pvmosttracksshared() : -1,
                               sv1.ntrackssharedwpvs() ? sv1.pvmosttracksshared() : -1,
                               w);

    

    bool shared_jet = std::find_first_of (sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()) != sv_track_which_jet[0].end();
    h_sv_shared_jets->Fill(shared_jet, w);
    if (shared_jet) {
      h_svdist2d_shared_jets->Fill(svdist2d, w);
      h_absdeltaphi01_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);
      
	  int nsharedjets=1;
	  std::vector<int> nsharedjet_jet_index;
	  std::vector<std::vector<int> > sv_track_which_jet_copy;
	  sv_track_which_jet_copy = sv_track_which_jet;
	  std::vector<int> nsharedjet_tracks_sv0;
          std::vector<int> nsharedjet_tracks_sv1;
          std::vector<std::vector<int> >sv0_sharedjet_which_idx;
          std::vector<std::vector<int> >sv1_sharedjet_which_idx;
          std::vector<int>::iterator it = std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end());
          int idx = std::distance(sv_track_which_jet_copy[0].begin(),it);
	  int jet_index = sv_track_which_jet_copy[0].at(idx);
	  nsharedjet_jet_index.push_back(jet_index);
	  sv_track_which_jet_copy[0].erase(std::remove(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), jet_index), sv_track_which_jet_copy[0].end());
	  sv_track_which_jet_copy[1].erase(std::remove(sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end(), jet_index), sv_track_which_jet_copy[1].end());	  
	  nsharedjet_tracks_sv0.push_back(std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index)); 
	  
	  std::vector<int> sv0_track_which_jet = sv_track_which_jet[0];
	  std::vector<int> sv0_track_which_idx = sv_track_which_idx[0];
	  std::vector<int> sv0_track_which_temp_idx;
	  std::multimap<int, size_t> sv0_m_nshj1;

	  for (size_t k = 0; k < sv0_track_which_jet.size(); k++) if (sv0_track_which_jet[k] == jet_index ) { sv0_m_nshj1.insert({ sv0_track_which_jet[k], k }); }

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
          sv1_sharedjet_which_idx.push_back(sv1_track_which_temp_idx);                                                                                                                        sv1_track_which_temp_idx = {};                                                                                                                                            
      // std::cout << "shared-jet #" << nsharedjets << " has shared-jet ntracks from sv#0 =" << std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index) << ", shared-jet ntracks from sv#1 =" << std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index) << std::endl;

	  while (std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end()) != sv_track_which_jet_copy[0].end()) {
		  nsharedjets ++;
		  it = std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end());
          idx = std::distance(sv_track_which_jet_copy[0].begin(),it);
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
                  sv0_sharedjet_which_idx.push_back(sv0_track_which_temp_idx);                                                                                                                        sv0_track_which_temp_idx = {};

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
                  sv1_sharedjet_which_idx.push_back(sv1_track_which_temp_idx);                                                                                                                        sv1_track_which_temp_idx = {};
       //  std::cout << "shared-jet #" << nsharedjets << " has shared-jet ntracks from sv#0 =" << std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index) << ", shared-jet ntracks from sv#1 =" << std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index) << std::endl;

	  }
      //std::cout << "# of shared jets are " << nsharedjets << std::endl;
      //std::cout << "sv#0 has ntracks="<< sv_track_which_jet[0].size() << ", sv#1 has ntracks="<< sv_track_which_jet[1].size() << std::endl;
      h_nsharedjets_shared_jets->Fill(nsharedjets);
	  for (int i = 0; i < nsharedjets; i++) {
		  if (nsharedjet_tracks_sv0[i] > nsharedjet_tracks_sv1[i]) {
			  h_ratio_ntracks_shared_jets->Fill(nsharedjet_tracks_sv0[i] / nsharedjet_tracks_sv1[i]);
		  }
		  if (nsharedjet_tracks_sv0[i] <= nsharedjet_tracks_sv1[i]) {
			  h_ratio_ntracks_shared_jets->Fill(nsharedjet_tracks_sv1[i] / nsharedjet_tracks_sv0[i]);
		  }
		  if (nsharedjets == 1) {
			  h_ntracks_sv0sv1_nsharedjet1->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]);
			  if (nsharedjet_tracks_sv1[i] == 1) { h_ntracks_nshj1_sv0_ntrack1_sv1->Fill(nsharedjet_tracks_sv0[i]); }
			  if (nsharedjet_tracks_sv1[i] == 2) { h_ntracks_nshj1_sv0_ntrack2_sv1->Fill(nsharedjet_tracks_sv0[i]); }
			  if (nsharedjet_tracks_sv1[i] == 3) { h_ntracks_nshj1_sv0_ntrack3_sv1->Fill(nsharedjet_tracks_sv0[i]); }
			  if (nsharedjet_tracks_sv1[i] == 4) { h_ntracks_nshj1_sv0_ntrack4_sv1->Fill(nsharedjet_tracks_sv0[i]); }
			  if (nsharedjet_tracks_sv1[i] == 5) { h_ntracks_nshj1_sv0_ntrack5_sv1->Fill(nsharedjet_tracks_sv0[i]); }

		  }
		  if (nsharedjets == 2) { h_ntracks_sv0sv1_nsharedjet2->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
		  if (nsharedjets == 3) { h_ntracks_sv0sv1_nsharedjet3->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
		  if (nsharedjets == 4) { h_ntracks_sv0sv1_nsharedjet4->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
		  if (nsharedjets == 5) { h_ntracks_sv0sv1_nsharedjet5->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
		  if (nsharedjet_tracks_sv1[i] == 1) { h_ntracks_sv0_ntrack1_sv1->Fill(nsharedjet_tracks_sv0[i]); }
		  if (nsharedjet_tracks_sv1[i] == 2) { h_ntracks_sv0_ntrack2_sv1->Fill(nsharedjet_tracks_sv0[i]); }
		  if (nsharedjet_tracks_sv1[i] == 3) { h_ntracks_sv0_ntrack3_sv1->Fill(nsharedjet_tracks_sv0[i]); }
		  if (nsharedjet_tracks_sv1[i] == 4) { h_ntracks_sv0_ntrack4_sv1->Fill(nsharedjet_tracks_sv0[i]); }
		  if (nsharedjet_tracks_sv1[i] == 5) { h_ntracks_sv0_ntrack5_sv1->Fill(nsharedjet_tracks_sv0[i]); }

	  }
      if (nsv==2){
		  h_nsharedjets_nsv2_shared_jets->Fill(nsharedjets);
		  std::vector<double> absdeltaphi_sv0_shared_jets;
		  std::vector<double> absdeltaphi_sv1_shared_jets;
		  std::vector<double> max_dphi_sv0_sv1;
		  for (int i = 0; i < nsharedjets; i++) {
			  //TODO: add histo with max of max dphi from sv#1 and sv#0		
				  //    : add histo with dphi of another sv and such a jet index

			  int ntracks_sv0 = nsharedjet_tracks_sv0[i];
			  int ntracks_sv1 = nsharedjet_tracks_sv1[i];
			  std::vector<double> absdeltaphi_sv_shared_jets;
			  for (int j = 0; j < ntracks_sv0; ++j) {
				  jet_index = nsharedjet_jet_index[i];
				  double absdelta_sv0 = double(fabs(reco::deltaPhi(phi0, mevent->jet_phi[jet_index])));
				  absdeltaphi_sv0_shared_jets.push_back(absdelta_sv0);
			  }
			  for (int j = 0; j < ntracks_sv1; ++j) {
				  jet_index = nsharedjet_jet_index[i];
				  double absdelta_sv1 = double(fabs(reco::deltaPhi(phi1, mevent->jet_phi[jet_index])));
				  absdeltaphi_sv1_shared_jets.push_back(absdelta_sv1);
			  }
                          

			  if (nsharedjet_tracks_sv0[i] > nsharedjet_tracks_sv1[i]) {
				  h_ratio_ntracks_nsv2_shared_jets->Fill(nsharedjet_tracks_sv0[i] / nsharedjet_tracks_sv1[i]);
			  }
			  if (nsharedjet_tracks_sv0[i] <= nsharedjet_tracks_sv1[i]) {
				  h_ratio_ntracks_nsv2_shared_jets->Fill(nsharedjet_tracks_sv1[i] / nsharedjet_tracks_sv0[i]);
			  }
                          if (nsharedjets == 1){
				  h_ntracks_sv0sv1_nsv2_nsharedjet1->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]);
				  if (nsharedjet_tracks_sv1[i] == 1) { h_ntracks_nshj1_sv0_ntrack1_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
				  if (nsharedjet_tracks_sv1[i] == 2) { h_ntracks_nshj1_sv0_ntrack2_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
				  if (nsharedjet_tracks_sv1[i] == 3) { h_ntracks_nshj1_sv0_ntrack3_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
				  if (nsharedjet_tracks_sv1[i] == 4) { h_ntracks_nshj1_sv0_ntrack4_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
				  if (nsharedjet_tracks_sv1[i] == 5) { h_ntracks_nshj1_sv0_ntrack5_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }

				  
			  }
			  if (nsharedjets == 2) { h_ntracks_sv0sv1_nsv2_nsharedjet2->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
			  if (nsharedjets == 3) { h_ntracks_sv0sv1_nsv2_nsharedjet3->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
			  if (nsharedjets == 4) { h_ntracks_sv0sv1_nsv2_nsharedjet4->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
			  if (nsharedjets == 5) { h_ntracks_sv0sv1_nsv2_nsharedjet5->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
			  if (nsharedjet_tracks_sv1[i] == 1) { h_ntracks_sv0_ntrack1_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
			  if (nsharedjet_tracks_sv1[i] == 2) { h_ntracks_sv0_ntrack2_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
			  if (nsharedjet_tracks_sv1[i] == 3) { h_ntracks_sv0_ntrack3_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
			  if (nsharedjet_tracks_sv1[i] == 4) { h_ntracks_sv0_ntrack4_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
			  if (nsharedjet_tracks_sv1[i] == 5) { h_ntracks_sv0_ntrack5_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }

		  }
		  double max_dphi_sv0 = *std::max_element(absdeltaphi_sv0_shared_jets.begin(), absdeltaphi_sv0_shared_jets.end());
		  int max_dphi_sv0_idx = std::max_element(absdeltaphi_sv0_shared_jets.begin(), absdeltaphi_sv0_shared_jets.end()) - absdeltaphi_sv0_shared_jets.begin();
		  double max_dphi_sv1 = *std::max_element(absdeltaphi_sv1_shared_jets.begin(), absdeltaphi_sv1_shared_jets.end());
		  int max_dphi_sv1_idx = std::max_element(absdeltaphi_sv1_shared_jets.begin(), absdeltaphi_sv1_shared_jets.end()) - absdeltaphi_sv1_shared_jets.begin();
		  max_dphi_sv0_sv1.push_back(max_dphi_sv0);
		  max_dphi_sv0_sv1.push_back(max_dphi_sv1);


		  if (nsharedjets == 1) {
			  h_max_absdeltaphi0_sv_nshj1_shared_jets->Fill(*std::max_element(max_dphi_sv0_sv1.begin(), max_dphi_sv0_sv1.end()), w);
			  if (max_dphi_sv0 > max_dphi_sv1) {
				  h_max_absdeltaphi1_sv_nshj1_shared_jets->Fill(absdeltaphi_sv1_shared_jets[max_dphi_sv0_idx], w);

			  }
			  else {
				  h_max_absdeltaphi1_sv_nshj1_shared_jets->Fill(absdeltaphi_sv0_shared_jets[max_dphi_sv1_idx], w);
			  }
			
		  }

		  h_max_absdeltaphi0_sv_shared_jets->Fill(*std::max_element(max_dphi_sv0_sv1.begin(), max_dphi_sv0_sv1.end()), w);
		  if (max_dphi_sv0 > max_dphi_sv1) {
			  h_max_absdeltaphi1_sv_shared_jets->Fill(absdeltaphi_sv1_shared_jets[max_dphi_sv0_idx], w);
		  }
		  else {
			  h_max_absdeltaphi1_sv_shared_jets->Fill(absdeltaphi_sv0_shared_jets[max_dphi_sv1_idx], w);
		  }

		  h_lspdist2d_nsv2_shared_jets->Fill(mevent->lspdist2d(), w);
		  h_lspdist3d_nsv2_shared_jets->Fill(mevent->lspdist3d(), w);
		  h_absdeltaphi01_genlsp_nsv2_shared_jets->Fill(std::abs(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1])), w);

		  
          h_absdeltaphi01_nsv2_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);
		  
          if (fabs(reco::deltaPhi(phi0, phi1)) <= 0.5){
             h_svdist2d_small_absdeltaphi01_nsv2_shared_jets->Fill(svdist2d,w);
             h_svdist3d_small_absdeltaphi01_nsv2_shared_jets->Fill(svdist3d,w);
			 
			 // start split vertex
			 h_nsharedjets_small_nsv2_shared_jets->Fill(nsharedjets);
			 std::vector<double> absdeltaphi_small_sv0_shared_jets;
			 std::vector<double> absdeltaphi_small_sv1_shared_jets;
			 std::vector<double> max_small_dphi_sv0_sv1;
			 for (int i = 0; i < nsharedjets; i++) {
				 //TODO: add histo with max of max dphi from sv#1 and sv#0		
					 //    : add histo with dphi of another sv and such a jet index

		                 jet_index = nsharedjet_jet_index[i];
		                 double absdelta_small_sv0 = double(fabs(reco::deltaPhi(phi0, mevent->jet_phi[jet_index])));
				 absdeltaphi_small_sv0_shared_jets.push_back(absdelta_small_sv0);
				 
			         double absdelta_small_sv1 = double(fabs(reco::deltaPhi(phi1, mevent->jet_phi[jet_index])));
			         absdeltaphi_small_sv1_shared_jets.push_back(absdelta_small_sv1);
                                 if (nsharedjets == 1){

				 if (nsharedjet_tracks_sv0[i] > nsharedjet_tracks_sv1[i]) {
					 h_ratio_ntracks_small_nsv2_nshj1_shared_jets->Fill(nsharedjet_tracks_sv0[i] / nsharedjet_tracks_sv1[i]);		   
				 }
				 if (nsharedjet_tracks_sv0[i] <= nsharedjet_tracks_sv1[i]) {
					 h_ratio_ntracks_small_nsv2_nshj1_shared_jets->Fill(nsharedjet_tracks_sv1[i] / nsharedjet_tracks_sv0[i]);			 
				 }
			
					 h_ntracks_sv0sv1_small_nsv2_nsharedjet1->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]);				
					 if (nsharedjet_tracks_sv1[i] == 1) { h_ntracks_small_nshj1_sv0_ntrack1_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
					 if (nsharedjet_tracks_sv1[i] == 2) { h_ntracks_small_nshj1_sv0_ntrack2_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
					 if (nsharedjet_tracks_sv1[i] == 3) { h_ntracks_small_nshj1_sv0_ntrack3_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
					 if (nsharedjet_tracks_sv1[i] == 4) { h_ntracks_small_nshj1_sv0_ntrack4_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
					 if (nsharedjet_tracks_sv1[i] == 5) { h_ntracks_small_nshj1_sv0_ntrack5_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }


				 }
				 if (nsharedjets == 2) { h_ntracks_sv0sv1_small_nsv2_nsharedjet2->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
				 if (nsharedjets == 3) { h_ntracks_sv0sv1_small_nsv2_nsharedjet3->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
				 if (nsharedjets == 4) { h_ntracks_sv0sv1_small_nsv2_nsharedjet4->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
				 if (nsharedjets == 5) { h_ntracks_sv0sv1_small_nsv2_nsharedjet5->Fill(nsharedjet_tracks_sv0[i], nsharedjet_tracks_sv1[i]); }
				 if (nsharedjet_tracks_sv1[i] == 1) { h_ntracks_small_sv0_ntrack1_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
				 if (nsharedjet_tracks_sv1[i] == 2) { h_ntracks_small_sv0_ntrack2_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
				 if (nsharedjet_tracks_sv1[i] == 3) { h_ntracks_small_sv0_ntrack3_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
				 if (nsharedjet_tracks_sv1[i] == 4) { h_ntracks_small_sv0_ntrack4_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }
				 if (nsharedjet_tracks_sv1[i] == 5) { h_ntracks_small_sv0_ntrack5_sv1_nsv2->Fill(nsharedjet_tracks_sv0[i]); }

			 }
			 double max_small_dphi_sv0 = *std::max_element(absdeltaphi_small_sv0_shared_jets.begin(), absdeltaphi_small_sv0_shared_jets.end());
			 int max_small_dphi_sv0_idx = std::max_element(absdeltaphi_small_sv0_shared_jets.begin(), absdeltaphi_small_sv0_shared_jets.end()) - absdeltaphi_small_sv0_shared_jets.begin();
			 double max_small_dphi_sv1 = *std::max_element(absdeltaphi_small_sv1_shared_jets.begin(), absdeltaphi_small_sv1_shared_jets.end());
			 int max_small_dphi_sv1_idx = std::max_element(absdeltaphi_small_sv1_shared_jets.begin(), absdeltaphi_small_sv1_shared_jets.end()) - absdeltaphi_small_sv1_shared_jets.begin();
			 max_small_dphi_sv0_sv1.push_back(max_small_dphi_sv0);
			 max_small_dphi_sv0_sv1.push_back(max_small_dphi_sv1);

			 if (nsharedjets == 1) {
				 h_max_absdeltaphi0_small_sv_nshj1_shared_jets->Fill(*std::max_element(max_small_dphi_sv0_sv1.begin(), max_small_dphi_sv0_sv1.end()), w);
				 if (max_small_dphi_sv0 > max_small_dphi_sv1) {
					 h_max_absdeltaphi1_small_sv_nshj1_shared_jets->Fill(absdeltaphi_small_sv1_shared_jets[max_small_dphi_sv0_idx], w);
				 }
				 else {
					 h_max_absdeltaphi1_small_sv_nshj1_shared_jets->Fill(absdeltaphi_small_sv0_shared_jets[max_small_dphi_sv1_idx], w);
				 }
			 }

			 h_max_absdeltaphi0_small_sv_shared_jets->Fill(*std::max_element(max_small_dphi_sv0_sv1.begin(), max_small_dphi_sv0_sv1.end()), w);
			 if (max_small_dphi_sv0 > max_small_dphi_sv1) {
				 h_max_absdeltaphi1_small_sv_shared_jets->Fill(absdeltaphi_small_sv1_shared_jets[max_small_dphi_sv0_idx], w);
			 }
			 else {
				 h_max_absdeltaphi1_small_sv_shared_jets->Fill(absdeltaphi_small_sv0_shared_jets[max_small_dphi_sv1_idx], w);
			 }

			 // end split-vertex 



            }
		  else {

			// start not split vertex
			 h_nsharedjets_large_nsv2_shared_jets->Fill(nsharedjets, w);
			 h_svdist2d_large_absdeltaphi01_nsv2_shared_jets->Fill(svdist2d, w);
			 h_svdist3d_large_absdeltaphi01_nsv2_shared_jets->Fill(svdist3d, w);
                                                                                                                                                                                   
                   std::vector<double> absdeltaphi_large_sv0_shared_jets;                                                                                                                              std::vector<double> absdeltaphi_large_sv1_shared_jets;                                                                                                                              
                   std::vector<double> max_large_dphi_sv0_sv1;                                                                                                                                         for (int i = 0; i < nsharedjets; i++) {                                                                                                                                                                                                                                                                                                                                        jet_index = nsharedjet_jet_index[i];                                                                                                                                                double absdelta_large_sv0 = double(fabs(reco::deltaPhi(phi0, mevent->jet_phi[jet_index])));                                                                                         absdeltaphi_large_sv0_shared_jets.push_back(absdelta_large_sv0);                                                                                                                                                                                                                                                                                                        double absdelta_large_sv1 = double(fabs(reco::deltaPhi(phi1, mevent->jet_phi[jet_index])));                                                                                         absdeltaphi_large_sv1_shared_jets.push_back(absdelta_large_sv1);
                   }
                   
                   double max_large_dphi_sv0 = *std::max_element(absdeltaphi_large_sv0_shared_jets.begin(), absdeltaphi_large_sv0_shared_jets.end());                                                  int max_large_dphi_sv0_idx = std::max_element(absdeltaphi_large_sv0_shared_jets.begin(), absdeltaphi_large_sv0_shared_jets.end()) - absdeltaphi_large_sv0_shared_jets.begin();                                                                                                                                                                                          double max_large_dphi_sv1 = *std::max_element(absdeltaphi_large_sv1_shared_jets.begin(), absdeltaphi_large_sv1_shared_jets.end());                                                  int max_large_dphi_sv1_idx = std::max_element(absdeltaphi_large_sv1_shared_jets.begin(), absdeltaphi_large_sv1_shared_jets.end()) - absdeltaphi_large_sv1_shared_jets.begin();                                                                                                                                                                                          max_large_dphi_sv0_sv1.push_back(max_large_dphi_sv0);                                                                                                                               max_large_dphi_sv0_sv1.push_back(max_large_dphi_sv1);
                   
                   if (nsharedjets == 1){
                         std::cout << "run " << run << " lumi " << lumi << " event " << evt << std::endl; 
                         h_max_absdeltaphi0_large_sv_nshj1_shared_jets->Fill(*std::max_element(max_large_dphi_sv0_sv1.begin(), max_large_dphi_sv0_sv1.end()), w);                                            if (max_large_dphi_sv0 > max_large_dphi_sv1) {                                                                                                                                              h_max_absdeltaphi1_large_sv_nshj1_shared_jets->Fill(absdeltaphi_large_sv1_shared_jets[max_large_dphi_sv0_idx], w);                                                          }                                                                                                                                                                                   else {                                                                                                                                                                                      h_max_absdeltaphi1_large_sv_nshj1_shared_jets->Fill(absdeltaphi_large_sv0_shared_jets[max_large_dphi_sv1_idx], w);                                                          }

                         if (nsharedjet_tracks_sv0[0] > nsharedjet_tracks_sv1[0]) {                                                                                                                                  h_ratio_ntracks_large_nsv2_nshj1_shared_jets->Fill(nsharedjet_tracks_sv0[0] / nsharedjet_tracks_sv1[0]);                                                                          }                                                                                                                                                                                   if (nsharedjet_tracks_sv0[0] <= nsharedjet_tracks_sv1[0]) {                                                                                                                                 h_ratio_ntracks_large_nsv2_nshj1_shared_jets->Fill(nsharedjet_tracks_sv1[0] / nsharedjet_tracks_sv0[0]);                                                                          }
                         		
                         if (max_dphi_sv0 > max_dphi_sv1) {
				 std::vector<double> absdeltaphi_min_sv1_shared_tracks;
                                 std::vector<int> sv1_nsharedjets1_which_idx = sv1_sharedjet_which_idx[0];
				 for (int j = 0; j < nsharedjet_tracks_sv1[0];j++){
					 int track_idx = sv1_nsharedjets1_which_idx[j];
					 double absdelta_min_sv1_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx]))); //phi1
                                         double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx])));
                                         h_absdeltaphi_large_jet_shared_tracks_nshj1_nsv2->Fill(absdelta_jet_track);
                                        std::cout << "sv0>sv1 : phi1 is " << phi1 << " with track phi " << sv1.track_phi[track_idx] << ", deltaPhi(jet,trk) is " << double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx]))) << std::endl;
					 absdeltaphi_min_sv1_shared_tracks.push_back(absdelta_min_sv1_track);
				 }
				 double min_dphi_sv1_track = *std::min_element(absdeltaphi_min_sv1_shared_tracks.begin(), absdeltaphi_min_sv1_shared_tracks.end());
				 int min_dphi_sv1_track_idx = std::min_element(absdeltaphi_min_sv1_shared_tracks.begin(), absdeltaphi_min_sv1_shared_tracks.end()) - absdeltaphi_min_sv1_shared_tracks.begin();

                                 int min_sv1_track_idx = sv1_nsharedjets1_which_idx[min_dphi_sv1_track_idx];
				 h_max_absdeltaphi1_large_sv_nshj1_shared_tracks->Fill(double(fabs(reco::deltaPhi(phi1, sv1.track_phi[min_sv1_track_idx]))), w);
				 h_max_pt_absdeltaphi1_large_sv_nshj1_shared_tracks->Fill(sv1.track_pt(min_sv1_track_idx), w);

				  std::cout << "sv0>sv1 : phi1 is " << phi1 << " with min track phi " << sv1.track_phi[min_sv1_track_idx] << ", shj is at" << mevent->jet_phi[nsharedjet_jet_index[0]] << ", deltaPhi(sv1,trk) is " << double(fabs(reco::deltaPhi(phi1, sv1.track_phi[min_sv1_track_idx]))) << std::endl;

                                 std::vector<int> sv0_nsharedjets1_which_idx = sv0_sharedjet_which_idx[0];
				 std::vector<double> absdeltaphi_max_sv0_shared_tracks;
				 for (int j = 0; j < nsharedjet_tracks_sv0[0]; j++) {
					 int track_idx = sv0_nsharedjets1_which_idx[j];
					 double absdelta_max_sv0_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx]))); //phi0
                                         double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx])));                                                                     h_absdeltaphi_large_jet_shared_tracks_nshj1_nsv2->Fill(absdelta_jet_track);
                                         std::cout << "sv0>sv1 : phi0 is " << phi0 << " with track phi " << sv0.track_phi[track_idx] << ", deltaPhi(jet,trk) is " << double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx]))) << std::endl;
					 absdeltaphi_max_sv0_shared_tracks.push_back(absdelta_max_sv0_track);
				 }
				 double max_dphi_sv0_track = *std::max_element(absdeltaphi_max_sv0_shared_tracks.begin(), absdeltaphi_max_sv0_shared_tracks.end());
				 int max_dphi_sv0_track_idx = std::max_element(absdeltaphi_max_sv0_shared_tracks.begin(), absdeltaphi_max_sv0_shared_tracks.end()) - absdeltaphi_max_sv0_shared_tracks.begin();

				 int max_sv0_track_idx = sv0_nsharedjets1_which_idx[max_dphi_sv0_track_idx];
                                 h_max_absdeltaphi0_large_sv_nshj1_shared_tracks->Fill(double(fabs(reco::deltaPhi(phi0, sv0.track_phi[max_sv0_track_idx]))), w); 
				 h_max_pt_absdeltaphi0_large_sv_nshj1_shared_tracks->Fill(sv0.track_pt(max_sv0_track_idx), w);

                                std::cout << "sv0>sv1 : phi0 is " << phi0 << " with max track phi " << sv0.track_phi[max_sv0_track_idx] << ", shj is at" << mevent->jet_phi[nsharedjet_jet_index[0]] << ", detaPhi(sv0,trk) is " << double(fabs(reco::deltaPhi(phi0, sv0.track_phi[max_sv0_track_idx]))) << std::endl; 
				 
				 AlgebraicVector3 mom_max(sv0.track_px[max_sv0_track_idx], sv0.track_py[max_sv0_track_idx], sv0.track_pz[max_sv0_track_idx]);
				 AlgebraicVector3 ref_max(sv0.track_vx[max_sv0_track_idx], sv0.track_vy[max_sv0_track_idx], sv0.track_vz[max_sv0_track_idx]);
				 AlgebraicVector3 mom_min(sv1.track_px[min_sv1_track_idx], sv1.track_py[min_sv1_track_idx], sv1.track_pz[min_sv1_track_idx]);
                 AlgebraicVector3 ref_min(sv1.track_vx[min_sv1_track_idx], sv1.track_vy[min_sv1_track_idx], sv1.track_vz[min_sv1_track_idx]);

				 Measurement1D miss_dist_max = miss_dist(sv1,ref_max,mom_max);
				 Measurement1D miss_dist_2D_max = miss_dist_2D(sv1, ref_max, mom_max);
				 h_miss_dist_absdeltaphi0_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_max.value(), w);
				 h_miss_dist_err_absdeltaphi0_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_max.error(), w);
                                 h_miss_dist_2D_absdeltaphi0_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_2D_max.value(), w);
				 h_miss_dist_sig_absdeltaphi0_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_max.significance(), w);
				 Measurement1D miss_dist_min = miss_dist(sv0,ref_min,mom_min); 
				 Measurement1D miss_dist_2D_min = miss_dist_2D(sv0, ref_min, mom_min);
				 h_miss_dist_absdeltaphi1_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_min.value(), w);
				 h_miss_dist_err_absdeltaphi1_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_min.error(), w);
                                 h_miss_dist_2D_absdeltaphi1_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_2D_min.value(), w);
				 h_miss_dist_sig_absdeltaphi1_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_min.significance(), w);
				 
				 Measurement1D miss_dist_max_own = miss_dist(sv0, ref_max, mom_max);
				 Measurement1D miss_dist_2D_max_own = miss_dist_2D(sv0, ref_max, mom_max);
				 h_miss_dist_absdeltaphi0_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_max_own.value(), w);
				 h_miss_dist_err_absdeltaphi0_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_max_own.error(), w);
                                 h_miss_dist_2D_absdeltaphi0_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_2D_max_own.value(), w);
				 h_miss_dist_sig_absdeltaphi0_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_max_own.significance(), w);
				 Measurement1D miss_dist_min_own = miss_dist(sv1, ref_min, mom_min);
				 Measurement1D miss_dist_2D_min_own = miss_dist_2D(sv1, ref_min, mom_min);
				 h_miss_dist_absdeltaphi1_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_min_own.value(), w);
				 h_miss_dist_err_absdeltaphi1_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_min_own.error(), w); 
                                 h_miss_dist_2D_absdeltaphi1_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_2D_min_own.value(), w);
				 h_miss_dist_sig_absdeltaphi1_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_min_own.significance(), w);

				 h_vertex_chi2dof_absdeltaphi0_large_sv2_nshj1->Fill(sv0.chi2dof(),w);



			 }
			 else {

				 std::vector<double> absdeltaphi_min_sv0_shared_tracks;
                                 std::vector<int> sv0_nsharedjets1_which_idx = sv0_sharedjet_which_idx[0];
				 for (int j = 0; j < nsharedjet_tracks_sv0[0]; j++) {
					 int track_idx = sv0_nsharedjets1_which_idx[j];
					 double absdelta_min_sv0_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx]))); //phi0
                                         double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx])));                                                                     h_absdeltaphi_large_jet_shared_tracks_nshj1_nsv2->Fill(absdelta_jet_track);
					 absdeltaphi_min_sv0_shared_tracks.push_back(absdelta_min_sv0_track);
                                         std::cout << "sv1>sv0 : phi0 is " << phi0 << " with track phi " << sv0.track_phi[track_idx] <<", deltaPhi(jet,trk) is " << double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv0.track_phi[track_idx]))) << std::endl;
				 }
				 double min_dphi_sv0_track = *std::min_element(absdeltaphi_min_sv0_shared_tracks.begin(), absdeltaphi_min_sv0_shared_tracks.end());
				 int min_dphi_sv0_track_idx = std::min_element(absdeltaphi_min_sv0_shared_tracks.begin(), absdeltaphi_min_sv0_shared_tracks.end()) - absdeltaphi_min_sv0_shared_tracks.begin();

				 int min_sv0_track_idx = sv0_nsharedjets1_which_idx[min_dphi_sv0_track_idx];
                                 h_max_absdeltaphi1_large_sv_nshj1_shared_tracks->Fill(double(fabs(reco::deltaPhi(phi0, sv0.track_phi[min_sv0_track_idx]))), w);
				 h_max_pt_absdeltaphi1_large_sv_nshj1_shared_tracks->Fill(sv0.track_pt(min_sv0_track_idx), w);
                                 
                                 std::cout << "sv1>sv0 : phi0 is " << phi0 << " with min track phi " << sv0.track_phi[min_sv0_track_idx] << ", shj is at" << mevent->jet_phi[nsharedjet_jet_index[0]] << ", deltaPhi(sv0,trk) is " << double(fabs(reco::deltaPhi(phi0, sv0.track_phi[min_sv0_track_idx]))) << std::endl;
                                 std::vector<int> sv1_nsharedjets1_which_idx = sv1_sharedjet_which_idx[0];
				 std::vector<double> absdeltaphi_max_sv1_shared_tracks;
				 for (int j = 0; j < nsharedjet_tracks_sv1[0]; j++) {
					 int track_idx = sv1_nsharedjets1_which_idx[j];
					 double absdelta_max_sv1_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx]))); //phi1
                                         double absdelta_jet_track = double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx])));                                                                     h_absdeltaphi_large_jet_shared_tracks_nshj1_nsv2->Fill(absdelta_jet_track);
					 absdeltaphi_max_sv1_shared_tracks.push_back(absdelta_max_sv1_track);
                                         std::cout << "sv1>sv0 : phi1 is " << phi1 << " with track phi " << sv1.track_phi[track_idx] <<", deltaPhi(jet,trk) is " << double(fabs(reco::deltaPhi(mevent->jet_phi[jet_index], sv1.track_phi[track_idx]))) << std::endl;
				 }
				 double max_dphi_sv1_track = *std::max_element(absdeltaphi_max_sv1_shared_tracks.begin(), absdeltaphi_max_sv1_shared_tracks.end());
				 int max_dphi_sv1_track_idx = std::max_element(absdeltaphi_max_sv1_shared_tracks.begin(), absdeltaphi_max_sv1_shared_tracks.end()) - absdeltaphi_max_sv1_shared_tracks.begin();
				 

				 int max_sv1_track_idx = sv1_nsharedjets1_which_idx[max_dphi_sv1_track_idx];
                                 h_max_absdeltaphi0_large_sv_nshj1_shared_tracks->Fill(double(fabs(reco::deltaPhi(phi1, sv1.track_phi[max_sv1_track_idx]))), w);
				 h_max_pt_absdeltaphi0_large_sv_nshj1_shared_tracks->Fill(sv1.track_pt(max_sv1_track_idx), w);
                                 
                                std::cout << "sv1>sv0 : phi1 is " << phi1 << " with max track phi " << sv1.track_phi[max_sv1_track_idx] << ", shj is at" << mevent->jet_phi[nsharedjet_jet_index[0]] <<", deltaPhi(sv1,trk) is " << double(fabs(reco::deltaPhi(phi1, sv1.track_phi[max_sv1_track_idx]))) << std::endl;
				 
				 AlgebraicVector3 mom_max(sv1.track_px[max_sv1_track_idx], sv1.track_py[max_sv1_track_idx], sv1.track_pz[max_sv1_track_idx]);
				 AlgebraicVector3 ref_max(sv1.track_vx[max_sv1_track_idx], sv1.track_vy[max_sv1_track_idx], sv1.track_vz[max_sv1_track_idx]); 
				 AlgebraicVector3 mom_min(sv0.track_px[min_sv0_track_idx], sv0.track_py[min_sv0_track_idx], sv0.track_pz[min_sv0_track_idx]);
                 AlgebraicVector3 ref_min(sv0.track_vx[min_sv0_track_idx], sv0.track_vy[min_sv0_track_idx], sv0.track_vz[min_sv0_track_idx]); 
                                 
			    Measurement1D miss_dist_max = miss_dist(sv0, ref_max, mom_max);
				Measurement1D miss_dist_2D_max = miss_dist_2D(sv0, ref_max, mom_max);
				h_miss_dist_absdeltaphi0_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_max.value(), w);
				h_miss_dist_err_absdeltaphi0_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_max.error(), w);
                                h_miss_dist_2D_absdeltaphi0_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_2D_max.value(), w);
				h_miss_dist_sig_absdeltaphi0_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_max.significance(), w);
                Measurement1D miss_dist_min = miss_dist(sv1, ref_min, mom_min);
				Measurement1D miss_dist_2D_min = miss_dist_2D(sv1, ref_min, mom_min);
				h_miss_dist_absdeltaphi1_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_min.value(), w);
			        h_miss_dist_err_absdeltaphi1_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_min.error(), w); 
                        	h_miss_dist_2D_absdeltaphi1_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_2D_min.value(), w);
				h_miss_dist_sig_absdeltaphi1_large_other_sv_nshj1_shared_tracks->Fill(miss_dist_min.significance(), w);

				Measurement1D miss_dist_max_own = miss_dist(sv1, ref_max, mom_max);
				Measurement1D miss_dist_2D_max_own = miss_dist_2D(sv1, ref_max, mom_max);
			    h_miss_dist_absdeltaphi0_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_max_own.value(), w);
		                h_miss_dist_err_absdeltaphi0_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_max_own.error(), w); 
                 		h_miss_dist_2D_absdeltaphi0_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_2D_max_own.value(), w);
				h_miss_dist_sig_absdeltaphi0_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_max_own.significance(), w);
				Measurement1D miss_dist_min_own = miss_dist(sv0, ref_min, mom_min);
				Measurement1D miss_dist_2D_min_own = miss_dist_2D(sv0, ref_min, mom_min);
			    h_miss_dist_absdeltaphi1_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_min_own.value(), w);
				h_miss_dist_err_absdeltaphi1_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_min_own.error(), w);
                                h_miss_dist_2D_absdeltaphi1_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_2D_min_own.value(), w);
				h_miss_dist_sig_absdeltaphi1_large_its_sv_nshj1_shared_tracks->Fill(miss_dist_min_own.significance(), w);

				h_vertex_chi2dof_absdeltaphi0_large_sv2_nshj1->Fill(sv1.chi2dof(), w);

				//std::cout << max_sv1_track_idx << "!=" << min_sv0_track_idx << std::endl;
			 }
			 
			 
                          
                        }    
			//end not split vertex 
		  }
      }
       
    } else {

      h_nsharedjets_shared_jets->Fill((int)0, w);
      h_svdist2d_no_shared_jets->Fill(svdist2d, w);
      h_absdeltaphi01_no_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);
	  h_lspdist2d_nsv2_no_shared_jets->Fill(mevent->lspdist2d(), w);
	  h_lspdist3d_nsv2_no_shared_jets->Fill(mevent->lspdist3d(), w);
	  h_absdeltaphi01_genlsp_nsv2_no_shared_jets->Fill(std::abs(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1])), w);
      if (nsv==2){
                 h_nsharedjets_nsv2_shared_jets->Fill((int)0, w);
                 h_absdeltaphi01_nsv2_no_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);                                                                                                     h_svdist2d_both_absdeltaphi01_nsv2_no_shared_jets->Fill(svdist2d,w);                                                                                                                h_svdist3d_both_absdeltaphi01_nsv2_no_shared_jets->Fill(svdist3d,w);                                                                                                          
                 if ((reco::deltaPhi(phi0, phi1)) <= 0.5) {h_nsharedjets_small_nsv2_shared_jets->Fill((int)0, w);}
                 else{h_nsharedjets_large_nsv2_shared_jets->Fill((int)0, w);}
         }
    }

	
 }  
 }
  }

DEFINE_FWK_MODULE(MFVVertexHistos);
