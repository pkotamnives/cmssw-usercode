#include "TH2.h"
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
#include <iostream>
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
  TH1F* h_max_absdeltaphi1_sv_nshj1_shared_jets;
  TH1F* h_max_absdeltaphi0_sv_shared_jets;
  TH1F* h_max_absdeltaphi1_sv_shared_jets;

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

  TH1F* h_svdist2d_small_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_svdist3d_small_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_svdist2d_large_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_svdist3d_large_absdeltaphi01_nsv2_shared_jets;
  TH1F* h_svdist2d_small_absdeltaphi01_no_shared_jets;
  TH1F* h_svdist3d_small_absdeltaphi01_no_shared_jets;
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
  h_max_absdeltaphi0_sv_nshj1_shared_jets = fs->make<TH1F>("h_max_absdeltaphi0_sv_nshj1_shared_jets", "nsv = 2, nsharedjets = 1;max( dphi(sv0, only-one shared jet), dphi(sv1, only-one shared jet) );arb. units", 316, 0, 3.16);
  h_max_absdeltaphi1_sv_nshj1_shared_jets = fs->make<TH1F>("h_max_absdeltaphi1_sv_nshj1_shared_jets", "nsv = 2, nsharedjets = 1;abs(delta(phi of max-dphi shared jet, phi of another sv ));arb. units", 316, 0, 3.16);
  h_max_absdeltaphi0_sv_shared_jets = fs->make<TH1F>("h_max_absdeltaphi0_sv_shared_jets", "nsv = 2;max( dphi(sv0, shared jets), dphi(sv1, shared jets) );arb. units", 316, 0, 3.16);
  h_max_absdeltaphi1_sv_shared_jets = fs->make<TH1F>("h_max_absdeltaphi1_sv_shared_jets", "nsv = 2;abs(delta(phi of max-dphi shared jet, phi of another sv ));arb. units", 316, 0, 3.16);
  
  h_lspdist2d_nsv2_shared_jets = fs->make<TH1F>("h_lspdist2d_nsv2_shared_jets", "nsv = 2;dist2d(gen vtx #0, #1) (cm)", 200, 0, 2);
  h_lspdist3d_nsv2_shared_jets = fs->make<TH1F>("h_lspdist3d_nsv2_shared_jets", " nsv = 2;dist3d(gen vtx #0, #1) (cm)", 200, 0, 2);
  h_absdeltaphi01_genlsp_nsv2_shared_jets = fs->make<TH1F>("h_absdeltaphi01_genlsp_nsv2_shared_jets", "nsv = 2;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);
  h_lspdist2d_nsv2_no_shared_jets = fs->make<TH1F>("h_lspdist2d_nsv2_no_shared_jets", "nsv = 2;dist2d(gen vtx #0, #1) (cm)", 200, 0, 2);
  h_lspdist3d_nsv2_no_shared_jets = fs->make<TH1F>("h_lspdist3d_nsv2_no_shared_jets", " nsv = 2;dist3d(gen vtx #0, #1) (cm)", 200, 0, 2);
  h_absdeltaphi01_genlsp_nsv2_no_shared_jets = fs->make<TH1F>("h_absdeltaphi01_genlsp_nsv2_no_shared_jets", "nsv = 2;abs(delta(gen lsp phi of sv #0, gen lsp phi of sv #1));arb. units", 316, 0, 3.16);

  h_nsharedjets_shared_jets = fs->make<TH1F>("h_nsharedjets_shared_jets", "nsv >= 2;# of shared jets;arb. units", 10, 0, 10);
  h_ratio_ntracks_shared_jets = fs->make<TH1F>("h_ratio_ntracks_shared_jets", "nsv >= 2;ratios of shared jets (>=1);arb. units", 50, 0, 10);
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
  h_ratio_ntracks_nsv2_shared_jets = fs->make<TH1F>("h_ratio_ntracks_nsv2_shared_jets", "nsv = 2;ratios of shared jets (>=1);arb. units", 50, 0, 10);
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

  
  h_svdist2d_small_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist2d_small_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5 ;dist2d(sv #0, #1) (cm);arb. units", 500, 0, 0.1); 
  h_svdist3d_small_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist3d_small_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 <= 0.5 ;dist3d(sv #0, #1) (cm);arb. units", 500, 0, 0.1); 
  h_svdist2d_large_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist2d_large_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 > 0.5 ;dist2d(sv #0, #1) (cm);arb. units", 500, 0, 0.1);
  h_svdist3d_large_absdeltaphi01_nsv2_shared_jets = fs->make<TH1F>("h_svdist3d_large_absdeltaphi01_nsv2_shared_jets", "nsv = 2, absdeltaphi01 > 0.5 ;dist3d(sv #0, #1) (cm);arb. units", 500, 0, 0.1);
  h_svdist2d_small_absdeltaphi01_no_shared_jets = fs->make<TH1F>("h_svdist2d_small_absdeltaphi01_no_shared_jets", "nsv >= 2, absdeltaphi01 <= 0.5 ;dist2d(sv #0, #1) (cm);arb. units", 500, 0, 0.1);    
  h_svdist3d_small_absdeltaphi01_no_shared_jets = fs->make<TH1F>("h_svdist3d_small_absdeltaphi01_no_shared_jets", "nsv >= 2, absdeltaphi01 <= 0.5 ;dist3d(sv #0, #1) (cm);arb. units", 500, 0, 0.1);   
}

void MFVVertexHistos::analyze(const edm::Event& event, const edm::EventSetup&) {
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
  std::vector<std::vector<int> > sv_track_which_jet;
  for (int isv = 0; isv < nsv; ++isv) {			//loop over vertices
	  const MFVVertexAux& aux = auxes->at(isv);
	  const int ntracks = aux.ntracks();

	  //double phin = atan2(aux.y - bsy, aux.x - bsx);

	  std::vector<int> track_which_jet;
	  std::vector<double> absdeltaphi_sv_jets;
	  for (int i = 0; i < ntracks; ++i) {		   //loop over tracks associated with a vertex
		  double match_threshold = 1.3;
		  int jet_index = 255;
		  double absdelta = 0;
                  for (unsigned j = 0; j < mevent->jet_track_which_jet.size(); ++j) {		   //loop over jets associated with all tracks
			  double a = fabs(aux.track_pt(i) - fabs(mevent->jet_track_qpt[j])) + 1;
			  double b = fabs(aux.track_eta[i] - mevent->jet_track_eta[j]) + 1;
			  double c = fabs(aux.track_phi[i] - mevent->jet_track_phi[j]) + 1;
			  if (a * b * c < match_threshold) {
				  match_threshold = a * b * c;
				  jet_index = mevent->jet_track_which_jet[j];
			  }
		  }
		  if (jet_index != 255) {
			  track_which_jet.push_back((int)jet_index);	  // get a valid set of jets for a track
			  //absdelta = double(fabs(reco::deltaPhi(phin, mevent->jet_phi[jet_index])));
			  //absdeltaphi_sv_jets.push_back(absdelta);
          	  }
	  }
	  sv_track_which_jet.push_back(track_which_jet);
	  /*
	  if (absdeltaphi_sv_jets.size() > 0 ){
              h_max_absdeltaphi_sv_jets->Fill(*max_element(absdeltaphi_sv_jets.begin(), absdeltaphi_sv_jets.end()),w);
              
	  }
	  */
	  
   }







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
          std::vector<int>::iterator it = std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end());
          int idx = std::distance(sv_track_which_jet_copy[0].begin(),it);
	  int jet_index = sv_track_which_jet_copy[0].at(idx);
	  nsharedjet_jet_index.push_back(jet_index);
	  sv_track_which_jet_copy[0].erase(std::remove(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), jet_index), sv_track_which_jet_copy[0].end());
	  sv_track_which_jet_copy[1].erase(std::remove(sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end(), jet_index), sv_track_which_jet_copy[1].end());	  
	  nsharedjet_tracks_sv0.push_back(std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index));                                                                 nsharedjet_tracks_sv1.push_back(std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index));
      // std::cout << "shared-jet #" << nsharedjets << " has shared-jet ntracks from sv#0 =" << std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index) << ", shared-jet ntracks from sv#1 =" << std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index) << std::endl;

	  while (std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end()) != sv_track_which_jet_copy[0].end()) {
		  nsharedjets ++;
		  it = std::find_first_of(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end());
          idx = std::distance(sv_track_which_jet_copy[0].begin(),it);
          jet_index = sv_track_which_jet_copy[0].at(idx);
		  nsharedjet_jet_index.push_back(jet_index);
	      sv_track_which_jet_copy[0].erase(std::remove(sv_track_which_jet_copy[0].begin(), sv_track_which_jet_copy[0].end(), jet_index), sv_track_which_jet_copy[0].end());
	      sv_track_which_jet_copy[1].erase(std::remove(sv_track_which_jet_copy[1].begin(), sv_track_which_jet_copy[1].end(), jet_index), sv_track_which_jet_copy[1].end());	  
              nsharedjet_tracks_sv0.push_back(std::count(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), jet_index));
              nsharedjet_tracks_sv1.push_back(std::count(sv_track_which_jet[1].begin(), sv_track_which_jet[1].end(), jet_index));	   
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
		  std::vector<double> max_dphi_sv0_sv1_nshj1;
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
			  if (nsharedjets == 1) {
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
		  double max_dphi_sv0 = *max_element(absdeltaphi_sv0_shared_jets.begin(), absdeltaphi_sv0_shared_jets.end());
		  int max_dphi_sv0_idx = max_element(absdeltaphi_sv0_shared_jets.begin(), absdeltaphi_sv0_shared_jets.end()) - absdeltaphi_sv0_shared_jets.begin();
		  double max_dphi_sv1 = *max_element(absdeltaphi_sv1_shared_jets.begin(), absdeltaphi_sv1_shared_jets.end());
		  int max_dphi_sv1_idx = max_element(absdeltaphi_sv1_shared_jets.begin(), absdeltaphi_sv1_shared_jets.end()) - absdeltaphi_sv1_shared_jets.begin();
		  max_dphi_sv0_sv1.push_back(max_dphi_sv0_nshj1);
		  max_dphi_sv0_sv1.push_back(max_dphi_sv1_nshj1);

		  if (nsharedjets == 1) {
			  h_max_absdeltaphi0_sv_nshj1_shared_jets->Fill(*max_element(max_dphi_sv0_sv1.begin(), max_dphi_sv0_sv1.end()), w);
			  if (max_dphi_sv0 > max_dphi_sv1) {
				  h_max_absdeltaphi1_sv_nshj1_shared_jets->Fill(absdeltaphi_sv1_shared_jets[max_dphi_sv0_idx], w);
			  }
			  else {
				  h_max_absdeltaphi1_sv_nshj1_shared_jets->Fill(absdeltaphi_sv0_shared_jets[max_dphi_sv1_idx], w);
			  }
		  }

		  h_max_absdeltaphi0_sv_shared_jets->Fill(*max_element(max_dphi_sv0_sv1.begin(), max_dphi_sv0_sv1.end()), w);
		  if (max_dphi_sv0 > max_dphi_sv1) {
			  h_max_absdeltaphi1_sv_shared_jets->Fill(absdeltaphi_sv1_shared_jets[max_dphi_sv0_idx], w);
		  }
		  else {
			  h_max_absdeltaphi1_sv_shared_jets->Fill(absdeltaphi_sv0_shared_jets[max_dphi_sv1_idx], w);
		  }

		  h_lspdist2d_nsv2_shared_jets->Fill(mevent->lspdist2d(), w);
		  h_lspdist3d_nsv2_shared_jets->Fill(mevent->lspdist3d(), w);
		  h_absdeltaphi01_genlsp_nsv2_shared_jets->Fill(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1]), w);

		  
          h_absdeltaphi01_nsv2_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);
		  
          if (fabs(reco::deltaPhi(phi0, phi1)) <= 0.5){
             h_svdist2d_small_absdeltaphi01_nsv2_shared_jets->Fill(svdist2d,w);
             h_svdist3d_small_absdeltaphi01_nsv2_shared_jets->Fill(svdist3d,w);
            }
		  else {
			 h_svdist2d_large_absdeltaphi01_nsv2_shared_jets->Fill(svdist2d, w);
			 h_svdist3d_large_absdeltaphi01_nsv2_shared_jets->Fill(svdist3d, w);
		  }
      }
       
    } else {
      h_svdist2d_no_shared_jets->Fill(svdist2d, w);
      h_absdeltaphi01_no_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);
	  h_lspdist2d_nsv2_no_shared_jets->Fill(mevent->lspdist2d(), w);
	  h_lspdist3d_nsv2_no_shared_jets->Fill(mevent->lspdist3d(), w);
	  h_absdeltaphi01_genlsp_nsv2_no_shared_jets->Fill(reco::deltaPhi(mevent->gen_lsp_phi[0], mevent->gen_lsp_phi[1]), w);
      if (nsv==2){                                                                                                                                                                            h_absdeltaphi01_nsv2_no_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);                                                                                           
          if (fabs(reco::deltaPhi(phi0, phi1)) <= 0.5){                                                                                                                                          h_svdist2d_small_absdeltaphi01_no_shared_jets->Fill(svdist2d,w);                                                                                                                    h_svdist3d_small_absdeltaphi01_no_shared_jets->Fill(svdist3d,w);                                                                                                                      } 
          }
    }

	
  }
}

DEFINE_FWK_MODULE(MFVVertexHistos);
