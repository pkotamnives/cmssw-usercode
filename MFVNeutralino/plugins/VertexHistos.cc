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
  void fill(PairwiseHistos* hs, const int, const PairwiseHistos::ValueMap& val, const double weight) const { hs[sv_all].Fill(val, -1, weight); }

  PairwiseHistos h_sv[sv_num_indices];

  TH1F* h_sv_jets_deltaphi[4][sv_num_indices];

  TH2F* h_sv_bs2derr_bsbs2ddist[sv_num_indices];
  TH2F* h_pvrho_bsbs2ddist[sv_num_indices];

  TH1F* h_w;
  TH1F* h_nsv;
  TH2F* h_sv_xy;
  TH2F* h_sv_yz;
  TH2F* h_sv_xz;
  TH2F* h_sv_rz;
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

  TH1F* h_sv_track_weight[sv_num_indices];
  TH1F* h_sv_track_q[sv_num_indices];
  TH1F* h_sv_track_pt[sv_num_indices];
  TH1F* h_sv_track_eta[sv_num_indices];
  TH1F* h_sv_track_phi[sv_num_indices];
  TH1F* h_sv_track_dxy[sv_num_indices];
  TH1F* h_sv_track_dz[sv_num_indices];
  TH1F* h_sv_track_pt_err[sv_num_indices];
  TH1F* h_sv_track_eta_err[sv_num_indices];
  TH1F* h_sv_track_phi_err[sv_num_indices];
  TH1F* h_sv_track_dxy_err[sv_num_indices];
  TH1F* h_sv_track_dz_err[sv_num_indices];
  TH1F* h_sv_track_nsigmadxy[sv_num_indices];
  TH1F* h_sv_track_chi2dof[sv_num_indices];
  TH1F* h_sv_track_npxhits[sv_num_indices];
  TH1F* h_sv_track_nsthits[sv_num_indices];
  TH1F* h_sv_track_nhitsbehind[sv_num_indices];
  TH1F* h_sv_track_nhitslost[sv_num_indices];
  TH1F* h_sv_track_nhits[sv_num_indices];
  TH1F* h_sv_track_injet[sv_num_indices];
  TH1F* h_sv_track_inpv[sv_num_indices];
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

  PairwiseHistos::HistoDefs hs;

  hs.add("x", "SV x (cm)", 100, -4, 4);
  hs.add("y", "SV y (cm)", 100, -4, 4);
  hs.add("z", "SV z (cm)", 100, -25, 25);
  hs.add("phi", "SV phi", 25, -3.15, 3.15);
  hs.add("phi_pv", "SV phi w.r.t. PV", 25, -3.15, 3.15);
  hs.add("cxx", "SV covariance xx (cm^{2})", 100, 0, 1e-5);
  hs.add("cxy", "SV covariance xy (cm^{2})", 100, -1e-5, 1e-5);
  hs.add("cxz", "SV covariance xz (cm^{2})", 100, -1e-5, 1e-5);
  hs.add("cyy", "SV covariance yy (cm^{2})", 100, 0, 1e-5);
  hs.add("cyz", "SV covariance yz (cm^{2})", 100, -1e-5, 1e-5);
  hs.add("czz", "SV covariance zz (cm^{2})", 100, 0, 1e-5);

  hs.add("rescale_chi2", "rescaled-fit SV x (cm)", 40, 0, 10);
  hs.add("rescale_x", "rescaled-fit SV x (cm)", 100, -4, 4);
  hs.add("rescale_y", "rescaled-fit SV y (cm)", 100, -4, 4);
  hs.add("rescale_z", "rescaled-fit SV z (cm)", 100, -25, 25);
  hs.add("rescale_cxx", "rescaled-fit SV covariance xx (cm^{2})", 100, 0, 1e-5);
  hs.add("rescale_cxy", "rescaled-fit SV covariance xy (cm^{2})", 100, -1e-5, 1e-5);
  hs.add("rescale_cxz", "rescaled-fit SV covariance xz (cm^{2})", 100, -1e-5, 1e-5);
  hs.add("rescale_cyy", "rescaled-fit SV covariance yy (cm^{2})", 100, 0, 1e-5);
  hs.add("rescale_cyz", "rescaled-fit SV covariance yz (cm^{2})", 100, -1e-5, 1e-5);
  hs.add("rescale_czz", "rescaled-fit SV covariance zz (cm^{2})", 100, 0, 1e-5);
  hs.add("rescale_dx", "rescaled-fit - nominal SV x (cm)", 100, -5e-4, 5e-4);
  hs.add("rescale_dy", "rescaled-fit - nominal SV y (cm)", 100, -5e-4, 5e-4);
  hs.add("rescale_dz", "rescaled-fit - nominal SV z (cm)", 100, -5e-4, 5e-4);
  hs.add("rescale_dx_big", "rescaled-fit - nominal SV x (cm)", 100, -4, 4);
  hs.add("rescale_dy_big", "rescaled-fit - nominal SV y (cm)", 100, -4, 4);
  hs.add("rescale_dz_big", "rescaled-fit - nominal SV z (cm)", 100, -4, 4);
  hs.add("rescale_d2", "rescaled-fit - nominal SV (2D) (cm)", 100, 0, 8e-4);
  hs.add("rescale_d2_big", "rescaled-fit - nominal SV (2D) (cm)", 100, 0, 4);
  hs.add("rescale_d3", "rescaled-fit - nominal SV (3D) (cm)", 100, 0, 1e-3);
  hs.add("rescale_d3_big", "rescaled-fit - nominal SV (3D) (cm)", 100, 0, 4);
  hs.add("rescale_bsbs2ddist", "rescaled-fit d_{BV} (cm)", 500, 0, 2.5);
  hs.add("rescale_bs2derr", "rescaled-fit #sigma(dist2d(SV, beamspot)) (cm)", 1000, 0, 0.05);

  hs.add("max_nm1_refit_dist3_wbad", "maximum n-1 refit distance (3D) (cm)", 1001, -0.001, 1);
  hs.add("max_nm1_refit_dist3", "maximum n-1 refit distance (3D) (cm)", 1000, 0, 1);
  hs.add("max_nm1_refit_dist2", "maximum n-1 refit distance (2D) (cm)", 1000, 0, 1);
  hs.add("max_nm1_refit_distz", "maximum n-1 refit z distance (cm)", 1000, 0, 1);

  hs.add("nlep", "# leptons", 10, 0, 10);

  hs.add("ntracks",                       "# of tracks/SV",                                                               40,    0,      40);
  hs.add("ntracksptgt3",                  "# of tracks/SV w/ p_{T} > 3 GeV",                                              40,    0,      40);
  hs.add("ntracksptgt10",                 "# of tracks/SV w/ p_{T} > 10 GeV",                                             40,    0,      40);
  hs.add("ntracksetagt1p5",               "# of tracks/SV w/ |#eta| > 1.5",                                               40,    0,      40);
  hs.add("trackminnhits",                 "min number of hits on track per SV",                                           40,    0,      40);
  hs.add("trackmaxnhits",                 "max number of hits on track per SV",                                           40,    0,      40);
  hs.add("njetsntks",                     "# of jets assoc. by tracks to SV",                                             10,    0,      10);
  hs.add("chi2dof",                       "SV #chi^2/dof",                                                                50,    0,       7);
  hs.add("chi2dofprob",                   "SV p(#chi^2, dof)",                                                            50,    0,       1.2);

  hs.add("tkonlyp",                       "SV tracks-only p (GeV)",                                                       50,    0,     500);
  hs.add("tkonlypt",                      "SV tracks-only p_{T} (GeV)",                                                   50,    0,     400);
  hs.add("tkonlyeta",                     "SV tracks-only #eta",                                                          50,   -4,       4);
  hs.add("tkonlyrapidity",                "SV tracks-only rapidity",                                                      50,   -4,       4);
  hs.add("tkonlyphi",                     "SV tracks-only #phi",                                                          50,   -3.15,    3.15);
  hs.add("tkonlymass",                    "SV tracks-only mass (GeV)",                                                   100,    0,    1000);

  hs.add("jetsntkp",                      "SV jets-by-ntracks -only p (GeV)",                                             50,    0,    1000);
  hs.add("jetsntkpt",                     "SV jets-by-ntracks -only p_{T} (GeV)",                                         50,    0,    1000);
  hs.add("jetsntketa",                    "SV jets-by-ntracks -only #eta",                                                50,   -4,       4);
  hs.add("jetsntkrapidity",               "SV jets-by-ntracks -only rapidity",                                            50,   -4,       4);
  hs.add("jetsntkphi",                    "SV jets-by-ntracks -only #phi",                                                50,   -3.15,    3.15);
  hs.add("jetsntkmass",                   "SV jets-by-ntracks -only mass (GeV)",                                          50,    0,    2000);

  hs.add("tksjetsntkp",                   "SV tracks-plus-jets-by-ntracks p (GeV)",                                       50,    0,    1000);
  hs.add("tksjetsntkpt",                  "SV tracks-plus-jets-by-ntracks p_{T} (GeV)",                                   50,    0,    1000);
  hs.add("tksjetsntketa",                 "SV tracks-plus-jets-by-ntracks #eta",                                          50,   -4,       4);
  hs.add("tksjetsntkrapidity",            "SV tracks-plus-jets-by-ntracks rapidity",                                      50,   -4,       4);
  hs.add("tksjetsntkphi",                 "SV tracks-plus-jets-by-ntracks #phi",                                          50,   -3.15,    3.15);
  hs.add("tksjetsntkmass",                "SV tracks-plus-jets-by-ntracks mass (GeV)",                                   100,    0,    5000);
				        
  hs.add("costhtkonlymombs",              "cos(angle(2-momentum (tracks-only), 2-dist to BS))",                           21,   -1,       1.1);
  hs.add("costhtkonlymompv2d",            "cos(angle(2-momentum (tracks-only), 2-dist to PV))",                           21,   -1,       1.1);
  hs.add("costhtkonlymompv3d",            "cos(angle(3-momentum (tracks-only), 3-dist to PV))",                           21,   -1,       1.1);

  hs.add("costhtksjetsntkmombs",          "cos(angle(2-momentum (tracks-plus-jets-by-ntracks), 2-dist to BS))",          21,   -1,       1.1);
  hs.add("costhtksjetsntkmompv2d",        "cos(angle(2-momentum (tracks-plus-jets-by-ntracks), 2-dist to PV))",          21,   -1,       1.1);
  hs.add("costhtksjetsntkmompv3d",        "cos(angle(3-momentum (tracks-plus-jets-by-ntracks), 3-dist to PV))",          21,   -1,       1.1);

  hs.add("missdisttkonlypv",              "miss dist. (tracks-only) of SV to PV (cm)",                                   100,    0,       2);
  hs.add("missdisttkonlypverr",           "#sigma(miss dist. (tracks-only) of SV to PV) (cm)",                           100,    0,       0.05);
  hs.add("missdisttkonlypvsig",           "N#sigma(miss dist. (tracks-only) of SV to PV) (cm)",                          100,    0,     100);

  hs.add("missdisttksjetsntkpv",          "miss dist. (tracks-plus-jets-by-ntracks) of SV to PV (cm)",                   100,    0,       2);
  hs.add("missdisttksjetsntkpverr",       "#sigma(miss dist. (tracks-plus-jets-by-ntracks) of SV to PV) (cm)",           100,    0,       0.05);
  hs.add("missdisttksjetsntkpvsig",       "N#sigma(miss dist. (tracks-plus-jets-by-ntracks) of SV to PV) (cm)",          100,    0,     100);
					  
  hs.add("sumpt2",                        "SV #Sigma p_{T}^{2} (GeV^2)",                                                  50,    0,    10000);

  hs.add("ntrackssharedwpv",  "number of tracks shared with the PV", 30, 0, 30);
  hs.add("ntrackssharedwpvs", "number of tracks shared with any PV", 30, 0, 30);
  hs.add("fractrackssharedwpv",  "fraction of tracks shared with the PV", 41, 0, 1.025);
  hs.add("fractrackssharedwpvs", "fraction of tracks shared with any PV", 41, 0, 1.025);
  hs.add("npvswtracksshared", "number of PVs having tracks shared",  30, 0, 30);
  
  hs.add("trackdxymin", "SV min{trk_{i} dxy(BS)} (cm)", 50, 0, 0.2);
  hs.add("trackdxymax", "SV max{trk_{i} dxy(BS)} (cm)", 50, 0, 2);
  hs.add("trackdxyavg", "SV avg{trk_{i} dxy(BS)} (cm)", 50, 0, 0.5);
  hs.add("trackdxyrms", "SV rms{trk_{i} dxy(BS)} (cm)", 50, 0, 0.5);

  hs.add("trackdzmin", "SV min{trk_{i} dz(PV)} (cm)", 50, 0, 0.5);
  hs.add("trackdzmax", "SV max{trk_{i} dz(PV)} (cm)", 50, 0, 2);
  hs.add("trackdzavg", "SV avg{trk_{i} dz(PV)} (cm)", 50, 0, 1);
  hs.add("trackdzrms", "SV rms{trk_{i} dz(PV)} (cm)", 50, 0, 0.5);

  hs.add("trackpterrmin", "SV min{frac. #sigma trk_{i} p_{T}}", 32, 0, 2);
  hs.add("trackpterrmax", "SV max{frac. #sigma trk_{i} p_{T}}", 32, 0, 2);
  hs.add("trackpterravg", "SV avg{frac. #sigma trk_{i} p_{T}}", 32, 0, 2);
  hs.add("trackpterrrms", "SV rms{frac. #sigma trk_{i} p_{T}}", 32, 0, 2);

  hs.add("tracketaerrmin", "SV min{frac. #sigma trk_{i} #eta}", 32, 0, 0.002);
  hs.add("tracketaerrmax", "SV max{frac. #sigma trk_{i} #eta}", 32, 0, 0.005);
  hs.add("tracketaerravg", "SV avg{frac. #sigma trk_{i} #eta}", 32, 0, 0.002);
  hs.add("tracketaerrrms", "SV rms{frac. #sigma trk_{i} #eta}", 32, 0, 0.002);

  hs.add("trackphierrmin", "SV min{frac. #sigma trk_{i} #phi}", 32, 0, 0.002);
  hs.add("trackphierrmax", "SV max{frac. #sigma trk_{i} #phi}", 32, 0, 0.005);
  hs.add("trackphierravg", "SV avg{frac. #sigma trk_{i} #phi}", 32, 0, 0.002);
  hs.add("trackphierrrms", "SV rms{frac. #sigma trk_{i} #phi}", 32, 0, 0.002);

  hs.add("trackdxyerrmin", "SV min{#sigma trk_{i} dxy(BS)} (cm)", 32, 0, 0.004);
  hs.add("trackdxyerrmax", "SV max{#sigma trk_{i} dxy(BS)} (cm)", 32, 0, 0.1);
  hs.add("trackdxyerravg", "SV avg{#sigma trk_{i} dxy(BS)} (cm)", 32, 0, 0.1);
  hs.add("trackdxyerrrms", "SV rms{#sigma trk_{i} dxy(BS)} (cm)", 32, 0, 0.1);

  hs.add("trackdzerrmin", "SV min{#sigma trk_{i} dz(PV)} (cm)", 32, 0, 0.01);
  hs.add("trackdzerrmax", "SV max{#sigma trk_{i} dz(PV)} (cm)", 32, 0, 0.1);
  hs.add("trackdzerravg", "SV avg{#sigma trk_{i} dz(PV)} (cm)", 32, 0, 0.1);
  hs.add("trackdzerrrms", "SV rms{#sigma trk_{i} dz(PV)} (cm)", 32, 0, 0.1);

  hs.add("trackpairdetamin", "SV min{#Delta #eta(i,j)}", 150,    0,       1.5);
  hs.add("trackpairdetamax", "SV max{#Delta #eta(i,j)}", 150,    0,       7);
  hs.add("trackpairdetaavg", "SV avg{#Delta #eta(i,j)}", 150,    0,       5);
  hs.add("trackpairdetarms", "SV rms{#Delta #eta(i,j)}", 150,    0,       3);

  hs.add("trackpairdphimax",   "SV max{|#Delta #phi(i,j)|}",   100, 0, 3.15);
  hs.add("trackpairdphimaxm1", "SV max-1{|#Delta #phi(i,j)|}", 100, 0, 3.15);
  hs.add("trackpairdphimaxm2", "SV max-2{|#Delta #phi(i,j)|}", 100, 0, 3.15);

  hs.add("trackpairdrmin", "SV min{#Delta R(i,j)}", 150, 0, 1.5);
  hs.add("trackpairdrmax", "SV max{#Delta R(i,j)}", 150, 0, 7);
  hs.add("trackpairdravg", "SV avg{#Delta R(i,j)}", 150, 0, 5);
  hs.add("trackpairdrrms", "SV rms{#Delta R(i,j)}", 150, 0, 3);

  hs.add("costhtkmomvtxdispmin", "SV min{cos(angle(trk_{i}, SV-PV))}", 50, -1, 1);
  hs.add("costhtkmomvtxdispmax", "SV max{cos(angle(trk_{i}, SV-PV))}", 50, -1, 1);
  hs.add("costhtkmomvtxdispavg", "SV avg{cos(angle(trk_{i}, SV-PV))}", 50, -1, 1);
  hs.add("costhtkmomvtxdisprms", "SV rms{cos(angle(trk_{i}, SV-PV))}", 50,  0, 1);

  hs.add("costhjetmomvtxdispmin", "SV min{cos(angle(jet_{i}, SV-PV))}", 50, -1, 1);
  hs.add("costhjetmomvtxdispmax", "SV max{cos(angle(jet_{i}, SV-PV))}", 50, -1, 1);
  hs.add("costhjetmomvtxdispavg", "SV avg{cos(angle(jet_{i}, SV-PV))}", 50, -1, 1);
  hs.add("costhjetmomvtxdisprms", "SV rms{cos(angle(jet_{i}, SV-PV))}", 50,  0, 1);

  hs.add("multipv_maxdz",    "max #Delta z of PV w tracks shared (cm)", 100, 0, 10);
  hs.add("multipvbyz_maxdz", "max #Delta z of PV w track-assoc-by-z (cm)", 100, 0, 10);

  hs.add("gen2ddist",                     "dist2d(SV, closest gen vtx) (cm)",                                            200,    0,       0.2);
  hs.add("gen2derr",                      "#sigma(dist2d(SV, closest gen vtx)) (cm)",                                    200,    0,       0.2);
  hs.add("gen2dsig",                      "N#sigma(dist2d(SV, closest gen vtx)) (cm)",                                   200,    0,     100);
  hs.add("gen3ddist",                     "dist3d(SV, closest gen vtx) (cm)",                                            200,    0,       0.2);
  hs.add("gen3derr",                      "#sigma(dist3d(SV, closest gen vtx)) (cm)",                                    200,    0,       0.2);
  hs.add("gen3dsig",                      "N#sigma(dist3d(SV, closest gen vtx)) (cm)",                                   200,    0,     100);
  hs.add("bs2ddist",                      "dist2d(SV, beamspot) (cm)",                                                   500,    0,      2.5);
  hs.add("bsbs2ddist",                    "dist2d(SV, beamspot) (cm)",                                                   500,    0,      2.5);
  hs.add("bs2derr",                       "#sigma(dist2d(SV, beamspot)) (cm)",                                           1000,    0,       0.05);
  hs.add("bs2dsig",                       "N#sigma(dist2d(SV, beamspot))",                                               100,    0,     100);
  hs.add("pv2ddist",                      "dist2d(SV, PV) (cm)",                                                         100,    0,       0.5);
  hs.add("pv2derr",                       "#sigma(dist2d(SV, PV)) (cm)",                                                 100,    0,       0.05);
  hs.add("pv2dsig",                       "N#sigma(dist2d(SV, PV))",                                                     100,    0,     100);
  hs.add("pv3ddist",                      "dist3d(SV, PV) (cm)",                                                         100,    0,       0.5);
  hs.add("pv3derr",                       "#sigma(dist3d(SV, PV)) (cm)",                                                 100,    0,       0.1);
  hs.add("pv3dsig",                       "N#sigma(dist3d(SV, PV))",                                                     100,    0,     100);
  hs.add("pvdz",                          "dz(SV, PV) (cm)",                                                             100,    0,       0.5);
  hs.add("pvdzerr",                       "#sigma(dz(SV, PV)) (cm)",                                                     100,    0,       0.1);
  hs.add("pvdzsig",                       "N#sigma(dz(SV, PV))",                                                         100,    0,     100);

  const char* lmt_ex[4] = {"", "loose b-", "medium b-", "tight b-"};
  for (int i = 0; i < 4; ++i) {
    hs.add(TString::Format("jet%d_deltaphi0", i), TString::Format("|#Delta#phi| to closest %sjet", lmt_ex[i]),      25, 0, 3.15);
    hs.add(TString::Format("jet%d_deltaphi1", i), TString::Format("|#Delta#phi| to next closest %sjet", lmt_ex[i]), 25, 0, 3.15);
  }

  for (int i = 0; i < max_ntrackplots; ++i) {
    hs.add(TString::Format("track%i_weight",        i), TString::Format("track%i weight",                      i),  21,  0,      1.05);
    hs.add(TString::Format("track%i_q",             i), TString::Format("track%i charge",                      i),   4, -2,      2);
    hs.add(TString::Format("track%i_pt",            i), TString::Format("track%i p_{T} (GeV)",                 i), 200,  0,    200);
    hs.add(TString::Format("track%i_eta",           i), TString::Format("track%i #eta",                        i),  50, -4,      4);
    hs.add(TString::Format("track%i_phi",           i), TString::Format("track%i #phi",                        i),  50, -3.15,   3.15);
    hs.add(TString::Format("track%i_dxy",           i), TString::Format("track%i dxy (cm)",                    i), 100,  0,      1);
    hs.add(TString::Format("track%i_dz",            i), TString::Format("track%i dz (cm)",                     i), 100,  0,      1);
    hs.add(TString::Format("track%i_pt_err",        i), TString::Format("track%i #sigma(p_{T})/p_{T}",         i), 200,  0,      2);
    hs.add(TString::Format("track%i_eta_err",       i), TString::Format("track%i #sigma(#eta)",                i), 200,  0,      0.02);
    hs.add(TString::Format("track%i_phi_err",       i), TString::Format("track%i #sigma(#phi)",                i), 200,  0,      0.02);
    hs.add(TString::Format("track%i_dxy_err",       i), TString::Format("track%i #sigma(dxy) (cm)",            i), 100,  0,      0.1);
    hs.add(TString::Format("track%i_dz_err",        i), TString::Format("track%i #sigma(dz) (cm)",             i), 100,  0,      0.1);
    hs.add(TString::Format("track%i_nsigmadxy",     i), TString::Format("track%i n#sigma(dxy)",                i), 400,  0,     40);
    hs.add(TString::Format("track%i_chi2dof",       i), TString::Format("track%i #chi^{2}/dof",                i), 100,  0,     10);
    hs.add(TString::Format("track%i_npxhits",       i), TString::Format("track%i number of pixel hits",        i),  12,  0,     12);
    hs.add(TString::Format("track%i_nsthits",       i), TString::Format("track%i number of strip hits",        i),  28,  0,     28);
    hs.add(TString::Format("track%i_nhitsbehind",   i), TString::Format("track%i number of hits behind",       i),  10,  0,     10);
    hs.add(TString::Format("track%i_nhitslost",     i), TString::Format("track%i number of hits lost",         i),  10,  0,     10);
    hs.add(TString::Format("track%i_nhits",         i), TString::Format("track%i number of hits",              i),  40,  0,     40);
    hs.add(TString::Format("track%i_injet",         i), TString::Format("track%i in-jet?",                     i),   2,  0,      2);
    hs.add(TString::Format("track%i_inpv",          i), TString::Format("track%i in-PV?",                      i),  10, -1,      9);
    hs.add(TString::Format("track%i_jet_deltaphi0", i), TString::Format("track%i |#Delta#phi| to closest jet", i),  25,  0,      3.15);
  }

  for (int j = 0; j < sv_num_indices; ++j) {
    const char* exc = sv_index_names[j];

    h_sv[j].Init("h_sv_" + std::string(exc), hs, true, do_scatterplots);

    for (int i = 0; i < 4; ++i)
      h_sv_jets_deltaphi[i][j] = fs->make<TH1F>(TString::Format("h_sv_%s_%sjets_deltaphi", exc, lmt_ex[i]), TString::Format(";%s SV #Delta#phi to %sjets;arb. units", exc, lmt_ex[i]), 50, -3.15, 3.15);

    h_sv_bs2derr_bsbs2ddist[j] = fs->make<TH2F>(TString::Format("h_sv_%s_bs2derr_bsbs2ddist", exc), TString::Format("%s SV;dist2d(SV, beamspot) (cm);#sigma(dist2d(SV, beamspot)) (cm)", exc), 500, 0, 2.5, 100, 0, 0.05);
    h_pvrho_bsbs2ddist[j] = fs->make<TH2F>(TString::Format("h_pvrho_sv_%s_bsbs2ddist", exc), TString::Format("%s SV;dist2d(SV, beamspot) (cm);dist2d(PV, beamspot)) (cm)", exc), 5000, 0, 2.5, 200, 0, 0.1);

    h_sv_track_weight[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_weight", exc), TString::Format(";%s SV tracks weight;arb. units", exc), 21, 0, 1.05);
    h_sv_track_q[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_q", exc), TString::Format(";%s SV tracks charge;arb. units.", exc), 4, -2, 2);
    h_sv_track_pt[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_pt", exc), TString::Format(";%s SV tracks p_{T} (GeV);arb. units", exc), 200, 0, 200);
    h_sv_track_eta[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_eta", exc), TString::Format(";%s SV tracks #eta;arb. units", exc), 50, -4, 4);
    h_sv_track_phi[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_phi", exc), TString::Format(";%s SV tracks #phi;arb. units", exc), 50, -3.15, 3.15);
    h_sv_track_dxy[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_dxy", exc), TString::Format(";%s SV tracks dxy (cm);arb. units", exc), 100, 0, 1);
    h_sv_track_dz[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_dz", exc), TString::Format(";%s SV tracks dz (cm);arb. units", exc), 100, 0, 1);
    h_sv_track_pt_err[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_pt_err", exc), TString::Format(";%s SV tracks #sigma(p_{T})/p_{T};arb. units", exc), 200, 0, 2);
    h_sv_track_eta_err[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_eta_err", exc), TString::Format(";%s SV tracks #sigma(#eta);arb. units", exc), 200, 0, 0.02);
    h_sv_track_phi_err[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_phi_err", exc), TString::Format(";%s SV tracks #sigma(#phi);arb. units", exc), 200, 0, 0.02);
    h_sv_track_dxy_err[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_dxy_err", exc), TString::Format(";%s SV tracks #sigma(dxy) (cm);arb. units", exc), 100, 0, 0.1);
    h_sv_track_dz_err[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_dz_err", exc), TString::Format(";%s SV tracks #sigma(dz) (cm);arb. units", exc), 100, 0, 0.1);
    h_sv_track_nsigmadxy[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_nsigmadxy", exc), TString::Format(";%s SV tracks n#sigma(dxy);arb. units", exc), 400, 0, 40);
    h_sv_track_chi2dof[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_chi2dof", exc), TString::Format(";%s SV tracks #chi^{2}/dof;arb. units", exc), 100, 0, 10);
    h_sv_track_npxhits[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_npxhits", exc), TString::Format(";%s SV tracks number of pixel hits;arb. units", exc), 12, 0, 12);
    h_sv_track_nsthits[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_nsthits", exc), TString::Format(";%s SV tracks number of strip hits;arb. units", exc), 28, 0, 28);
    h_sv_track_nhitsbehind[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_nhitsbehind", exc), TString::Format(";%s SV tracks number of hits behind;arb. units", exc), 10, 0, 10);
    h_sv_track_nhitslost[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_nhitslost", exc), TString::Format(";%s SV tracks number of hits lost;arb. units", exc), 10, 0, 10);
    h_sv_track_nhits[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_nhits", exc), TString::Format(";%s SV tracks number of hits", exc), 40, 0, 40);
    h_sv_track_injet[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_injet", exc), TString::Format(";%s SV tracks in-jet?", exc), 2, 0, 2);
    h_sv_track_inpv[j] = fs->make<TH1F>(TString::Format("h_sv_%s_track_inpv", exc), TString::Format(";%s SV tracks in-PV?", exc), 10, -1, 9);
  }

  h_sv_xy = fs->make<TH2F>("h_sv_xy", ";SV x (cm);SV y (cm)", 100, -4, 4, 100, -4, 4);
  h_sv_xz = fs->make<TH2F>("h_sv_xz", ";SV x (cm);SV z (cm)", 100, -4, 4, 100, -25, 25);
  h_sv_yz = fs->make<TH2F>("h_sv_yz", ";SV y (cm);SV z (cm)", 100, -4, 4, 100, -25, 25);
  h_sv_rz = fs->make<TH2F>("h_sv_rz", ";SV r (cm);SV z (cm)", 100, -4, 4, 100, -25, 25);
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
  h_absdeltaphi01_shared_jets = fs->make<TH1F>("h_absdeltaphi01_shared_jets", ";abs(delta(phi of sv #0, phi of sv #1));arb. units", 316, 0, 3.16);
  h_absdeltaphi01_no_shared_jets = fs->make<TH1F>("h_absdeltaphi01_no_shared_jets", ";abs(delta(phi of sv #0, phi of sv #1));arb. units", 316, 0, 3.16);
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
  if (0.00150 < mevent->lspdist2d() < 2) {
	  //////////////////////////////////////////////////////////////////////

	  if (nsv >= 2) {
		  const MFVVertexAux& sv0 = auxes->at(0);
		  const MFVVertexAux& sv1 = auxes->at(1);
		  if (sv0.gen3ddist < 0.0085 || sv1.gen3ddist < 0.0085) {
			  double svdist2d = mag(sv0.x - sv1.x, sv0.y - sv1.y);
			  double svdist3d = mag(sv0.x - sv1.x, sv0.y - sv1.y, sv0.z - sv1.z);
			  h_svdist2d->Fill(svdist2d, w);
			  h_svdist3d->Fill(svdist3d, w);
			  h_sv0pvdz_v_sv1pvdz->Fill(sv0.pvdz(), sv1.pvdz(), w);
			  h_sv0pvdzsig_v_sv1pvdzsig->Fill(sv0.pvdzsig(), sv1.pvdzsig(), w);
			  double phi0 = atan2(sv0.y - bsy, sv0.x - bsx);
			  double phi1 = atan2(sv1.y - bsy, sv1.x - bsx);
			  h_absdeltaphi01->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);

			  h_fractrackssharedwpv01->Fill(float(sv0.ntrackssharedwpv() + sv1.ntrackssharedwpv()) / (sv0.ntracks() + sv1.ntracks()), w);
			  h_fractrackssharedwpvs01->Fill(float(sv0.ntrackssharedwpvs() + sv1.ntrackssharedwpvs()) / (sv0.ntracks() + sv1.ntracks()), w);
			  h_pvmosttracksshared->Fill(sv0.ntrackssharedwpvs() ? sv0.pvmosttracksshared() : -1,
				  sv1.ntrackssharedwpvs() ? sv1.pvmosttracksshared() : -1,
				  w);

			  std::vector<std::vector<int> > sv_track_which_jet;
			  for (int isv = 0; isv < nsv; ++isv) {
				  const MFVVertexAux& aux = auxes->at(isv);
				  const int ntracks = aux.ntracks();

				  std::vector<int> track_which_jet;
				  for (int i = 0; i < ntracks; ++i) {
					  double match_threshold = 1.3;
					  int jet_index = 255;
					  for (unsigned j = 0; j < mevent->jet_track_which_jet.size(); ++j) {
						  double a = fabs(aux.track_pt(i) - fabs(mevent->jet_track_qpt[j])) + 1;
						  double b = fabs(aux.track_eta[i] - mevent->jet_track_eta[j]) + 1;
						  double c = fabs(aux.track_phi[i] - mevent->jet_track_phi[j]) + 1;
						  if (a * b * c < match_threshold) {
							  match_threshold = a * b * c;
							  jet_index = mevent->jet_track_which_jet[j];
						  }
					  }
					  if (jet_index != 255) {
						  track_which_jet.push_back((int)jet_index);
					  }
				  }
				  sv_track_which_jet.push_back(track_which_jet);
			  }

			  bool shared_jet = std::find_first_of(sv_track_which_jet[0].begin(), sv_track_which_jet[0].end(), sv_track_which_jet[1].begin(), sv_track_which_jet[1].end()) != sv_track_which_jet[0].end();
			  h_sv_shared_jets->Fill(shared_jet, w);
			  if (shared_jet) {
				  h_svdist2d_shared_jets->Fill(svdist2d, w);
				  h_absdeltaphi01_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);
			  }
			  else {
				  h_svdist2d_no_shared_jets->Fill(svdist2d, w);
				  h_absdeltaphi01_no_shared_jets->Fill(fabs(reco::deltaPhi(phi0, phi1)), w);
			  }
		  }
	  }
  }
}

DEFINE_FWK_MODULE(MFVVertexHistos);
