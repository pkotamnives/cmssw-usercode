#include "TH2.h"
#include "TMath.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "JMTucker/MFVNeutralinoFormats/interface/VertexerPairEff.h"
#include "JMTucker/MFVNeutralino/interface/VertexerParams.h"
#include "JMTucker/Tools/interface/Utilities.h"

class MFVVertexer : public edm::EDProducer {
public:
  MFVVertexer(const edm::ParameterSet&);
  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  typedef std::set<reco::TrackRef> track_set;
  typedef std::vector<reco::TrackRef> track_vec;

  void finish(edm::Event&, const std::vector<reco::TransientTrack>&, std::unique_ptr<reco::VertexCollection>, std::unique_ptr<VertexerPairEffs>, const std::vector<std::pair<track_set, track_set>>&);

  template <typename T>
  void print_track_set(const T& ts) const {
    for (auto r : ts)
      printf(" %u", r.key());
  }

  template <typename T>
  void print_track_set(const T& ts, const reco::Vertex& v) const {
    for (auto r : ts)
      printf(" %u%s", r.key(), (v.trackWeight(r) < mfv::track_vertex_weight_min ? "!" : ""));
  }

  void print_track_set(const reco::Vertex& v) const {
    for (auto r = v.tracks_begin(), re = v.tracks_end(); r != re; ++r)
      printf(" %lu%s", r->key(), (v.trackWeight(*r) < mfv::track_vertex_weight_min ? "!" : ""));
  }

  bool is_track_subset(const track_set& a, const track_set& b) const {
    bool is_subset = true;
    const track_set& smaller = a.size() <= b.size() ? a : b;
    const track_set& bigger  = a.size() <= b.size() ? b : a;
    
    for (auto t : smaller)
      if (bigger.count(t) < 1) {
        is_subset = false;
        break;
      }

    return is_subset;
  }

  track_set vertex_track_set(const reco::Vertex& v, const double min_weight = mfv::track_vertex_weight_min) const {
    track_set result;

    for (auto it = v.tracks_begin(), ite = v.tracks_end(); it != ite; ++it) {
      const double w = v.trackWeight(*it);
      const bool use = w >= min_weight;
      assert(use);
      //if (verbose) ("trk #%2i pt %6.3f eta %6.3f phi %6.3f dxy %6.3f dz %6.3f w %5.3f  use? %i\n", int(it-v.tracks_begin()), (*it)->pt(), (*it)->eta(), (*it)->phi(), (*it)->dxy(), (*it)->dz(), w, use);
      if (use)
        result.insert(it->castTo<reco::TrackRef>());
    }

    return result;
  }

  track_vec vertex_track_vec(const reco::Vertex& v, const double min_weight = mfv::track_vertex_weight_min) const {
    track_set s = vertex_track_set(v, min_weight);
    return track_vec(s.begin(), s.end());
  }

  Measurement1D vertex_dist(const reco::Vertex& v0, const reco::Vertex& v1) const {
    if (use_2d_vertex_dist)
      return vertex_dist_2d.distance(v0, v1);
    else
      return vertex_dist_3d.distance(v0, v1);
  }

  std::pair<bool, Measurement1D> track_dist(const reco::TransientTrack& t, const reco::Vertex& v) const {
    if (use_2d_track_dist)
      return IPTools::absoluteTransverseImpactParameter(t, v);
    else
      return IPTools::absoluteImpactParameter3D(t, v);
  }

  VertexDistanceXY vertex_dist_2d;
  VertexDistance3D vertex_dist_3d;
  std::unique_ptr<KalmanVertexFitter> kv_reco;

  std::vector<TransientVertex> kv_reco_dropin(std::vector<reco::TransientTrack>& ttks) {
    if (ttks.size() < 2)
      return std::vector<TransientVertex>();
    std::vector<TransientVertex> v(1, kv_reco->vertex(ttks));
    if (v[0].normalisedChiSquared() > 5)
      return std::vector<TransientVertex>();
    return v;
  }

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_token;
  const edm::EDGetTokenT<std::vector<reco::TrackRef>> seed_tracks_token;
  const int n_tracks_per_seed_vertex;
  const double max_seed_vertex_chi2;
  const bool use_2d_vertex_dist;
  const bool use_2d_track_dist;
  const double merge_anyway_dist;
  const double merge_anyway_sig;
  const double merge_shared_dist;
  const double merge_shared_sig;
  const double max_track_vertex_dist;
  const double max_track_vertex_sig;
  const double min_track_vertex_sig_to_remove;
  const bool remove_one_track_at_a_time;
  const double max_nm1_refit_dist3;
  const double max_nm1_refit_distz;
  const int max_nm1_refit_count;
  const bool histos;
  const bool verbose;
  const std::string module_label;

  TH1F* h_n_seed_vertices;
  TH1F* h_seed_vertex_track_weights;
  TH1F* h_seed_vertex_chi2;
  TH1F* h_seed_vertex_ndof;
  TH1F* h_seed_vertex_x;
  TH1F* h_seed_vertex_y;
  TH1F* h_seed_vertex_rho;
  TH1F* h_seed_vertex_phi;
  TH1F* h_seed_vertex_z;
  TH1F* h_seed_vertex_r;
  TH1F* h_seed_vertex_paird2d;
  TH1F* h_seed_vertex_pairdphi;
  TH1F* h_n_resets;
  TH1F* h_n_onetracks;
  TH1F* h_noshare_vertex_tkvtxdist;
  TH1F* h_noshare_vertex_tkvtxdisterr;
  TH1F* h_noshare_vertex_tkvtxdistsig;
  TH1F* h_n_noshare_vertices;
  TH1F* h_noshare_vertex_ntracks;
  TH1F* h_noshare_vertex_track_weights;
  TH1F* h_noshare_vertex_chi2;
  TH1F* h_noshare_vertex_ndof;
  TH1F* h_noshare_vertex_x;
  TH1F* h_noshare_vertex_y;
  TH1F* h_noshare_vertex_rho;
  TH1F* h_noshare_vertex_phi;
  TH1F* h_noshare_vertex_z;
  TH1F* h_noshare_vertex_r;
  TH1F* h_noshare_vertex_paird2d;
  TH1F* h_noshare_vertex_pairdphi;
  TH1F* h_noshare_track_multiplicity;
  TH1F* h_max_noshare_track_multiplicity;
  TH1F* h_n_output_vertices;

  TH2F* h_2D_close_dvv_its_significance_before_merge;
  TH2F* h_2D_close_dvv_its_significance_passed_merge_pairs;
  TH2F* h_2D_close_dvv_its_significance_failed_merge_pairs;
  TH2F* h_2D_close_dvv_its_significance_after_merge;
  TH1F* h_merged_vertex_chi2;
  TH1F* h_non_merged_vertex_chi2;
  TH1F* h_merged_vertex_ntracks;
  TH1F* h_non_merged_vertex_ntracks;
  TH1F* h_merged_vertex_tkvtxdist;
  TH1F* h_non_merged_vertex_tkvtxdist;
  TH1F* h_merged_vertex_tkvtxdistsig;
  TH1F* h_non_merged_vertex_tkvtxdistsig;
  TH1F* h_merged_vertex_mass;
  TH1F* h_non_merged_vertex_mass;
  TH1F* h_merged_vertex_dBV;
  TH1F* h_non_merged_vertex_dBV;
  TH1F* h_merged_vertex_bs2derr;
  TH1F* h_non_merged_vertex_bs2derr;
};

MFVVertexer::MFVVertexer(const edm::ParameterSet& cfg)
  : kv_reco(new KalmanVertexFitter(cfg.getParameter<edm::ParameterSet>("kvr_params"), cfg.getParameter<edm::ParameterSet>("kvr_params").getParameter<bool>("doSmoothing"))),
    beamspot_token(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot_src"))),
    seed_tracks_token(consumes<std::vector<reco::TrackRef>>(cfg.getParameter<edm::InputTag>("seed_tracks_src"))),
    n_tracks_per_seed_vertex(cfg.getParameter<int>("n_tracks_per_seed_vertex")),
    max_seed_vertex_chi2(cfg.getParameter<double>("max_seed_vertex_chi2")),
    use_2d_vertex_dist(cfg.getParameter<bool>("use_2d_vertex_dist")),
    use_2d_track_dist(cfg.getParameter<bool>("use_2d_track_dist")),
    merge_anyway_dist(cfg.getParameter<double>("merge_anyway_dist")),
    merge_anyway_sig(cfg.getParameter<double>("merge_anyway_sig")),
    merge_shared_dist(cfg.getParameter<double>("merge_shared_dist")),
    merge_shared_sig(cfg.getParameter<double>("merge_shared_sig")),
    max_track_vertex_dist(cfg.getParameter<double>("max_track_vertex_dist")),
    max_track_vertex_sig(cfg.getParameter<double>("max_track_vertex_sig")),
    min_track_vertex_sig_to_remove(cfg.getParameter<double>("min_track_vertex_sig_to_remove")),
    remove_one_track_at_a_time(cfg.getParameter<bool>("remove_one_track_at_a_time")),
    max_nm1_refit_dist3(cfg.getParameter<double>("max_nm1_refit_dist3")),
    max_nm1_refit_distz(cfg.getParameter<double>("max_nm1_refit_distz")),
    max_nm1_refit_count(cfg.getParameter<int>("max_nm1_refit_count")),
    histos(cfg.getUntrackedParameter<bool>("histos", false)),
    verbose(cfg.getUntrackedParameter<bool>("verbose", false)),
    module_label(cfg.getParameter<std::string>("@module_label"))
{
  if (n_tracks_per_seed_vertex < 2 || n_tracks_per_seed_vertex > 5)
    throw cms::Exception("MFVVertexer", "n_tracks_per_seed_vertex must be one of 2,3,4,5");

  produces<reco::VertexCollection>();
  produces<VertexerPairEffs>();
  produces<reco::TrackCollection>("seed"); // JMTBAD remove me
  produces<reco::TrackCollection>("inVertices");

  if (histos) {
    edm::Service<TFileService> fs;

    h_n_seed_vertices                = fs->make<TH1F>("h_n_seed_vertices",                "",  50,   0,    200);
    h_seed_vertex_track_weights      = fs->make<TH1F>("h_seed_vertex_track_weights",      "",  21,   0,      1.05);
    h_seed_vertex_chi2               = fs->make<TH1F>("h_seed_vertex_chi2",               "",  20,   0, max_seed_vertex_chi2);
    h_seed_vertex_ndof               = fs->make<TH1F>("h_seed_vertex_ndof",               "",  10,   0,     20);
    h_seed_vertex_x                  = fs->make<TH1F>("h_seed_vertex_x",                  "", 100,  -1,      1);
    h_seed_vertex_y                  = fs->make<TH1F>("h_seed_vertex_y",                  "", 100,  -1,      1);
    h_seed_vertex_rho                = fs->make<TH1F>("h_seed_vertex_rho",                "", 100,   0,      2);
    h_seed_vertex_phi                = fs->make<TH1F>("h_seed_vertex_phi",                "",  50,  -3.15,   3.15);
    h_seed_vertex_z                  = fs->make<TH1F>("h_seed_vertex_z",                  "",  40, -20,     20);
    h_seed_vertex_r                  = fs->make<TH1F>("h_seed_vertex_r",                  "", 100,   0,      2);
    h_seed_vertex_paird2d            = fs->make<TH1F>("h_seed_vertex_paird2d",            "", 100,   0,      0.2);
    h_seed_vertex_pairdphi           = fs->make<TH1F>("h_seed_vertex_pairdphi",           "", 100,  -3.14,   3.14);

    h_n_resets                       = fs->make<TH1F>("h_n_resets",                       "", 50,   0,   500);
    h_n_onetracks                    = fs->make<TH1F>("h_n_onetracks",                    "",  5,   0,     5);

    h_n_noshare_vertices             = fs->make<TH1F>("h_n_noshare_vertices",             "", 50,   0,    50);
    h_noshare_vertex_tkvtxdist       = fs->make<TH1F>("h_noshare_vertex_tkvtxdist",       "", 100,  0,   0.1);
    h_noshare_vertex_tkvtxdisterr    = fs->make<TH1F>("h_noshare_vertex_tkvtxdisterr",    "", 100,  0,   0.1);
    h_noshare_vertex_tkvtxdistsig    = fs->make<TH1F>("h_noshare_vertex_tkvtxdistsig",    "", 100,  0,     6);
    h_noshare_vertex_ntracks         = fs->make<TH1F>("h_noshare_vertex_ntracks",         "",  30,  0, 30);
    h_noshare_vertex_track_weights   = fs->make<TH1F>("h_noshare_vertex_track_weights",   "",  21,   0,      1.05);
    h_noshare_vertex_chi2            = fs->make<TH1F>("h_noshare_vertex_chi2",            "", 20,   0, max_seed_vertex_chi2);
    h_noshare_vertex_ndof            = fs->make<TH1F>("h_noshare_vertex_ndof",            "", 10,   0,     20);
    h_noshare_vertex_x               = fs->make<TH1F>("h_noshare_vertex_x",               "", 100,  -1,      1);
    h_noshare_vertex_y               = fs->make<TH1F>("h_noshare_vertex_y",               "", 100,  -1,      1);
    h_noshare_vertex_rho             = fs->make<TH1F>("h_noshare_vertex_rho",             "", 100,   0,      2);
    h_noshare_vertex_phi             = fs->make<TH1F>("h_noshare_vertex_phi",             "", 50,  -3.15,   3.15);
    h_noshare_vertex_z               = fs->make<TH1F>("h_noshare_vertex_z",               "", 40, -20,     20);
    h_noshare_vertex_r               = fs->make<TH1F>("h_noshare_vertex_r",               "", 100,   0,      2);
    h_noshare_vertex_paird2d         = fs->make<TH1F>("h_noshare_vertex_paird2d",            "", 100,   0,      0.2);
    h_noshare_vertex_pairdphi        = fs->make<TH1F>("h_noshare_vertex_pairdphi",           "", 100,  -3.15,   3.15);
    h_noshare_track_multiplicity     = fs->make<TH1F>("h_noshare_track_multiplicity",     "",  40,   0,     40);
    h_max_noshare_track_multiplicity = fs->make<TH1F>("h_max_noshare_track_multiplicity", "",  40,   0,     40);
    h_n_output_vertices           = fs->make<TH1F>("h_n_output_vertices",           "", 50, 0, 50);

	h_2D_close_dvv_its_significance_before_merge = fs->make<TH2F>("h_2D_close_dvv_its_significance_before_merge", "Before merging by significance<4: dPhi(SV0,SV1)<0.5; svdist3d (cm); svdist3d significance(cm)", 50, 0, 0.1, 100, 0, 30);
	h_2D_close_dvv_its_significance_passed_merge_pairs = fs->make<TH2F>("h_2D_close_dvv_its_significance_passed_merge_pairs", "Only passed merging pairs by significance<4: dPhi(SV0,SV1)<0.5; svdist3d (cm); svdist3d significance(cm)", 50, 0, 0.1, 100, 0, 30);
	h_2D_close_dvv_its_significance_failed_merge_pairs = fs->make<TH2F>("h_2D_close_dvv_its_significance_failed_merge_pairs", "Only failed merging pairs by significance<4: dPhi(SV0,SV1)<0.5; svdist3d (cm); svdist3d significance(cm)", 50, 0, 0.1, 100, 0, 30);
	h_2D_close_dvv_its_significance_after_merge = fs->make<TH2F>("h_2D_close_dvv_its_significance_after_merge", "After merging by significance<4: dPhi(SV0,SV1)<0.5; svdist3d (cm); svdist3d significance(cm) ", 50, 0, 0.1, 100, 0, 30);
	h_merged_vertex_chi2 = fs->make<TH1F>("h_merged_vertex_chi2", "After merging by sigma<4: merged vertices; chi2/dof ", 20, 0, max_seed_vertex_chi2);
	h_non_merged_vertex_chi2 = fs->make<TH1F>("h_non_merged_vertex_chi2", "After merging by sigma<4: non-merged vertices; chi2/dof ", 20, 0, max_seed_vertex_chi2);
	h_merged_vertex_ntracks = fs->make<TH1F>("h_merged_vertex_ntracks", "After merging by sigma<4: merged vertices; # of tracks/vtx ", 30, 0, 30);
	h_non_merged_vertex_ntracks = fs->make<TH1F>("h_non_merged_vertex_ntracks", "After merging by sigma<4: non-merged vertices; # of tracks/vtx ", 30, 0, 30);
	h_merged_vertex_tkvtxdist = fs->make<TH1F>("h_merged_vertex_tkvtxdist", "After merging by sigma<4: merged vertices; miss dist 3d (cm)", 50, 0, 0.1);
	h_non_merged_vertex_tkvtxdist = fs->make<TH1F>("h_non_merged_vertex_tkvtxdist", "After merging by sigma<4: non-merged vertices; miss dist 3d (cm)", 50, 0, 0.1);
	h_merged_vertex_tkvtxdistsig = fs->make<TH1F>("h_merged_vertex_tkvtxdistsig", "After merging by sigma<4: merged vertices; miss dist significance 3d ", 50, 0, 10);
	h_non_merged_vertex_tkvtxdistsig = fs->make<TH1F>("h_non_merged_vertex_tkvtxdistsig", "After merging by sigma<4: non-merged vertices; miss dist significance 3d ", 50, 0, 10);
	h_merged_vertex_mass = fs->make<TH1F>("h_merged_vertex_mass", "After merging by sigma<4: merged vertices; vtx mass (GeV)", 50, 0, 2000);
	h_non_merged_vertex_mass = fs->make<TH1F>("h_non_merged_vertex_mass", "After merging by sigma<4: non-merged vertices; vtx mass (GeV)", 50, 0, 2000);
	h_merged_vertex_dBV = fs->make<TH1F>("h_merged_vertex_dBV", "After merging by sigma<4: merged vertices; bvdist2d (cm)", 50, 0, 0.1);
	h_non_merged_vertex_dBV = fs->make<TH1F>("h_non_merged_vertex_dBV", "After merging by sigma<4: non-merged vertices; bvdist2d (cm)", 50, 0, 0.1);
	h_merged_vertex_bs2derr = fs->make<TH1F>("h_merged_vertex_bs2derr", "After merging by sigma<4: merged vertices; bvdist2d error (cm)", 10, 0, 0.05);
	h_non_merged_vertex_bs2derr = fs->make<TH1F>("h_non_merged_vertex_bs2derr", "After merging by sigma<4: non-merged vertices; bvdist2d error (cm)", 10, 0, 0.05);

  }
}

void MFVVertexer::finish(edm::Event& event, const std::vector<reco::TransientTrack>& seed_tracks, std::unique_ptr<reco::VertexCollection> vertices, std::unique_ptr<VertexerPairEffs> vpeffs, const std::vector<std::pair<track_set, track_set>>& vpeffs_tracks) {
  std::unique_ptr<reco::TrackCollection> tracks_seed      (new reco::TrackCollection);
  std::unique_ptr<reco::TrackCollection> tracks_inVertices(new reco::TrackCollection);

  if (verbose) printf("finish:\nseed tracks:\n");

  std::map<std::pair<unsigned, unsigned>, unsigned char> seed_track_ref_map;
  unsigned char itk = 0;
  for (const reco::TransientTrack& ttk : seed_tracks) {
    tracks_seed->push_back(ttk.track());
    const reco::TrackBaseRef& tk(ttk.trackBaseRef());
    seed_track_ref_map[std::make_pair(tk.id().id(), tk.key())] = uint2uchar_clamp(itk++);

    if (verbose) printf("id: %i key: %lu pt: %f\n", tk.id().id(), tk.key(), tk->pt());
  }

  assert(vpeffs->size() == vpeffs_tracks.size());
  for (size_t i = 0, ie = vpeffs->size(); i < ie; ++i) {
    for (auto tk : vpeffs_tracks[i].first)  (*vpeffs)[i].tracks_push_back(0, seed_track_ref_map[std::make_pair(tk.id().id(), tk.key())]);
    for (auto tk : vpeffs_tracks[i].second) (*vpeffs)[i].tracks_push_back(1, seed_track_ref_map[std::make_pair(tk.id().id(), tk.key())]);
  }

  if (verbose) printf("vertices:\n");

  for (const reco::Vertex& v : *vertices) {
    if (verbose) printf("x: %f y %f z %f\n", v.x(), v.y(), v.z());
    for (auto it = v.tracks_begin(), ite = v.tracks_end(); it != ite; ++it) {
      reco::TrackRef tk = it->castTo<reco::TrackRef>();
      if (verbose) printf("id: %i key: %u <%f,%f,%f,%f,%f>\n", tk.id().id(), tk.key(), tk->charge()*tk->pt(), tk->eta(), tk->phi(), tk->dxy(), tk->dz());
      tracks_inVertices->push_back(*tk);
    }
  }

  if (verbose)
    printf("n_output_vertices: %lu\n", vertices->size());
  if (histos)
    h_n_output_vertices->Fill(vertices->size());

  event.put(std::move(vertices));
  event.put(std::move(vpeffs));
  event.put(std::move(tracks_seed),       "seed");
  event.put(std::move(tracks_inVertices), "inVertices");
}

void MFVVertexer::produce(edm::Event& event, const edm::EventSetup& setup) {
  if (verbose)
    std::cout << "MFVVertexer " << module_label << " run " << event.id().run() << " lumi " << event.luminosityBlock() << " event " << event.id().event() << "\n";

  edm::Handle<reco::BeamSpot> beamspot;
  event.getByToken(beamspot_token, beamspot);
  const double bsx = beamspot->position().x();
  const double bsy = beamspot->position().y();
  const double bsz = beamspot->position().z();

  edm::ESHandle<TransientTrackBuilder> tt_builder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", tt_builder);

  edm::Handle<std::vector<reco::TrackRef>> seed_track_refs;
  event.getByToken(seed_tracks_token, seed_track_refs);

  std::vector<reco::TransientTrack> seed_tracks;
  std::map<reco::TrackRef, size_t> seed_track_ref_map;

  for (const reco::TrackRef& tk : *seed_track_refs) {
    seed_tracks.push_back(tt_builder->build(tk));
    seed_track_ref_map[tk] = seed_tracks.size() - 1;
  }

  const size_t ntk = seed_tracks.size();
  if (verbose)
    printf("n_seed_tracks: %5lu\n", ntk);

  //////////////////////////////////////////////////////////////////////
  // Form seed vertices from all pairs of tracks whose vertex fit
  // passes cuts.
  //////////////////////////////////////////////////////////////////////

  std::unique_ptr<reco::VertexCollection> vertices(new reco::VertexCollection);
  std::unique_ptr<VertexerPairEffs> vpeffs(new VertexerPairEffs);
  std::vector<std::pair<track_set, track_set>> vpeffs_tracks;

  if (ntk == 0) {
    if (verbose)
      printf("no seed tracks -> putting empty vertex collection into event\n");
    finish(event, seed_tracks, std::move(vertices), std::move(vpeffs), vpeffs_tracks);
    return;
  }

  std::vector<size_t> itks(n_tracks_per_seed_vertex, 0);

  auto try_seed_vertex = [&]() {
    std::vector<reco::TransientTrack> ttks(n_tracks_per_seed_vertex);
    for (int i = 0; i < n_tracks_per_seed_vertex; ++i)
      ttks[i] = seed_tracks[itks[i]];	// pk: how can we have the value of itks[i]?(the loop below) where seed_tracks is a vector of "TransientTrack"

    TransientVertex seed_vertex = kv_reco->vertex(ttks);
    if (seed_vertex.isValid() && seed_vertex.normalisedChiSquared() < max_seed_vertex_chi2) {
      vertices->push_back(reco::Vertex(seed_vertex));	   // pk: collect seed vertices

      if (verbose || histos) {
        const reco::Vertex& v = vertices->back();
        const double vchi2 = v.normalizedChi2();
        const double vndof = v.ndof();
        const double vx = v.position().x() - bsx;
        const double vy = v.position().y() - bsy;
        const double vz = v.position().z() - bsz;
        const double phi = atan2(vy, vx);
        const double rho = mag(vx, vy);
        const double r = mag(vx, vy, vz);
        if (verbose) {
          printf("from tracks");
          for (auto itk : itks)
            printf(" %lu", itk);
          printf(": vertex #%3lu: chi2/dof: %7.3f dof: %7.3f pos: <%7.3f, %7.3f, %7.3f>  rho: %7.3f  phi: %7.3f  r: %7.3f\n", vertices->size()-1, vchi2, vndof, vx, vy, vz, rho, phi, r);
        }
        if (histos) {
          for (auto it = v.tracks_begin(), ite = v.tracks_end(); it != ite; ++it)
            h_seed_vertex_track_weights->Fill(v.trackWeight(*it));
          h_seed_vertex_chi2->Fill(vchi2);
          h_seed_vertex_ndof->Fill(vndof);
          h_seed_vertex_x->Fill(vx);
          h_seed_vertex_y->Fill(vy);
          h_seed_vertex_rho->Fill(rho);
          h_seed_vertex_phi->Fill(phi);
          h_seed_vertex_z->Fill(vz);
          h_seed_vertex_r->Fill(r);
        }
      }
    }
  };

  // ha
  for (size_t itk = 0; itk < ntk; ++itk) {	 // pk: loop trough all possible 2-track or 3-track or 4-track seed vertices in ntk but the order is stil preserved?
    itks[0] = itk;
    for (size_t jtk = itk+1; jtk < ntk; ++jtk) {
      itks[1] = jtk;
      if (n_tracks_per_seed_vertex == 2) { try_seed_vertex(); continue; }
      for (size_t ktk = jtk+1; ktk < ntk; ++ktk) {
        itks[2] = ktk;
        if (n_tracks_per_seed_vertex == 3) { try_seed_vertex(); continue; }
        for (size_t ltk = ktk+1; ltk < ntk; ++ltk) {
          itks[3] = ltk;
          if (n_tracks_per_seed_vertex == 4) { try_seed_vertex(); continue; }
          for (size_t mtk = ltk+1; mtk < ntk; ++mtk) {
            itks[4] = mtk;
            try_seed_vertex();
          }
        }
      }
    }
  }

  if (histos) {
    for (std::vector<reco::Vertex>::const_iterator v0 = vertices->begin(); v0 != vertices->end(); ++v0) {
      const double v0x = v0->position().x() - bsx;
      const double v0y = v0->position().y() - bsy;
      const double phi0 = atan2(v0y, v0x);
      for (std::vector<reco::Vertex>::const_iterator v1 = v0 + 1; v1 != vertices->end(); ++v1) {
        const double v1x = v1->position().x() - bsx;
        const double v1y = v1->position().y() - bsy;
        const double phi1 = atan2(v1y, v1x);
        h_seed_vertex_paird2d ->Fill(mag(v0x - v1x, v0y - v1y));
        h_seed_vertex_pairdphi->Fill(reco::deltaPhi(phi0, phi1));
      }
    }
  }

  if (verbose)
    printf("n_seed_vertices: %lu\n", vertices->size());
  if (histos)
    h_n_seed_vertices->Fill(vertices->size());

  //////////////////////////////////////////////////////////////////////
  // Take care of track sharing. If a track is in two vertices, and
  // the vertices are "close", refit the tracks from the two together
  // as one vertex. If the vertices are not close, keep the track in
  // the vertex to which it is "closer".
  //////////////////////////////////////////////////////////////////////

  if (verbose)
    printf("fun time!\n");

  track_set discarded_tracks;
  int n_resets = 0;
  int n_onetracks = 0;
  std::vector<reco::Vertex>::iterator v[2];	   // pk: v[2] is an iterator?
  std::vector<reco::Vertex>::iterator nv[2];
  size_t ivtx[2];
  for (v[0] = vertices->begin(); v[0] != vertices->end(); ++v[0]) {
    track_set tracks[2];			// pk: tracks[2] isn't defined anywhere	and what is track_set
    ivtx[0] = v[0] - vertices->begin();
    tracks[0] = vertex_track_set(*v[0]);

    if (tracks[0].size() < 2) {
      if (verbose)
        printf("track-sharing: vertex-0 #%lu is down to one track, junking it\n", ivtx[0]);
      v[0] = vertices->erase(v[0]) - 1;
      ++n_onetracks;
      continue;
    }

    bool duplicate = false;
    bool merge = false;
    bool refit = false;
    track_set tracks_to_remove_in_refit[2];
    VertexerPairEff* vpeff = 0;
    const size_t max_vpeffs_size = 20000; // enough for 200 vertices to share tracks

    for (v[1] = v[0] + 1; v[1] != vertices->end(); ++v[1]) {
      ivtx[1] = v[1] - vertices->begin();
      tracks[1] = vertex_track_set(*v[1]);

      if (tracks[1].size() < 2) {
        if (verbose)
          printf("track-sharing: vertex-1 #%lu is down to one track, junking it\n", ivtx[1]);
        v[1] = vertices->erase(v[1]) - 1;
        ++n_onetracks;
        continue;
      }

      if (verbose) {
        printf("track-sharing: # vertices = %lu. considering vertices #%lu (chi2/dof %.3f prob %.2e, track set", vertices->size(), ivtx[0], v[0]->chi2()/v[0]->ndof(), TMath::Prob(v[0]->chi2(), int(v[0]->ndof())));
        print_track_set(tracks[0], *v[0]);
        printf(") and #%lu (chi2/dof %.3f prob %.2e, track set", ivtx[1], v[1]->chi2()/v[1]->ndof(), TMath::Prob(v[1]->chi2(), int(v[1]->ndof())));
        print_track_set(tracks[1], *v[1]);
        printf("):\n");
      }

      if (is_track_subset(tracks[0], tracks[1])) {
        if (verbose)
          printf("   subset/duplicate vertices %lu and %lu, erasing second and starting over\n", ivtx[0], ivtx[1]);
        duplicate = true;
        break;
      }

      if (vpeffs->size() < max_vpeffs_size) {
        std::pair<track_set, track_set> vpeff_tracks(tracks[0], tracks[1]);
        auto it = std::find(vpeffs_tracks.begin(), vpeffs_tracks.end(), vpeff_tracks);
        if (it != vpeffs_tracks.end()) {
          vpeffs->at(it - vpeffs_tracks.begin()).inc_weight();			 // pk: what is vertexerpaireff? why is it zero when equal? like shared jets alg.?
          vpeff = 0;
        }
        else {
          vpeffs->push_back(VertexerPairEff());		
          vpeff = &vpeffs->back();
          vpeff->set_vertices(*v[0], *v[1]);
          vpeffs_tracks.push_back(vpeff_tracks);
        }
      }
      else
        vpeff = 0;

      reco::TrackRefVector shared_tracks;
      for (auto tk : tracks[0])
        if (tracks[1].count(tk) > 0)
          shared_tracks.push_back(tk);

      if (verbose) {
        if (shared_tracks.size()) {
          printf("   shared tracks are: ");
          print_track_set(shared_tracks);
          printf("\n");
        }
        else
          printf("   no shared tracks\n");
      }  

      if (shared_tracks.size() > 0) {
        if (vpeff)
          vpeff->kind(VertexerPairEff::share);

        Measurement1D v_dist = vertex_dist(*v[0], *v[1]);
        if (verbose)
          printf("   vertex dist (2d? %i) %7.3f  sig %7.3f\n", use_2d_vertex_dist, v_dist.value(), v_dist.significance());

        if (v_dist.value() < merge_shared_dist || v_dist.significance() < merge_shared_sig) {
          if (verbose) printf("          dist < %7.3f || sig < %7.3f, will try using merge result first before arbitration\n", merge_shared_dist, merge_shared_sig);
          merge = true;
        }
        else
          refit = true;

        if (verbose) printf("   checking for arbitration refit:\n");
        for (auto tk : shared_tracks) {
          const reco::TransientTrack& ttk = seed_tracks[seed_track_ref_map[tk]];
          std::pair<bool, Measurement1D> t_dist_0 = track_dist(ttk, *v[0]);
          std::pair<bool, Measurement1D> t_dist_1 = track_dist(ttk, *v[1]);
          if (verbose) {
            printf("      track-vertex0 dist (2d? %i) calc success? %i  dist %7.3f  sig %7.3f\n", use_2d_track_dist, t_dist_0.first, t_dist_0.second.value(), t_dist_0.second.significance());
            printf("      track-vertex1 dist (2d? %i) calc success? %i  dist %7.3f  sig %7.3f\n", use_2d_track_dist, t_dist_1.first, t_dist_1.second.value(), t_dist_1.second.significance());
          }

          t_dist_0.first = t_dist_0.first && (t_dist_0.second.value() < max_track_vertex_dist || t_dist_0.second.significance() < max_track_vertex_sig);
          t_dist_1.first = t_dist_1.first && (t_dist_1.second.value() < max_track_vertex_dist || t_dist_1.second.significance() < max_track_vertex_sig);
          bool remove_from_0 = !t_dist_0.first;
          bool remove_from_1 = !t_dist_1.first;
          if (t_dist_0.second.significance() < min_track_vertex_sig_to_remove && t_dist_1.second.significance() < min_track_vertex_sig_to_remove) {
            if (tracks[0].size() > tracks[1].size())
              remove_from_1 = true;
            else
              remove_from_0 = true;
          }
          else if (t_dist_0.second.significance() < t_dist_1.second.significance())
            remove_from_1 = true;
          else
            remove_from_0 = true;

          if (verbose) {
            printf("   for tk %u:\n", tk.key());
            printf("      track-vertex0 dist < %7.3f || sig < %7.3f ? %i  remove? %i\n", max_track_vertex_dist, max_track_vertex_sig, t_dist_0.first, remove_from_0);
            printf("      track-vertex1 dist < %7.3f || sig < %7.3f ? %i  remove? %i\n", max_track_vertex_dist, max_track_vertex_sig, t_dist_1.first, remove_from_1);
          }

          if (remove_from_0) tracks_to_remove_in_refit[0].insert(tk);
          if (remove_from_1) tracks_to_remove_in_refit[1].insert(tk);

          if (remove_one_track_at_a_time) {
            if (verbose)
              printf("   arbitrate only one track at a time\n");
            break;
          }
        }

        if (verbose)
          printf("   breaking to refit\n");
        
        break;
      }

      if (verbose) printf("   moving on to next vertex pair.\n");
    }

    if (duplicate) {
      vertices->erase(v[1]);
    }
    else if (merge) {
      if (verbose)
        printf("      before merge, # total vertices = %lu\n", vertices->size());

      track_set tracks_to_fit;
      for (int i = 0; i < 2; ++i)
        for (auto tk : tracks[i])
          tracks_to_fit.insert(tk);

      if (verbose) {
        printf("   merging vertices %lu and %lu with these tracks:", ivtx[0], ivtx[1]);
        print_track_set(tracks_to_fit);
        printf("\n");
      }

      std::vector<reco::TransientTrack> ttks;
      for (auto tk : tracks_to_fit)
        ttks.push_back(seed_tracks[seed_track_ref_map[tk]]);
      
      reco::VertexCollection new_vertices;
      for (const TransientVertex& tv : kv_reco_dropin(ttks))
        new_vertices.push_back(reco::Vertex(tv));
      
      if (verbose) {
        printf("      got %lu new vertices out of the av fit\n", new_vertices.size());
        printf("      these (chi2/dof : prob | track sets):");
        for (const auto& nv : new_vertices) {
          printf(" (%.3f : %.2e | ", nv.chi2()/nv.ndof(), TMath::Prob(nv.chi2(), int(nv.ndof())));
          print_track_set(nv);
          printf(" ),");
        }
        printf("\n");
      }

      // If we got two new vertices, maybe it took A B and A C D and made a better one from B C D, and left a broken one A B! C! D!.
      // If we get one that is truly the merger of the track lists, great. If it is just something like A B , A C -> A B C!, or we get nothing, then default to arbitration.
      if (new_vertices.size() > 1) {
        if (verbose)
          printf("   jiggled again?\n");   
        assert(new_vertices.size() == 2);
        *v[1] = reco::Vertex(new_vertices[1]);
        *v[0] = reco::Vertex(new_vertices[0]);
      }
      else if (new_vertices.size() == 1 && vertex_track_set(new_vertices[0], 0) == tracks_to_fit) {
        if (verbose)
          printf("   merge worked!\n");   

        if (vpeff)
          vpeff->kind(VertexerPairEff::merge);

        vertices->erase(v[1]);
        *v[0] = reco::Vertex(new_vertices[0]); // ok to use v[0] after the erase(v[1]) because v[0] is by construction before v[1]
      }
      else {
        if (verbose)
          printf("   merge didn't work, trying arbitration refits\n");   
        refit = true;
      }

      if (verbose)
        printf("   vertices size is now %lu\n", vertices->size());
    }

    if (refit) {
      bool erase[2] = { false };
      reco::Vertex vsave[2] = { *v[0], *v[1] };

      for (int i = 0; i < 2; ++i) {
        if (tracks_to_remove_in_refit[i].empty())
          continue;

        if (verbose) {
          printf("   refit vertex%i %lu with these tracks:", i, ivtx[i]);
          print_track_set(tracks[i]);
          printf("   but skip these:");
          print_track_set(tracks_to_remove_in_refit[i]);
          printf("\n");
        }

        std::vector<reco::TransientTrack> ttks;
        for (auto tk : tracks[i])
          if (tracks_to_remove_in_refit[i].count(tk) == 0) 
            ttks.push_back(seed_tracks[seed_track_ref_map[tk]]);

        reco::VertexCollection new_vertices;
        for (const TransientVertex& tv : kv_reco_dropin(ttks))
          new_vertices.push_back(reco::Vertex(tv));
        if (verbose) {
          printf("      got %lu new vertices out of the av fit for v%i\n", new_vertices.size(), i);
          printf("      these track sets:");
          for (const auto& nv : new_vertices) {
            printf(" (");
            print_track_set(nv);
            printf(" ),");
          }
          printf("\n");
        }
        if (new_vertices.size() == 1)
          *v[i] = new_vertices[0];
        else
          erase[i] = true;
      }

      if (vpeff && (erase[0] || erase[1]))
        vpeff->kind(VertexerPairEff::erase);

      if (erase[1]) vertices->erase(v[1]);
      if (erase[0]) vertices->erase(v[0]);

      if (verbose)
        printf("      vertices size is now %lu\n", vertices->size());
    }

    // If we changed the vertices at all, start loop over completely.
    if (duplicate || merge || refit) {
      v[0] = vertices->begin() - 1;  // -1 because about to ++sv
      ++n_resets;
      if (verbose) printf("   resetting from vertices %lu and %lu. # of resets: %i\n", ivtx[0], ivtx[1], n_resets);
      
      //if (n_resets == 3000)
      //  throw "I'm dumb";
    }
  }

  if (verbose)
    printf("n_resets: %i  n_onetracks: %i  n_noshare_vertices: %lu\n", n_resets, n_onetracks, vertices->size());
  if (histos) {
    h_n_resets->Fill(n_resets);
    h_n_onetracks->Fill(n_onetracks);
    h_n_noshare_vertices->Fill(vertices->size());
  }

  if (histos || verbose) {
    std::map<reco::TrackRef, int> track_use;
    for (size_t i = 0, ie = vertices->size(); i < ie; ++i) {
      const reco::Vertex& v = vertices->at(i);
      const int ntracks = v.nTracks();
      const double vchi2 = v.normalizedChi2();
      const double vndof = v.ndof();
      const double vx = v.position().x() - bsx;
      const double vy = v.position().y() - bsy;
      const double vz = v.position().z() - bsz;
      const double rho = mag(vx, vy);
      const double phi = atan2(vy, vx);
      const double r = mag(vx, vy, vz);
      for (const auto& r : vertex_track_set(v)) {
        if (track_use.find(r) != track_use.end())
          track_use[r] += 1;
        else
          track_use[r] = 1;
      }

      if (verbose)
        printf("no-share vertex #%3lu: ntracks: %i chi2/dof: %7.3f dof: %7.3f pos: <%7.3f, %7.3f, %7.3f>  rho: %7.3f  phi: %7.3f  r: %7.3f\n", i, ntracks, vchi2, vndof, vx, vy, vz, rho, phi, r);

      if (histos) {
        h_noshare_vertex_ntracks->Fill(ntracks);
        for (auto it = v.tracks_begin(), ite = v.tracks_end(); it != ite; ++it) {
	  h_noshare_vertex_track_weights->Fill(v.trackWeight(*it));

	  reco::TransientTrack seed_track;
	  seed_track = tt_builder->build(*it.operator*());
	  std::pair<bool, Measurement1D> tk_vtx_dist = track_dist(seed_track, v);
	  h_noshare_vertex_tkvtxdist->Fill(tk_vtx_dist.second.value());
	  h_noshare_vertex_tkvtxdisterr->Fill(tk_vtx_dist.second.error());
	  h_noshare_vertex_tkvtxdistsig->Fill(tk_vtx_dist.second.significance());
	}
        h_noshare_vertex_chi2->Fill(vchi2);
        h_noshare_vertex_ndof->Fill(vndof);
        h_noshare_vertex_x->Fill(vx);
        h_noshare_vertex_y->Fill(vy);
        h_noshare_vertex_rho->Fill(rho);
        h_noshare_vertex_phi->Fill(phi);
        h_noshare_vertex_z->Fill(vz);
        h_noshare_vertex_r->Fill(r);

        for (size_t j = i+1, je = vertices->size(); j < je; ++j) {
          const reco::Vertex& vj = vertices->at(j);
          const double vjx = vj.position().x() - bsx;
          const double vjy = vj.position().y() - bsy;
          const double phij = atan2(vjy, vjx);
          h_noshare_vertex_paird2d->Fill(mag(vx - vjx, vy - vjy));
          h_noshare_vertex_pairdphi->Fill(reco::deltaPhi(phi, phij));
        }
      }
    }
    
    if (verbose)
      printf("track multiple uses:\n");

    int max_noshare_track_multiplicity = 0;
    for (const auto& p : track_use) {
      if (verbose && p.second > 1)
        printf("track %3u used %3i times\n", p.first.key(), p.second);
      if (histos)
        h_noshare_track_multiplicity->Fill(p.second);
      if (p.second > max_noshare_track_multiplicity)
        max_noshare_track_multiplicity = p.second;
    }
    if (histos)
      h_max_noshare_track_multiplicity->Fill(max_noshare_track_multiplicity);
  }

  //////////////////////////////////////////////////////////////////////
  // Merge vertices that are still "close". JMTBAD this doesn't do anything currently, only run in verbose mode
  //////////////////////////////////////////////////////////////////////

  if (verbose)
    printf("fun2! before merge loop, # vertices = %lu\n", vertices->size());
  
  
  std::cout << __LINE__ << std::endl;
  if (merge_anyway_sig > 0 || merge_anyway_dist > 0){
	  double v0x;
	  double v0y;
	  //double v0z;
	  double phi0;
	  std::cout << __LINE__ << std::endl;
    for (v[0] = vertices->begin(); v[0] != vertices->end(); ++v[0]) {
      ivtx[0] = v[0] - vertices->begin();
	  std::cout << __LINE__ << std::endl;
	  double v1x;
	  double v1y;
	  //double v1z;
	  double phi1;

      bool merge = false;
      for (v[1] = v[0] + 1; v[1] != vertices->end(); ++v[1]) {
		std::cout << __LINE__ << std::endl;
		std::cout << vertices.size() << std::endl;
        ivtx[1] = v[1] - vertices->begin();
		
		std::cout << __LINE__ << std::endl;
		std::cout << v[1] - vertices->begin() << "in" << vertices.size() << std::endl;
        if (verbose)
          printf("close-merge: # vertices = %lu. considering vertices #%lu (ntk = %i) and #%lu (ntk = %i):", vertices->size(), ivtx[0], v[0]->nTracks(), ivtx[1], v[1]->nTracks());

        Measurement1D v_dist = vertex_dist(*v[0], *v[1]);
        if (verbose)
          printf("   vertex dist (2d? %i) %7.3f  sig %7.3f\n", use_2d_vertex_dist, v_dist.value(), v_dist.significance());
		
		v0x = v[0]->x() - bsx;
		v0y = v[0]->y() - bsy;
		//v0z = v[0]->z() - bsz;
		phi0 = atan2(v0y, v0x);
		v1x = v[1]->x() - bsx;
		v1y = v[1]->y() - bsy;
		//v1z = v[1]->z() - bsz;
		phi1 = atan2(v1y, v1x);

		if (reco::deltaPhi(phi0, phi1) < 0.5)
			h_2D_close_dvv_its_significance_before_merge->Fill(v_dist.value(), v_dist.significance());

        if (v_dist.value() < merge_anyway_dist || v_dist.significance() < merge_anyway_sig) {
          if (verbose)
            printf("          dist < %7.3f || sig < %7.3f, breaking to merge\n", merge_anyway_dist, merge_anyway_sig);
          merge = true;
		  
		  std::cout << __LINE__ << std::endl;
		  
			  std::vector<reco::TransientTrack> ttks;
			  
			  for (int i = 0; i < 2; ++i) {
				  std::cout << __LINE__ << std::endl;
				  std::cout << i << ": " << vertex_track_set(*v[i]).size() << std::endl;
				  std::cout << v[i] - vertices->begin() << "in" << vertices.size() << std::endl;
				  std::cout << __LINE__ << std::endl;
				  for (auto tk : vertex_track_set(*v[i])) {
					  
					  ttks.push_back(tt_builder->build(tk));
					  
				  }

			  }

			  std::cout << __LINE__ << std::endl;
			  reco::VertexCollection new_vertices;
			  for (const TransientVertex& tv : kv_reco_dropin(ttks))
			  {
				  new_vertices.push_back(reco::Vertex(tv));
				  h_merged_vertex_chi2->Fill(double(new_vertices[0].normalizedChi2()));
				  h_merged_vertex_ntracks->Fill(double(new_vertices[0].nTracks()));
				  h_merged_vertex_mass->Fill(double(new_vertices[0].p4().mass()));

				  const reco::Vertex fake_bs_vtx(beamspot->position(), beamspot->covariance3D());
				  Measurement1D dBV_Meas1D = vertex_dist_2d.distance(new_vertices[0], fake_bs_vtx); // where vtx is your reco::Vertex, which maybe means *v[0] but I don't remember offhand. make sure you use the 2D distance here, since that's what we actually use for dBV!!
				  double dBV = dBV_Meas1D.value();
				  double bs2derr = dBV_Meas1D.error();

				  h_merged_vertex_dBV->Fill(dBV);
				  h_merged_vertex_bs2derr->Fill(bs2derr);
				  std::cout << __LINE__ << std::endl;
				  for (auto it = new_vertices[0].tracks_begin(), ite = new_vertices[0].tracks_end(); it != ite; ++it) {

					  reco::TransientTrack seed_track;
					  seed_track = tt_builder->build(*it.operator*());
					  std::pair<bool, Measurement1D> tk_vtx_dist = track_dist(seed_track, new_vertices[0]);
					  h_merged_vertex_tkvtxdist->Fill(tk_vtx_dist.second.value());
					  h_merged_vertex_tkvtxdistsig->Fill(tk_vtx_dist.second.significance());
				  }

				  std::cout << __LINE__ << std::endl;
			  }


			  if (verbose) {
				  printf("      got %lu new vertices out of the av fit\n", new_vertices.size());
				  printf("      these track sets:");
				  for (const auto& nv : new_vertices) {
					  printf(" (");
					  print_track_set(nv);
					  printf(" ),");
				  }
				  printf("\n");
			  }
			  // pk: change vertices
			  std::cout << __LINE__ << std::endl;
			  if (new_vertices.size() == 1)
			  {
				  if (reco::deltaPhi(phi0, phi1) < 0.5) {
					  h_2D_close_dvv_its_significance_passed_merge_pairs->Fill(v_dist.value(), v_dist.significance());
				  }

				  std::cout << "check no mem out of ranges (before) : " << v[1] - vertices->begin() << std::endl;
				  *v[0] = new_vertices[0];
				  std::cout << "check no mem out of ranges (after) : " << v[1] - vertices->begin() << std::endl;
				  std::cout << __LINE__ << std::endl;
				  std::cout << vertices.size() << std::endl;
				  vertices->erase(v[1]);
				  std::cout << vertices.size() << std::endl;
				  std::cout << __LINE__ << std::endl;
			  }
			  // pk: change vertices
			  else {	// pk: size is zero 
				  if (reco::deltaPhi(phi0, phi1) < 0.5)
					  h_2D_close_dvv_its_significance_failed_merge_pairs->Fill(v_dist.value(), v_dist.significance());
			  }
			  std::cout << __LINE__ << std::endl;

		  
        }
      }

	  

      if (!merge) { //until the last one in vertices 
		  h_non_merged_vertex_chi2->Fill(double(v[0]->normalizedChi2()));
		  h_non_merged_vertex_ntracks->Fill(double(v[0]->nTracks()));
		  h_non_merged_vertex_mass->Fill(double(v[0]->p4().mass()));
		  

		  const reco::Vertex fake_bs_vtx(beamspot->position(), beamspot->covariance3D());
		  Measurement1D dBV_Meas1D = vertex_dist_2d.distance(*v[0], fake_bs_vtx); // where vtx is your reco::Vertex, which maybe means *v[0] but I don't remember offhand. make sure you use the 2D distance here, since that's what we actually use for dBV!!
		  double dBV = dBV_Meas1D.value();
		  double bs2derr = dBV_Meas1D.error();

		  h_non_merged_vertex_dBV->Fill(dBV);
		  h_non_merged_vertex_bs2derr->Fill(bs2derr);

		  
		  for (auto it = v[0]->tracks_begin(), ite = v[0]->tracks_end(); it != ite; ++it) {
			  
			  reco::TransientTrack seed_track;
			  seed_track = tt_builder->build(*it.operator*());
			  std::pair<bool, Measurement1D> tk_vtx_dist = track_dist(seed_track, *v[0]);
			  
			  h_non_merged_vertex_tkvtxdist->Fill(tk_vtx_dist.second.value());
			  h_non_merged_vertex_tkvtxdistsig->Fill(tk_vtx_dist.second.significance());
			  
		  }
	  }
	  
    }

	double nv0x;
	double nv0y;
	//double nv0z;
	double nvphi0;
	
	for (nv[0] = vertices->begin(); nv[0] != vertices->end(); ++nv[0]) {
		
		double nv1x;
		double nv1y;
		//double nv1z;
		double nvphi1;

		for (nv[1] = nv[0] + 1; nv[1] != vertices->end(); ++nv[1]) {

			
			Measurement1D nv_dist = vertex_dist(*nv[0], *nv[1]);
			if (verbose)
				printf("  new vertex dist (2d? %i) %7.3f  sig %7.3f\n", use_2d_vertex_dist, nv_dist.value(), nv_dist.significance());

			nv0x = nv[0]->x() - bsx;
			nv0y = nv[0]->y() - bsy;
			//nv0z = nv[0]->z() - bsz;
			nvphi0 = atan2(nv0y, nv0x);
			nv1x = nv[1]->x() - bsx;
			nv1y = nv[1]->y() - bsy;
			//nv1z = nv[1]->z() - bsz;
			nvphi1 = atan2(nv1y, nv1x);
			
			if (reco::deltaPhi(nvphi0, nvphi1) < 0.5)
				h_2D_close_dvv_its_significance_after_merge->Fill(nv_dist.value(), nv_dist.significance());
		}
		
	}
  }

  std::cout << __LINE__ << std::endl;
  //////////////////////////////////////////////////////////////////////
  // Drop tracks that "move" the vertex too much by refitting without each track.
  //////////////////////////////////////////////////////////////////////

  if (max_nm1_refit_dist3 > 0 || max_nm1_refit_distz > 0) {
    std::vector<int> refit_count(vertices->size(), 0);

    int iv = 0;
    for (v[0] = vertices->begin(); v[0] != vertices->end(); ++v[0], ++iv) {
      if (max_nm1_refit_count > 0 && refit_count[iv] >= max_nm1_refit_count)
        continue;

      const track_vec tks = vertex_track_vec(*v[0]);
      const size_t ntks = tks.size();
      if (ntks < 3)
        continue;

      if (verbose) {
        printf("doing n-%i refit on vertex at %7.4f %7.4f %7.4f with %lu tracks\n", refit_count[iv]+1, v[0]->x(), v[0]->y(), v[0]->z(), ntks);
        for (size_t i = 0; i < ntks; ++i)
          printf("  refit %lu will drop tk pt %7.4f +- %7.4f eta %7.4f +- %7.4f phi %7.4f +- %7.4f dxy %7.4f +- %7.4f dz %7.4f +- %7.4f\n", i, tks[i]->pt(), tks[i]->ptError(), tks[i]->eta(), tks[i]->etaError(), tks[i]->phi(), tks[i]->phiError(), tks[i]->dxy(), tks[i]->dxyError(), tks[i]->dz(), tks[i]->dzError());
      }

      std::vector<reco::TransientTrack> ttks(ntks-1);
      for (size_t i = 0; i < ntks; ++i) {
        for (size_t j = 0; j < ntks; ++j)
          if (j != i)
            ttks[j-(j>=i)] = tt_builder->build(tks[j]);

        reco::Vertex vnm1(TransientVertex(kv_reco->vertex(ttks)));
        const double dist3_2 = mag2(vnm1.x() - v[0]->x(), vnm1.y() - v[0]->y(), vnm1.z() - v[0]->z());
        const double distz = mag(vnm1.z() - v[0]->z());
        if (verbose) printf("  refit %lu chi2 %7.4f vtx %7.4f %7.4f %7.4f dist3 %7.4f distz %7.4f\n", i, vnm1.chi2(), vnm1.x(), vnm1.y(), vnm1.z(), sqrt(dist3_2), distz);

        if (vnm1.chi2() < 0 ||
            (max_nm1_refit_dist3 > 0 && mag2(vnm1.x() - v[0]->x(), vnm1.y() - v[0]->y(), vnm1.z() - v[0]->z()) > pow(max_nm1_refit_dist3, 2)) ||
            (max_nm1_refit_distz > 0 && distz > max_nm1_refit_distz)) {
          if (verbose) {
            printf("    replacing");
            if (refit_count[iv] < max_nm1_refit_count - 1)
              printf(" and reconsidering");
            printf("\n");
          }

          *v[0] = vnm1;
          ++refit_count[iv];
          --v[0], --iv;
          break;
        }
      }
    }

  }
  
  // Peace's scratch 
  /*
  edm::Handle<pat::JetCollection> jets;
  typedef std::vector<reco::TrackRef> jet_which_track_ref;
  std::vector<int> jet_which_jet_idx;
  for (int jjet = 0, njjet = int(jets->size()); jjet < njjet; ++jjet) {
	  const pat::Jet& jet = jets->at(jjet);

	  const size_t ijet = int(jets->size()) - 1;
	  assert(ijet <= 255);  
	  jet_which_jet_idx.push_back(jjet);
	  std::vector<int> jet_jjet_which_track_ref;
	
	  for (size_t idau = 0, ndau = jet.numberOfDaughters(); idau < ndau; ++idau) {
		  // handle both regular aod and miniaod: in the latter
		  // getPFConstituents() doesn't work because the daughters are
		  // pat::PackedCandidates. Since we don't care about track
		  // identities e.g. to compare with vertices we don't use
		  // TrackRefGetter here, but could
		  const reco::Candidate* dau = jet.daughter(idau);
		  if (dau->charge() == 0)
			  continue;

		  const reco::Track * tk = 0;
		  const reco::PFCandidate * pf = dynamic_cast<const reco::PFCandidate*>(dau);
		  if (pf) {
			  const reco::TrackRef& r = pf->trackRef();
			  if (r.isNonnull())
				  tk = &*r;
		  }
		  else {
			  const pat::PackedCandidate* pk = dynamic_cast<const pat::PackedCandidate*>(dau);
			  if (pk && pk->charge() && pk->hasTrackDetails())
				  tk = &pk->pseudoTrack();
		  }
		  if (tk) {
			  assert(abs(tk->charge()) == 1);
			  jet_jjet_which_track_ref.push_back(*tk);	  // not sure how to store address of tk
			  }
	  }
	  jet_which_track_ref.push_back(jet_jjet_which_track_ref);
	  
  }

  for (size_t ivt = 0, nvt = vertices.size(); ivt < nvt; ++ivt) {
	  const reco::Vertex& v = vertices->at(ivt);
	  tracks[ivt] = vertex_track_set(*v[ivt]);
	  for (auto tk : tracks[ivt])
		  if (tk.key() ==

	  

  }
  */


  //////////////////////////////////////////////////////////////////////
  // Put the output.
  //////////////////////////////////////////////////////////////////////

  finish(event, seed_tracks, std::move(vertices), std::move(vpeffs), vpeffs_tracks);
}

DEFINE_FWK_MODULE(MFVVertexer);
