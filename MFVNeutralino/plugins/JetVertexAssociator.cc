#include "TH2F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class MFVJetVertexAssociator : public edm::EDProducer {
public:
  MFVJetVertexAssociator(const edm::ParameterSet&);

  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  typedef edm::AssociationMap<edm::OneToMany<reco::VertexCollection, pat::JetCollection> > JetVertexAssociation;

  const edm::InputTag jet_src;
  const edm::InputTag vertex_src;
  const std::string tag_info_name;
  const double min_vertex_track_weight;
  const int min_tracks_shared;
  const double min_track_pt;
  const int min_hits_shared;
  const double max_cos_angle_diff;
  const double max_miss_dist;
  const double max_miss_sig;
  const bool histos;
  const bool verbose;

  TH2F* h_n_jets_v_vertices;
  TH1F* h_n_vertex_tracks;
  TH1F* h_n_jet_tracks;

  TH1F* h_ntracks;
  TH2F* h_best_ntracks;
  TH1F* h_ntracks_ptmin;
  TH2F* h_best_ntracks_ptmin;
  TH1F* h_sum_nhits;
  TH2F* h_best_sum_nhits;
  TH1F* h_cos_angle;
  TH2F* h_best_cos_angle;
  TH1F* h_miss_dist;
  TH2F* h_best_miss_dist;
  TH1F* h_miss_dist_err;
  TH2F* h_best_miss_dist_err;
  TH1F* h_miss_dist_sig;
  TH2F* h_best_miss_dist_sig;
  TH2F* h_miss_dist_err_v;
  TH2F* h_best_miss_dist_err_v;

  TH2F* h_n_matchedjets_v_jets;
  TH2F* h_n_matchedjets_v_vertices;
  TH2F* h_n_matchedvertices_v_jets;
  TH2F* h_n_matchedvertices_v_vertices;
};

MFVJetVertexAssociator::MFVJetVertexAssociator(const edm::ParameterSet& cfg)
  : jet_src(cfg.getParameter<edm::InputTag>("jet_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    tag_info_name(cfg.getParameter<std::string>("tag_info_name")),
    min_vertex_track_weight(cfg.getParameter<double>("min_vertex_track_weight")),
    min_tracks_shared(cfg.getParameter<int>("min_tracks_shared")),
    min_track_pt(cfg.getParameter<double>("min_track_pt")),
    min_hits_shared(cfg.getParameter<int>("min_hits_shared")),
    max_cos_angle_diff(cfg.getParameter<double>("max_cos_angle_diff")),
    max_miss_dist(cfg.getParameter<double>("max_miss_dist")),
    max_miss_sig(cfg.getParameter<double>("max_miss_sig")),
    histos(cfg.getUntrackedParameter<bool>("histos", false)),
    verbose(cfg.getUntrackedParameter<bool>("verbose", false))
{
  produces<JetVertexAssociation>();

  if (histos) {
    edm::Service<TFileService> fs;
    
    h_n_jets_v_vertices = fs->make<TH2F>("h_n_jets_v_vertices", ";# of vertices;# of jets", 20, 0, 20, 20, 0, 20);
    h_n_vertex_tracks = fs->make<TH1F>("h_n_vertex_tracks", ";# tracks per vertex;arb. units", 50, 0, 50);
    h_n_jet_tracks = fs->make<TH1F>("h_n_jet_tracks", ";# tracks per jet;arb. units", 50, 0, 50);

    h_ntracks = fs->make<TH1F>("h_ntracks", ";# tracks shared w. a vertex;arb. units", 20, 0, 20);
    h_best_ntracks = fs->make<TH2F>("h_best_ntracks", ";# tracks shared w. 2nd-best vertex;# tracks shared w. best vertex", 20, 0, 20, 20, 0, 20);
    h_ntracks_ptmin = fs->make<TH1F>("h_ntracks_ptmin", ";# tracks (p_{T} cut) shared w. a vertex;arb. units", 20, 0, 20);
    h_best_ntracks_ptmin = fs->make<TH2F>("h_best_ntracks_ptmin", ";# tracks (p_{T} cut) shared w. 2nd-best vertex;# tracks (p_{T} cut) shared w. best vertex", 20, 0, 20, 20, 0, 20);
    h_sum_nhits = fs->make<TH1F>("h_sum_nhits", ";# tracks' hits shared w. a vertex;arb. units", 100, 0, 100);
    h_best_sum_nhits = fs->make<TH2F>("h_best_sum_nhits", ";# tracks' hits shared w. 2nd-best vertex;# tracks' hits shared w. best vertex", 100, 0, 100, 100, 0, 100);
    h_cos_angle = fs->make<TH1F>("h_cos_angle", ";cos(angle between jet mom. and TV-SV);arb. units", 100, -1, 1);
    h_best_cos_angle = fs->make<TH2F>("h_best_cos_angle", ";cos(angle between jet mom. and TV-SV) for 2nd-best vertex;cos(angle between jet mom. and TV-SV) for best vertex", 100, 0, 100, 100, 0, 100);
    h_miss_dist = fs->make<TH1F>("h_miss_dist", ";jet miss distance (cm);arb. units", 100, 0, 2);
    h_best_miss_dist = fs->make<TH2F>("h_best_miss_dist", ";jet miss distance to 2nd-best vertex (cm);jet miss distance to best vertex", 100, 0, 2, 100, 0, 2);
    h_miss_dist_err = fs->make<TH1F>("h_miss_dist_err", ";#sigma(jet miss distance) (cm);arb. units", 100, 0, 2);
    h_best_miss_dist_err = fs->make<TH2F>("h_best_miss_dist_err", ";#sigma(jet miss distance to 2nd-best vertex) (cm);#sigma(jet miss distance to best vertex) (cm)", 100, 0, 2, 100, 0, 2);
    h_miss_dist_sig = fs->make<TH1F>("h_miss_dist_sig", ";N#sigma(jet miss distance);arb. units", 100, 0, 20);
    h_best_miss_dist_sig = fs->make<TH2F>("h_best_miss_dist_sig", ";N#sigma(jet miss distance to 2nd-best vertex);N#sigma(jet miss distance to best vertex)", 100, 0, 20, 100, 0, 20);
    h_miss_dist_err_v = fs->make<TH2F>("h_miss_dist_err_v", ";jet miss distance to a vertex (cm);#sigma(jet miss distance to a vertex) (cm)", 100, 0, 2, 100, 0, 2);
    h_best_miss_dist_err_v = fs->make<TH2F>("h_best_miss_dist_err_v", ";jet miss distance to best vertex (cm);#sigma(jet miss distance to best vertex) (cm)", 100, 0, 2, 100, 0, 2);

    h_n_matchedjets_v_jets = fs->make<TH2F>("h_n_matchedjets_v_jets", ";# of jets;# of matched jets", 20, 0, 20, 20, 0, 20);
    h_n_matchedjets_v_vertices = fs->make<TH2F>("h_n_matchedjets_v_vertices", ";# of vertices;# of matched jets", 20, 0, 20, 20, 0, 20);
    h_n_matchedvertices_v_jets = fs->make<TH2F>("h_n_matchedvertices_v_jets", ";# of jets;# of matched vertices", 20, 0, 20, 20, 0, 20);
    h_n_matchedvertices_v_vertices = fs->make<TH2F>("h_n_matchedvertices_v_vertices", ";# of vertices;# of matched vertices", 20, 0, 20, 20, 0, 20);
  }
}

void MFVJetVertexAssociator::produce(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jet_src, jets);

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel(vertex_src, vertices);

  const size_t n_jets = jets->size();
  const size_t n_vertices = vertices->size();

  if (histos) {
    h_n_jets_v_vertices->Fill(n_vertices, n_jets);
    for (const reco::Vertex& vtx : *vertices)
      h_n_vertex_tracks->Fill(vtx.nTracks(min_vertex_track_weight));
  }

  if (verbose) {
    for (size_t ivtx = 0; ivtx < n_vertices; ++ivtx) {
      const reco::Vertex& vtx = vertices->at(ivtx);
      printf("ivtx %lu ntracks %i mass %f\n", ivtx, vtx.nTracks(min_vertex_track_weight), vtx.p4().mass());
    }
  }

  // Associate jets to vertices. For each jet, find the vertex that
  // shares the most tracks, or the most tracks with given ptmin, or
  // the most tracks' hits.
  //
  // When considering #tracks, do not use fractions of tracks so that
  // low-ntrack vertices are not preferentially picked.
  //
  // Also, especially for the case when a jet does not share tracks
  // with any vertex, try to find a vertex (SV) that it points back to
  // in a straight line. The two points that form the line are the SV
  // position and the jet's vertex. The jet's vertex is taken from the
  // b-tag software "SV" (= "TV" here) tag for now. (JMTBAD try
  // fitting jet's tracks for vertices using our fitter, and taking
  // the one with e.g. the largest number of tracks, or pT.)
  //
  // Try taking closest in cos(angle between jet momentum and
  // (TV-SV)), as well as the distance of closest approach ("miss
  // distance").

  std::vector<int> index_by_ntracks(n_jets, -1);
  std::vector<int> best_ntracks(n_jets, 0);
  std::vector<int> second_best_ntracks(n_jets, 0);

  std::vector<int> index_by_ntracks_ptmin(n_jets, -1);
  std::vector<int> best_ntracks_ptmin(n_jets, 0);
  std::vector<int> second_best_ntracks_ptmin(n_jets, 0);

  std::vector<int> index_by_sum_nhits(n_jets, -1);
  std::vector<int> best_sum_nhits(n_jets, 0);
  std::vector<int> second_best_sum_nhits(n_jets, 0);

  std::vector<int> index_by_cos_angle(n_jets, -1);
  std::vector<double> best_cos_angle(n_jets, 1e9);
  std::vector<double> second_best_cos_angle(n_jets, 1e9);

  std::vector<int> index_by_miss_dist(n_jets);
  std::vector<Measurement1D> best_miss_dist(n_jets, Measurement1D(1e9, 1e9));
  std::vector<Measurement1D> second_best_miss_dist(n_jets, Measurement1D(1e9, 1e9));

  for (size_t ijet = 0; ijet < n_jets; ++ijet) {
    const pat::Jet& jet = jets->at(ijet);
    std::set<reco::TrackRef> jet_tracks;
    for (const reco::PFCandidatePtr& pfcand : jet.getPFConstituents()) {
      const reco::TrackRef& tk = pfcand->trackRef();
      if (tk.isNonnull())
        jet_tracks.insert(tk);
    }

    const size_t n_jet_tracks = jet_tracks.size();
    if (histos)
      h_n_jet_tracks->Fill(n_jet_tracks);

    if (n_jet_tracks == 0)
      continue;

    for (size_t ivtx = 0; ivtx < n_vertices; ++ivtx) {
      const reco::Vertex& vtx = vertices->at(ivtx);
      int ntracks = 0;
      int ntracks_ptmin = 0;
      int sum_nhits = 0;
      double cos_angle = 1e9;
      Measurement1D miss_dist(1e9, 1e9);

      for (auto itk = vtx.tracks_begin(), itke = vtx.tracks_end(); itk != itke; ++itk) {
        if (vtx.trackWeight(*itk) >= min_vertex_track_weight) {
          reco::TrackRef tk = itk->castTo<reco::TrackRef>();
          if (jet_tracks.count(tk) > 0) {
            ++ntracks;
            if (tk->pt() > min_track_pt)
              ++ntracks_ptmin;
            sum_nhits += tk->hitPattern().numberOfValidHits();
          }
        }
      }

      if (ntracks > best_ntracks[ijet]) {
        second_best_ntracks[ijet] = best_ntracks[ijet];
        best_ntracks[ijet] = ntracks;
        if (ntracks >= min_tracks_shared)
          index_by_ntracks[ijet] = ivtx;
      }

      if (ntracks_ptmin > best_ntracks_ptmin[ijet]) {
        second_best_ntracks_ptmin[ijet] = best_ntracks_ptmin[ijet];
        best_ntracks_ptmin[ijet] = ntracks_ptmin;
        if (ntracks_ptmin >= min_tracks_shared)
          index_by_ntracks_ptmin[ijet] = ivtx;
      }

      if (sum_nhits > best_sum_nhits[ijet]) {
        second_best_sum_nhits[ijet] = best_sum_nhits[ijet];
        best_sum_nhits[ijet] = sum_nhits;
        if (sum_nhits >= min_hits_shared)
          index_by_sum_nhits[ijet] = ivtx;
      }
        
      if (fabs(cos_angle - 1) < fabs(best_cos_angle[ijet] - 1)) {
        second_best_cos_angle[ijet] = best_cos_angle[ijet];
        best_cos_angle[ijet] = cos_angle;
        if (fabs(cos_angle - 1) <= max_cos_angle_diff)
          index_by_cos_angle[ijet] = ivtx;
      }

      if (miss_dist.value() < best_miss_dist[ijet].value()) {
        second_best_miss_dist[ijet] = best_miss_dist[ijet];
        best_miss_dist[ijet] = miss_dist;
        if (miss_dist.value() <= max_miss_dist && miss_dist.significance() <= max_miss_sig)
          index_by_miss_dist[ijet] = ivtx;
      }

      if (histos) {
        h_ntracks->Fill(ntracks);
        h_ntracks_ptmin->Fill(ntracks_ptmin);
        h_sum_nhits->Fill(sum_nhits);
        h_cos_angle->Fill(cos_angle);
        h_miss_dist->Fill(miss_dist.value());
        h_miss_dist_err->Fill(miss_dist.error());
        h_miss_dist_sig->Fill(miss_dist.significance());
        h_miss_dist_err_v->Fill(miss_dist.value(), miss_dist.error());
      }        
    }        

    if (histos) {
      h_best_ntracks->Fill(second_best_ntracks[ijet], best_ntracks[ijet]);
      h_best_ntracks_ptmin->Fill(second_best_ntracks_ptmin[ijet], best_ntracks_ptmin[ijet]);
      h_best_sum_nhits->Fill(second_best_sum_nhits[ijet], best_sum_nhits[ijet]);
      h_best_cos_angle->Fill(second_best_cos_angle[ijet], best_cos_angle[ijet]);
      h_best_miss_dist->Fill(second_best_miss_dist[ijet].value(), best_miss_dist[ijet].value());
      h_best_miss_dist_sig->Fill(second_best_miss_dist[ijet].significance(), best_miss_dist[ijet].significance());
      h_best_miss_dist_err_v->Fill(best_miss_dist[ijet].value(), best_miss_dist[ijet].error());
    }        

    //if (verbose)
    //  printf("ijet %lu pt %f eta %f phi %f  assoc ivtx %i\n", ijet, jet.pt(), jet.eta(), jet.phi(), indices[ijet]);
  }


  std::auto_ptr<JetVertexAssociation> assoc(new JetVertexAssociation);
  int n_matchedvertices = 0;
  int n_matchedjets = 0;

  for (size_t ivtx = 0; ivtx < n_vertices; ++ivtx) {
    reco::VertexRef vtxref(vertices, ivtx);
    int these_n_matchedjets = 0;

    for (size_t ijet = 0; ijet < n_jets; ++ijet) {
      pat::JetRef jetref(jets, ijet);
      if (index_by_ntracks[ijet] == int(ivtx)) {
        assoc->insert(vtxref, jetref);
        ++these_n_matchedjets;
      }
    }

    n_matchedjets += these_n_matchedjets;
    if (these_n_matchedjets > 0)
      ++n_matchedvertices;
  }

  if (histos) {
    h_n_matchedjets_v_jets->Fill(n_jets, n_matchedjets);
    h_n_matchedjets_v_vertices->Fill(n_vertices, n_matchedjets);
    h_n_matchedvertices_v_jets->Fill(n_jets, n_matchedvertices);
    h_n_matchedvertices_v_vertices->Fill(n_vertices, n_matchedvertices);
  }

  event.put(assoc);
}

DEFINE_FWK_MODULE(MFVJetVertexAssociator);
