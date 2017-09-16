// g++ -I $CMSSW_BASE/src -g -Wall `root-config --cflags --libs --glibs` ../../src/MiniNtuple.cc 2v_from_jets.cc -o 2v_from_jets.exe && ./2v_from_jets.exe

#include <cstdlib>
#include <math.h>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVector2.h"
#include "JMTucker/MFVNeutralino/interface/MiniNtuple.h"

const char* file_path = "/uscms_data/d2/tucker/crab_dirs/MiniTreeV15";

int dvv_nbins = 40;
double dvv_bin_width = 0.01;

float ht(int njets, float* jet_pt) {
  double sum = 0;
  for (int i = 0; i < njets; ++i) {
    sum += jet_pt[i];
  }
  return sum;
}

struct ConstructDvvcParameters {
  bool is_mc_;
  bool inject_signal_;
  std::string year_;
  int ntracks_;
  bool correct_bquarks_;
  int bquarks_;
  bool vary_dphi_;
  bool clearing_from_eff_;
  bool vary_eff_;
  int min_npu_;
  int max_npu_;

  ConstructDvvcParameters()
    : is_mc_(true),
      inject_signal_(false),
      year_("2015p6"),
      ntracks_(5),
      correct_bquarks_(true),
      bquarks_(-1),
      vary_dphi_(false),
      clearing_from_eff_(true),
      vary_eff_(false),
      min_npu_(0),
      max_npu_(255)
  {
  }

  bool is_mc() const { return is_mc_; }
  bool inject_signal() const { return inject_signal_; }
  std::string year() const { return year_; }
  int ntracks() const { return ntracks_; }
  bool correct_bquarks() const { return correct_bquarks_; }
  int bquarks() const { return bquarks_; }
  bool vary_dphi() const { return vary_dphi_; }
  bool clearing_from_eff() const { return clearing_from_eff_; }
  bool vary_eff() const { return vary_eff_; }
  int min_npu() const { return min_npu_; }
  int max_npu() const { return max_npu_; }

  ConstructDvvcParameters is_mc(bool x)             { ConstructDvvcParameters y(*this); y.is_mc_             = x; return y; }
  ConstructDvvcParameters inject_signal(bool x)     { ConstructDvvcParameters y(*this); y.inject_signal_     = x; return y; }
  ConstructDvvcParameters year(std::string x)       { ConstructDvvcParameters y(*this); y.year_              = x; return y; }
  ConstructDvvcParameters ntracks(int x)            { ConstructDvvcParameters y(*this); y.ntracks_           = x; return y; }
  ConstructDvvcParameters correct_bquarks(bool x)   { ConstructDvvcParameters y(*this); y.correct_bquarks_   = x; return y; }
  ConstructDvvcParameters bquarks(int x)            { ConstructDvvcParameters y(*this); y.bquarks_           = x; return y; }
  ConstructDvvcParameters vary_dphi(bool x)         { ConstructDvvcParameters y(*this); y.vary_dphi_         = x; return y; }
  ConstructDvvcParameters clearing_from_eff(bool x) { ConstructDvvcParameters y(*this); y.clearing_from_eff_ = x; return y; }
  ConstructDvvcParameters vary_eff(bool x)          { ConstructDvvcParameters y(*this); y.vary_eff_          = x; return y; }
  ConstructDvvcParameters min_npu(int x)            { ConstructDvvcParameters y(*this); y.min_npu_           = x; return y; }
  ConstructDvvcParameters max_npu(int x)            { ConstructDvvcParameters y(*this); y.max_npu_           = x; return y; }

  void print() const {
    printf("is_mc = %d, inject_signal = %d, year = %s, ntracks = %d, correct_bquarks = %d, bquarks = %d, vary_dphi = %d, clearing_from_eff = %d, vary_eff = %d", is_mc(), inject_signal(), year_.c_str(), ntracks(), correct_bquarks(), bquarks(), vary_dphi(), clearing_from_eff(), vary_eff());
  }
};

void construct_dvvc(ConstructDvvcParameters p, const char* out_fn) {
  p.print(); printf(", out_fn = %s\n", out_fn);

  const int nbkg = 20;
  const char* samples[nbkg];
  float       weights[nbkg];
  samples[0]  = "mfv_neu_tau01000um_M0800_2015"; weights[0]  = 0.0002683;
  samples[1]  = "qcdht1000sum_2015";             weights[1]  = 0.21105;
  samples[2]  = "qcdht1500sum_2015";             weights[2]  = 0.02736;
  samples[3]  = "qcdht2000sum_2015";             weights[3]  = 0.01132;
  samples[4]  = "ttbar_2015";                    weights[4]  = 0.05799;
  samples[5]  = "qcdht1000sum";                  weights[5]  = 2.84372;
  samples[6]  = "qcdht1500sum";                  weights[6]  = 0.36354;
  samples[7]  = "qcdht2000sum";                  weights[7]  = 0.15026;
  samples[8]  = "ttbar";                         weights[8]  = 0.68346;
  samples[9]  = "mfv_neu_tau01000um_M0800";      weights[9]  = 0.00035985;
  samples[10] = "JetHT2015C";                    weights[10] = 1;
  samples[11] = "JetHT2015D";                    weights[11] = 1;
  samples[12] = "JetHT2016B3";                   weights[12] = 1;
  samples[13] = "JetHT2016C";                    weights[13] = 1;
  samples[14] = "JetHT2016D";                    weights[14] = 1;
  samples[15] = "JetHT2016E";                    weights[15] = 1;
  samples[16] = "JetHT2016F";                    weights[16] = 1;
  samples[17] = "JetHT2016G";                    weights[17] = 1;
  samples[18] = "JetHT2016H2";                   weights[18] = 1;
  samples[19] = "JetHT2016H3";                   weights[19] = 1;


  int ibkg_begin;
  int ibkg_end;
  double dphi_pdf_c;
  double dphi_pdf_e = 2; // JMTBAD could be moved into params
  double dphi_pdf_a;
  const char* eff_file;
  const char* eff_hist = "maxtk3";

  if (p.is_mc()) {
    if (p.year() == "2015") {
      ibkg_begin = 1;
      if (p.inject_signal()) ibkg_begin = 0;
      ibkg_end = 4;
      dphi_pdf_c = 1.31;
      dphi_pdf_a = 4.04;
      eff_file = "vpeffs_2015_v15.root";
    } else if (p.year() == "2016") {
      ibkg_begin = 5;
      ibkg_end = 8;
      if (p.inject_signal()) ibkg_end = 9;
      dphi_pdf_c = 1.34;
      dphi_pdf_a = 3.86;
      eff_file = "vpeffs_2016_v15.root";
    } else if (p.year() == "2015p6") {
      ibkg_begin = 1;
      if (p.inject_signal()) ibkg_begin = 0;
      ibkg_end = 8;
      if (p.inject_signal()) ibkg_end = 9;
      dphi_pdf_c = 1.34;
      dphi_pdf_a = 3.87;
      eff_file = "vpeffs_2015p6_v15.root";
    } else {
      fprintf(stderr, "bad year"); exit(1);
    }
  } else {
    if (p.year() == "2015") {
      ibkg_begin = 10;
      ibkg_end = 11;
      dphi_pdf_c = 1.33;
      dphi_pdf_a = 4.35;
      eff_file = "vpeffs_data_2015_v15.root";
    } else if (p.year() == "2016") {
      ibkg_begin = 12;
      ibkg_end = 19;
      dphi_pdf_c = 1.32;
      dphi_pdf_a = 4.65;
      eff_file = "vpeffs_data_2016_v15.root";
    } else if (p.year() == "2015p6") {
      ibkg_begin = 10;
      ibkg_end = 19;
      dphi_pdf_c = 1.32;
      dphi_pdf_a = 4.64;
      eff_file = "vpeffs_data_2015p6_v15.root";
    } else if (p.year() == "2016BCD") {
      ibkg_begin = 12;
      ibkg_end = 14;
      dphi_pdf_c = 1.36;
      dphi_pdf_a = 4.59;
      eff_file = "vpeffs_data_2016BCD_v15.root";
    } else if (p.year() == "2016EF") {
      ibkg_begin = 15;
      ibkg_end = 16;
      dphi_pdf_c = 1.21;
      dphi_pdf_a = 4.83;
      eff_file = "vpeffs_data_2016EF_v15.root";
    } else if (p.year() == "2016G") {
      ibkg_begin = 17;
      ibkg_end = 17;
      dphi_pdf_c = 1.38;
      dphi_pdf_a = 4.38;
      eff_file = "vpeffs_data_2016G_v15.root";
    } else if (p.year() == "2016H") {
      ibkg_begin = 18;
      ibkg_end = 19;
      dphi_pdf_c = 1.24;
      dphi_pdf_a = 4.99;
      eff_file = "vpeffs_data_2016H_v15.root";
    } else {
      fprintf(stderr, "bad year"); exit(1);
    }
  }

  const char* tree_path;
  double bquark_correction[3] = {1,1,1};
  int min_ntracks0 = 0;
  int max_ntracks0 = 1000000;
  int min_ntracks1 = 0;
  int max_ntracks1 = 1000000;

  if (p.ntracks() == 3) {
    tree_path = "mfvMiniTreeNtk3/t";
    if (p.year() == "2015")        { bquark_correction[0] = 0.93; bquark_correction[1] = 1.06; bquark_correction[2] = 1.09; }
    else if (p.year() == "2015p6") { bquark_correction[0] = 0.93; bquark_correction[1] = 1.07; bquark_correction[2] = 1.10; }
    else                           { bquark_correction[0] = 0.94; bquark_correction[1] = 1.07; bquark_correction[2] = 1.11; }
  } else if (p.ntracks() == 4) {
    tree_path = "mfvMiniTreeNtk4/t";
    if (p.year() == "2015")        { bquark_correction[0] = 0.93; bquark_correction[1] = 1.10; bquark_correction[2] = 1.12; }
    else if (p.year() == "2015p6") { bquark_correction[0] = 0.93; bquark_correction[1] = 1.11; bquark_correction[2] = 1.20; }
    else                           { bquark_correction[0] = 0.93; bquark_correction[1] = 1.11; bquark_correction[2] = 1.19; }
  } else if (p.ntracks() == 5) {
    tree_path = "mfvMiniTree/t";
    if (p.year() == "2015")        { bquark_correction[0] = 0.90; bquark_correction[1] = 1.31; bquark_correction[2] = 1.26; }
    else if (p.year() == "2015p6") { bquark_correction[0] = 0.92; bquark_correction[1] = 1.25; bquark_correction[2] = 1.57; }
    else                           { bquark_correction[0] = 0.92; bquark_correction[1] = 1.25; bquark_correction[2] = 1.54; }
  } else if (p.ntracks() == 7) {
    tree_path = "mfvMiniTreeNtk3or4/t";
    if (p.year() == "2015")        { bquark_correction[0] = 0.93; bquark_correction[1] = 1.07; bquark_correction[2] = 1.11; }
    else if (p.year() == "2015p6") { bquark_correction[0] = 0.93; bquark_correction[1] = 1.09; bquark_correction[2] = 1.14; }
    else                           { bquark_correction[0] = 0.93; bquark_correction[1] = 1.09; bquark_correction[2] = 1.14; }
    min_ntracks0 = 4;
    max_ntracks0 = 4;
    min_ntracks1 = 3;
    max_ntracks1 = 3;
  } else {
    fprintf(stderr, "bad ntracks"); exit(1);
  }

  if (!p.correct_bquarks()) {
    bquark_correction[0] = 1; bquark_correction[1] = 1; bquark_correction[2] = 1;
  }

  if (p.vary_eff()) {
    if (p.is_mc()) {
      if (p.year() == "2015") {
        eff_file = "vpeffs_2015_v15_ntkseeds.root";
      } else if (p.year() == "2016") {
        eff_file = "vpeffs_2016_v15_ntkseeds.root";
      } else if (p.year() == "2015p6") {
        eff_file = "vpeffs_2015p6_v15_ntkseeds.root";
      }
    } else {
      if (p.year() == "2015") {
        eff_file = "vpeffs_data_2015_v15_ntkseeds.root";
      } else if (p.year() == "2016") {
        eff_file = "vpeffs_data_2016_v15_ntkseeds.root";
      } else if (p.year() == "2015p6") {
        eff_file = "vpeffs_data_2015p6_v15_ntkseeds.root";
      } else if (p.year() == "2016BCD") {
        eff_file = "vpeffs_data_2016BCD_v15_ntkseeds.root";
      } else if (p.year() == "2016EF") {
        eff_file = "vpeffs_data_2016EF_v15_ntkseeds.root";
      } else if (p.year() == "2016G") {
        eff_file = "vpeffs_data_2016G_v15_ntkseeds.root";
      } else if (p.year() == "2016H") {
        eff_file = "vpeffs_data_2016H_v15_ntkseeds.root";
      }
    }

    if (p.ntracks() == 3) {
      eff_hist = "maxtk3";
    } else if (p.ntracks() == 4) {
      eff_hist = "maxtk4";
    } else if (p.ntracks() == 5) {
      eff_hist = "maxtk5";
    } else if (p.ntracks() == 7) {
      eff_hist = "maxtk3";
    }
  }

  printf("\tdphi_pdf_c = %.2f, dphi_pdf_e = %.2f, dphi_pdf_a = %.2f, eff_file = %s, eff_hist = %s, bquark_correction = {%.2f, %.2f, %.2f}\n", dphi_pdf_c, dphi_pdf_e, dphi_pdf_a, eff_file, eff_hist, bquark_correction[0], bquark_correction[1], bquark_correction[2]);

  TH1::SetDefaultSumw2();
  gRandom->SetSeed(12191982);

  //fill only-one-vertex dBV distribution
  TH1D* h_1v_dbv = new TH1D("h_1v_dbv", "only-one-vertex events;d_{BV} (cm);events", 1250, 0, 2.5);
  TH1D* h_1v_dbv0 = new TH1D("h_1v_dbv0", "only-one-vertex events;d_{BV}^{0} (cm);events", 1250, 0, 2.5);
  TH1D* h_1v_dbv1 = new TH1D("h_1v_dbv1", "only-one-vertex events;d_{BV}^{1} (cm);events", 1250, 0, 2.5);
  TH1F* h_1v_phiv = new TH1F("h_1v_phiv", "only-one-vertex events;vertex #phi;events", 50, -3.15, 3.15);
  TH1D* h_1v_npu = new TH1D("h_1v_npu", "only-one-vertex events;# PU interactions;events", 100, 0, 100);
  TH1F* h_1v_njets = new TH1F("h_1v_njets", "only-one-vertex events;number of jets;events", 20, 0, 20);
  TH1F* h_1v_ht = new TH1F("h_1v_ht", "only-one-vertex events;#Sigma H_{T} of jets (GeV);events", 200, 0, 5000);
  TH1F* h_1v_phij = new TH1F("h_1v_phij", "only-one-vertex events;jets #phi;jets", 50, -3.15, 3.15);
  TH1F* h_1v_dphijj = new TH1F("h_1v_dphijj", "only-one-vertex events;#Delta#phi_{JJ};jet pairs", 100, -3.1416, 3.1416);
  TH1F* h_1v_dphijv = new TH1F("h_1v_dphijv", "only-one-vertex events;#Delta#phi_{JV};jet-vertex pairs", 100, -3.1416, 3.1416);
  TH1F* h_1v_dphijvpt = new TH1F("h_1v_dphijvpt", "only-one-vertex events;p_{T}-weighted #Delta#phi_{JV};jet-vertex pairs", 100, -3.1416, 3.1416);
  TH1F* h_1v_dphijvmin = new TH1F("h_1v_dphijvmin", "only-one-vertex events;#Delta#phi_{JV}^{min};events", 50, 0, 3.1416);
  TH1F* h_2v_dbv = new TH1F("h_2v_dbv", "two-vertex events;d_{BV} (cm);vertices", 1250, 0, 2.5);
  TH2F* h_2v_dbv1_dbv0 = new TH2F("h_2v_dbv1_dbv0", "two-vertex events;d_{BV}^{0} (cm);d_{BV}^{1} (cm)", 1250, 0, 2.5, 1250, 0, 2.5);
  TH1F* h_2v_dvv = new TH1F("h_2v_dvv", "two-vertex events;d_{VV} (cm);events", dvv_nbins, 0, dvv_nbins * dvv_bin_width);
  TH1F* h_2v_dphivv = new TH1F("h_2v_dphivv", "two-vertex events;#Delta#phi_{VV};events", 10, -3.15, 3.15);
  TH1F* h_2v_absdphivv = new TH1F("h_2v_absdphivv", "two-vertex events;|#Delta#phi_{VV}|;events", 5, 0, 3.15);
  TH1D* h_2v_npu = new TH1D("h_2v_npu", "two-vertex events;# PU interactions;events", 100, 0, 100);

  for (int i = ibkg_begin; i <= ibkg_end; ++i) {
    mfv::MiniNtuple nt;
    TFile* f = TFile::Open(TString::Format("%s/%s.root", file_path, samples[i]));
    if (!f || !f->IsOpen()) { fprintf(stderr, "bad file"); exit(1); }

    TTree* t = (TTree*)f->Get(tree_path);
    if (!t) { fprintf(stderr, "bad tree"); exit(1); }

    mfv::read_from_tree(t, nt);
    for (int j = 0, je = t->GetEntries(); j < je; ++j) {
      if (t->LoadTree(j) < 0) break;
      if (t->GetEntry(j) <= 0) continue;

      if ((p.bquarks() == 0 && nt.gen_flavor_code == 2) || (p.bquarks() == 1 && nt.gen_flavor_code != 2)) continue;
      if (nt.npu < p.min_npu() || nt.npu > p.max_npu()) continue;

      const float w = weights[i] * nt.weight;
      if (nt.nvtx == 1) {
        h_1v_dbv->Fill(sqrt(nt.x0*nt.x0 + nt.y0*nt.y0), w);
        if (nt.ntk0 >= min_ntracks0 && nt.ntk0 <= max_ntracks0) h_1v_dbv0->Fill(sqrt(nt.x0*nt.x0 + nt.y0*nt.y0), w);
        if (nt.ntk0 >= min_ntracks1 && nt.ntk0 <= max_ntracks1) h_1v_dbv1->Fill(sqrt(nt.x0*nt.x0 + nt.y0*nt.y0), w);
        h_1v_phiv->Fill(atan2(nt.y0,nt.x0), w);
        h_1v_npu->Fill(nt.npu, w);
        h_1v_njets->Fill(nt.njets, w);
        h_1v_ht->Fill(ht(nt.njets, nt.jet_pt), w);
        double dphijvmin = M_PI;
        for (int k = 0; k < nt.njets; ++k) {
          h_1v_phij->Fill(nt.jet_phi[k], w);
          h_1v_dphijv->Fill(TVector2::Phi_mpi_pi(atan2(nt.y0,nt.x0) - nt.jet_phi[k]), w);
          h_1v_dphijvpt->Fill(TVector2::Phi_mpi_pi(atan2(nt.y0,nt.x0) - nt.jet_phi[k]), w * (nt.jet_pt[k]/ht(nt.njets, nt.jet_pt)));
          if (fabs(TVector2::Phi_mpi_pi(atan2(nt.y0,nt.x0) - nt.jet_phi[k])) < dphijvmin) dphijvmin = fabs(TVector2::Phi_mpi_pi(atan2(nt.y0,nt.x0) - nt.jet_phi[k]));
          for (int l = k+1; l < nt.njets; ++l) {
            h_1v_dphijj->Fill(TVector2::Phi_mpi_pi(nt.jet_phi[k] - nt.jet_phi[l]), w);
          }
        }
        h_1v_dphijvmin->Fill(dphijvmin, w);
      }

      if (nt.nvtx >= 2 && nt.ntk0 >= min_ntracks0 && nt.ntk0 <= max_ntracks0 && nt.ntk1 >= min_ntracks1 && nt.ntk1 <= max_ntracks1) {
        double dbv0 = sqrt(nt.x0*nt.x0 + nt.y0*nt.y0);
        double dbv1 = sqrt(nt.x1*nt.x1 + nt.y1*nt.y1);
        h_2v_dbv->Fill(dbv0, w);
        h_2v_dbv->Fill(dbv1, w);
        h_2v_dbv1_dbv0->Fill(dbv0, dbv1, w);
        double dvv = sqrt((nt.x0-nt.x1)*(nt.x0-nt.x1) + (nt.y0-nt.y1)*(nt.y0-nt.y1));
        if (dvv > dvv_nbins * dvv_bin_width - 0.5*dvv_bin_width) dvv = dvv_nbins * dvv_bin_width - 0.5*dvv_bin_width;
        h_2v_dvv->Fill(dvv, w);
        double dphi = TVector2::Phi_mpi_pi(atan2(nt.y0,nt.x0)-atan2(nt.y1,nt.x1));
        h_2v_dphivv->Fill(dphi, w);
        h_2v_absdphivv->Fill(fabs(dphi), w);
        h_2v_npu->Fill(nt.npu, w);
      }
    }

    f->Close();
    delete f;
  }

  //construct dvvc
  TH1F* h_c1v_dbv = new TH1F("h_c1v_dbv", "constructed from only-one-vertex events;d_{BV} (cm);vertices", 1250, 0, 2.5);
  TH1F* h_c1v_dvv = new TH1F("h_c1v_dvv", "constructed from only-one-vertex events;d_{VV} (cm);events", dvv_nbins, 0, dvv_nbins * dvv_bin_width);
  TH1F* h_c1v_absdphivv = new TH1F("h_c1v_absdphivv", "constructed from only-one-vertex events;|#Delta#phi_{VV}|;events", 5, 0, 3.15);
  TH1F* h_c1v_dbv0 = new TH1F("h_c1v_dbv0", "constructed from only-one-vertex events;d_{BV}^{0} (cm);events", 1250, 0, 2.5);
  TH1F* h_c1v_dbv1 = new TH1F("h_c1v_dbv1", "constructed from only-one-vertex events;d_{BV}^{1} (cm);events", 1250, 0, 2.5);
  TH2F* h_c1v_dbv1_dbv0 = new TH2F("h_c1v_dbv1_dbv0", "constructed from only-one-vertex events;d_{BV}^{0} (cm);d_{BV}^{1} (cm)", 1250, 0, 2.5, 1250, 0, 2.5);

  TF1* f_dphi = new TF1("f_dphi", "(abs(x)-[0])**[1] + [2]", 0, M_PI);
  f_dphi->SetParameters(dphi_pdf_c, dphi_pdf_e, dphi_pdf_a);

  TF1* i_dphi = 0;
  TF1* i_dphi2 = 0;
  if (p.vary_dphi()) {
    i_dphi = new TF1("i_dphi", "((1/([1]+1))*(x-[0])**([1]+1) + [2]*x - (1/([1]+1))*(-[0])**([1]+1)) / ((1/([1]+1))*(3.14159-[0])**([1]+1) + [2]*3.14159 - (1/([1]+1))*(-[0])**([1]+1))", 0, M_PI);
    i_dphi->SetParameters(dphi_pdf_c, dphi_pdf_e, dphi_pdf_a);
    i_dphi2 = new TF1("i_dphi2", "x/3.14159", 0, M_PI);
  }

  TH1F* h_eff = 0;
  if (p.clearing_from_eff()) {
    h_eff = (TH1F*)TFile::Open(eff_file)->Get(eff_hist);
    h_eff->SetBinContent(h_eff->GetNbinsX()+1, h_eff->GetBinContent(h_eff->GetNbinsX()));
  }

  int bin1 = 0;
  int bin2 = 0;
  int bin3 = 0;
  int intobin1 = 0;
  int intobin2 = 0;
  int intobin3 = 0;
  int outofbin1 = 0;
  int outofbin2 = 0;
  int outofbin3 = 0;

  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < int(h_1v_dbv->GetEntries()); ++j) {
      double dbv0 = h_1v_dbv0->GetRandom();
      double dbv1 = h_1v_dbv1->GetRandom();
      h_c1v_dbv->Fill(dbv0);
      h_c1v_dbv->Fill(dbv1);

      double dphi = f_dphi->GetRandom();

      double dvvc = sqrt(dbv0*dbv0 + dbv1*dbv1 - 2*dbv0*dbv1*cos(dphi));

      if (p.vary_dphi()) {
        double dphi2 = i_dphi2->GetX(i_dphi->Eval(dphi), 0, M_PI);
        double dvvc2 = sqrt(dbv0*dbv0 + dbv1*dbv1 - 2*dbv0*dbv1*cos(dphi2));
        if (dvvc < 0.04) ++bin1;
        if (dvvc >= 0.04 && dvvc < 0.07) ++bin2;
        if (dvvc >= 0.07) ++bin3;
        if (!(dvvc < 0.04) && (dvvc2 < 0.04)) ++intobin1;
        if (!(dvvc >= 0.04 && dvvc < 0.07) && (dvvc2 >= 0.04 && dvvc2 < 0.07)) ++intobin2;
        if (!(dvvc >= 0.07) && (dvvc2 >= 0.07)) ++intobin3;
        if ((dvvc < 0.04) && !(dvvc2 < 0.04)) ++outofbin1;
        if ((dvvc >= 0.04 && dvvc < 0.07) && !(dvvc2 >= 0.04 && dvvc2 < 0.07)) ++outofbin2;
        if ((dvvc >= 0.07) && !(dvvc2 >= 0.07)) ++outofbin3;
        dphi = dphi2;
        dvvc = dvvc2;
      }

      double prob = 1;
      if (p.clearing_from_eff()) {
        prob = h_eff->GetBinContent(h_eff->FindBin(dvvc));
      }

      if (dvvc > dvv_nbins * dvv_bin_width - 0.5*dvv_bin_width) dvvc = dvv_nbins * dvv_bin_width - 0.5*dvv_bin_width;
      h_c1v_dvv->Fill(dvvc, prob);
      h_c1v_absdphivv->Fill(fabs(dphi), prob);
      h_c1v_dbv0->Fill(dbv0, prob);
      h_c1v_dbv1->Fill(dbv1, prob);
      h_c1v_dbv1_dbv0->Fill(dbv0, dbv1, prob);
    }
  }

  for (int i = 1; i <= h_c1v_dvv->GetNbinsX(); ++i) {
    if (h_c1v_dvv->GetBinLowEdge(i) < 0.04) {
      h_c1v_dvv->SetBinContent(i, h_c1v_dvv->GetBinContent(i) * bquark_correction[0]);
    } else if (h_c1v_dvv->GetBinLowEdge(i) < 0.07) {
      h_c1v_dvv->SetBinContent(i, h_c1v_dvv->GetBinContent(i) * bquark_correction[1]);
    } else {
      h_c1v_dvv->SetBinContent(i, h_c1v_dvv->GetBinContent(i) * bquark_correction[2]);
    }
  }

  if (p.vary_dphi()) {
    printf("bin1 = %d, bin2 = %d, bin3 = %d, intobin1 = %d, intobin2 = %d, intobin3 = %d, outofbin1 = %d, outofbin2 = %d, outofbin3 = %d\n", bin1, bin2, bin3, intobin1, intobin2, intobin3, outofbin1, outofbin2, outofbin3);
    printf("uncorrelated variation / default (bin 1): %f +/- %f\n", 1 + (intobin1 - outofbin1) / (1.*bin1), sqrt(bin1 + bin1 + intobin1 - outofbin1) / bin1);
    printf("  correlated variation / default (bin 1): %f +/- %f\n", 1 + (intobin1 - outofbin1) / (1.*bin1), sqrt(intobin1 + outofbin1) / bin1);
    printf("uncertainty correlated / uncorrelated (bin 1): %f\n", sqrt(intobin1 + outofbin1) / sqrt(bin1 + bin1 + intobin1 - outofbin1));
    printf("uncorrelated variation / default (bin 2): %f +/- %f\n", 1 + (intobin2 - outofbin2) / (1.*bin2), sqrt(bin2 + bin2 + intobin2 - outofbin2) / bin2);
    printf("  correlated variation / default (bin 2): %f +/- %f\n", 1 + (intobin2 - outofbin2) / (1.*bin2), sqrt(intobin2 + outofbin2) / bin2);
    printf("uncertainty correlated / uncorrelated (bin 2): %f\n", sqrt(intobin2 + outofbin2) / sqrt(bin2 + bin2 + intobin2 - outofbin2));
    printf("uncorrelated variation / default (bin 3): %f +/- %f\n", 1 + (intobin3 - outofbin3) / (1.*bin3), sqrt(bin3 + bin3 + intobin3 - outofbin3) / bin3);
    printf("  correlated variation / default (bin 3): %f +/- %f\n", 1 + (intobin3 - outofbin3) / (1.*bin3), sqrt(intobin3 + outofbin3) / bin3);
    printf("uncertainty correlated / uncorrelated (bin 3): %f\n", sqrt(intobin3 + outofbin3) / sqrt(bin3 + bin3 + intobin3 - outofbin3));
  }

  TFile* fh = TFile::Open(out_fn, "recreate");

  h_1v_dbv->Write();
  h_1v_dbv0->Write();
  h_1v_dbv1->Write();
  h_1v_phiv->Write();
  h_1v_npu->Write();
  h_1v_njets->Write();
  h_1v_ht->Write();
  h_1v_phij->Write();
  h_1v_dphijj->Write();
  h_1v_dphijv->Write();
  h_1v_dphijvpt->Write();
  h_1v_dphijvmin->Write();
  h_2v_dbv->Write();
  h_2v_dbv1_dbv0->Write();
  h_2v_dvv->Write();
  h_2v_dphivv->Write();
  h_2v_absdphivv->Write();
  h_2v_npu->Write();

  h_c1v_dbv->Write();
  h_c1v_dvv->Scale(1./h_c1v_dvv->Integral());
  h_c1v_dvv->Write();
  h_c1v_absdphivv->Write();
  h_c1v_dbv0->Write();
  h_c1v_dbv1->Write();
  h_c1v_dbv1_dbv0->Write();

  TCanvas* c_dvv = new TCanvas("c_dvv", "c_dvv", 700, 700);
  TLegend* l_dvv = new TLegend(0.35,0.75,0.85,0.85);
  h_2v_dvv->SetTitle(";d_{VV} (cm);events");
  h_2v_dvv->SetLineColor(kBlue);
  h_2v_dvv->SetLineWidth(3);
  h_2v_dvv->Scale(1./h_2v_dvv->Integral());
  h_2v_dvv->SetStats(0);
  h_2v_dvv->Draw();
  l_dvv->AddEntry(h_2v_dvv, "two-vertex events");
  h_c1v_dvv->SetLineColor(kRed);
  h_c1v_dvv->SetLineWidth(3);
  h_c1v_dvv->Scale(1./h_c1v_dvv->Integral());
  h_c1v_dvv->SetStats(0);
  h_c1v_dvv->Draw("sames");
  l_dvv->AddEntry(h_c1v_dvv, "constructed from only-one-vertex events");
  l_dvv->SetFillColor(0);
  l_dvv->Draw();
  c_dvv->SetTickx();
  c_dvv->SetTicky();
  c_dvv->Write();

  TCanvas* c_absdphivv = new TCanvas("c_absdphivv", "c_absdphivv", 700, 700);
  TLegend* l_absdphivv = new TLegend(0.25,0.75,0.75,0.85);
  h_2v_absdphivv->SetTitle(";|#Delta#phi_{VV}|;events");
  h_2v_absdphivv->SetLineColor(kBlue);
  h_2v_absdphivv->SetLineWidth(3);
  h_2v_absdphivv->Scale(1./h_2v_absdphivv->Integral());
  h_2v_absdphivv->SetStats(0);
  h_2v_absdphivv->Draw();
  l_absdphivv->AddEntry(h_2v_absdphivv, "two-vertex events");
  h_c1v_absdphivv->SetLineColor(kRed);
  h_c1v_absdphivv->SetLineWidth(3);
  h_c1v_absdphivv->Scale(1./h_c1v_absdphivv->Integral());
  h_c1v_absdphivv->SetStats(0);
  h_c1v_absdphivv->Draw("sames");
  l_absdphivv->AddEntry(h_c1v_absdphivv, "constructed from only-one-vertex events");
  l_absdphivv->SetFillColor(0);
  l_absdphivv->Draw();
  c_absdphivv->SetTickx();
  c_absdphivv->SetTicky();
  c_absdphivv->Write();

  f_dphi->Write();
  if (p.clearing_from_eff()) {
    h_eff->SetName("h_eff");
    h_eff->Write();
  }
  if (p.vary_dphi()) {
    i_dphi->Write();
    i_dphi2->Write();
  }

  fh->Close();

  delete h_1v_dbv;
  delete h_1v_dbv0;
  delete h_1v_dbv1;
  delete h_1v_phiv;
  delete h_1v_npu;
  delete h_1v_njets;
  delete h_1v_ht;
  delete h_1v_phij;
  delete h_1v_dphijj;
  delete h_1v_dphijv;
  delete h_1v_dphijvpt;
  delete h_1v_dphijvmin;
  delete h_2v_dbv;
  delete h_2v_dbv1_dbv0;
  delete h_2v_dvv;
  delete h_2v_dphivv;
  delete h_2v_absdphivv;
  delete h_2v_npu;
  delete c_dvv;
  delete c_absdphivv;
  delete h_c1v_dbv;
  delete h_c1v_dvv;
  delete h_c1v_absdphivv;
  delete h_c1v_dbv0;
  delete h_c1v_dbv1;
  delete h_c1v_dbv1_dbv0;
}

int main(int argc, const char* argv[]) {
  const bool only_default = strcmp(getenv("USER"), "tucker") == 0 || (argc >= 2 && strcmp(argv[1], "only_default"));
  ConstructDvvcParameters pars;
  pars.print(); printf("\n");
  ConstructDvvcParameters pars2 = pars.year("ekjafdskjds").ntracks(1359);
  pars2.print(); printf("\n");
  return 0;

  if (only_default) {
    construct_dvvc(pars, "2v_from_jets_2015p6_5track_default_v15.root");
    return 0;
  }

  for (const char* year : {"2015", "2016", "2015p6"}) {
    for (int ntracks : {3, 4, 5, 7}) {
      ConstructDvvcParameters pars2 = pars.year(year).ntracks(ntracks);
      construct_dvvc(pars2.correct_bquarks(false),              TString::Format("2v_from_jets_%s_%dtrack_bquark_uncorrected_v15.root", year, ntracks));
      construct_dvvc(pars2.correct_bquarks(false).bquarks(1),   TString::Format("2v_from_jets_%s_%dtrack_bquarks_v15.root", year, ntracks));
      construct_dvvc(pars2.correct_bquarks(false).bquarks(0),   TString::Format("2v_from_jets_%s_%dtrack_nobquarks_v15.root", year, ntracks));
      construct_dvvc(pars2,                                     TString::Format("2v_from_jets_%s_%dtrack_default_v15.root", year, ntracks));
      construct_dvvc(pars2.vary_dphi(true),                     TString::Format("2v_from_jets_%s_%dtrack_vary_dphi_v15.root", year, ntracks));
      construct_dvvc(pars2.clearing_from_eff(false),            TString::Format("2v_from_jets_%s_%dtrack_noclearing_v15.root", year, ntracks));
      construct_dvvc(pars2.vary_eff(true),                      TString::Format("2v_from_jets_%s_%dtrack_vary_eff_v15.root", year, ntracks));
      construct_dvvc(pars2.inject_signal(true),                 TString::Format("2v_from_jets_%s_%dtrack_inject_signal_v15.root", year, ntracks));
    }
  }
  for (const char* year : {"2015", "2016", "2015p6", "2016BCD", "2016EF", "2016G", "2016H"}) {
    for (int ntracks : {3, 4, 5, 7}) {
      ConstructDvvcParameters pars2 = pars.year(year).ntracks(ntracks).is_mc(false);
      construct_dvvc(pars2,                 TString::Format("2v_from_jets_data_%s_%dtrack_default_v15.root", year, ntracks));
      //construct_dvvc(pars2.vary_dphi(true), TString::Format("2v_from_jets_data_%s_%dtrack_vary_dphi_v15.root", year, ntracks));
      //construct_dvvc(pars2.vary_eff (true), TString::Format("2v_from_jets_data_%s_%dtrack_vary_eff_v15.root", year, ntracks));
    }
  }
}
