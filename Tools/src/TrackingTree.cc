#include "JMTucker/Tools/interface/TrackingTree.h"
#include <cassert>
#include "TTree.h"

TrackingTree::TrackingTree() {
  clear();
  p_pv_x_ = 0;
  p_pv_y_ = 0;
  p_pv_z_ = 0;
  p_pv_chi2dof_ = 0;
  p_pv_ndof_ = 0;
  p_pv_score_ = 0;
  p_pv_cxx_ = 0;
  p_pv_cxy_ = 0;
  p_pv_cxz_ = 0;
  p_pv_cyy_ = 0;
  p_pv_cyz_ = 0;
  p_pv_czz_ = 0;
  p_tk_qpt_ = 0;
  p_tk_eta_ = 0;
  p_tk_phi_ = 0;
  p_tk_dxybs_ = 0;
  p_tk_dxypv_ = 0;
  p_tk_dzpv_ = 0;
  p_tk_vx_ = 0;
  p_tk_vy_ = 0;
  p_tk_vz_ = 0;
  p_tk_err_pt_ = 0;
  p_tk_err_eta_ = 0;
  p_tk_err_phi_ = 0;
  p_tk_err_dxy_ = 0;
  p_tk_err_dz_ = 0;
  p_tk_chi2dof_ = 0;
  p_tk_hp_ = 0;
  p_tk_minhit_ = 0;
  p_tk_maxhit_ = 0;
  p_tk_maxpxhit_ = 0;
}

void TrackingTree::clear() {
  run_ = 0;
  lumi_ = 0;
  event_ = 0;
  npu_ = 0;
  bs_x_ = 0;
  bs_y_ = 0;
  bs_z_ = 0;
  bs_sigmaz_ = 0;
  bs_dxdz_ = 0;
  bs_dydz_ = 0;
  bs_width_ = 0;
  bs_err_x_ = 0;
  bs_err_y_ = 0;
  bs_err_z_ = 0;
  bs_err_sigmaz_ = 0;
  bs_err_dxdz_ = 0;
  bs_err_dydz_ = 0;
  bs_err_width_ = 0;
  pv_x_.clear();
  pv_y_.clear();
  pv_z_.clear();
  pv_chi2dof_.clear();
  pv_ndof_.clear();
  pv_score_.clear();
  pv_cxx_.clear();
  pv_cxy_.clear();
  pv_cxz_.clear();
  pv_cyy_.clear();
  pv_cyz_.clear();
  pv_czz_.clear();
  tk_qpt_.clear();
  tk_eta_.clear();
  tk_phi_.clear();
  tk_dxybs_.clear();
  tk_dxypv_.clear();
  tk_dzpv_.clear();
  tk_vx_.clear();
  tk_vy_.clear();
  tk_vz_.clear();
  tk_err_pt_.clear();
  tk_err_eta_.clear();
  tk_err_phi_.clear();
  tk_err_dxy_.clear();
  tk_err_dz_.clear();
  tk_chi2dof_.clear();
  tk_hp_.clear();
  tk_minhit_.clear();
  tk_maxhit_.clear();
  tk_maxpxhit_.clear();
}

void TrackingTree::write_to_tree(TTree* t) {
  t->SetAlias("npvs", "pv_x@.size()");
  t->SetAlias("ntks", "tk_qpt@.size()");
  t->SetAlias("tk_q", "tk_qpt > 0 ? 1 : -1");
  t->SetAlias("tk_pt", "abs(tk_qpt)");
  t->SetAlias("tk_npxhits", "tk_hp & 0x7");
  t->SetAlias("tk_nsthits", "(tk_hp >> 3) & 0x1f");
  t->SetAlias("tk_npxlayers", "(tk_hp >> 8) & 0x7");
  t->SetAlias("tk_nstlayers", "(tk_hp >> 11) & 0x1f");
  t->SetAlias("tk_nhits", "tk_npxhits + tk_nsthits");
  t->SetAlias("tk_nlayers", "tk_npxlayers + tk_nstlayers");
  t->SetAlias("tk_min_r", "tk_minhit & 0xf");
  t->SetAlias("tk_min_z", "tk_minhit >> 4");
  t->SetAlias("tk_max_r", "tk_maxhit & 0xf");
  t->SetAlias("tk_max_z", "tk_maxhit >> 4");
  t->SetAlias("tk_maxpx_r", "tk_maxpxhit & 0xf");
  t->SetAlias("tk_maxpx_z", "tk_maxpxhit >> 4");

  t->Branch("run", &run_);
  t->Branch("lumi", &lumi_);
  t->Branch("event", &event_);
  t->Branch("npu", &npu_);
  t->Branch("bs_x", &bs_x_);
  t->Branch("bs_y", &bs_y_);
  t->Branch("bs_z", &bs_z_);
  t->Branch("bs_sigmaz", &bs_sigmaz_);
  t->Branch("bs_dxdz", &bs_dxdz_);
  t->Branch("bs_dydz", &bs_dydz_);
  t->Branch("bs_width", &bs_width_);
  t->Branch("bs_err_x", &bs_err_x_);
  t->Branch("bs_err_y", &bs_err_y_);
  t->Branch("bs_err_z", &bs_err_z_);
  t->Branch("bs_err_sigmaz", &bs_err_sigmaz_);
  t->Branch("bs_err_dxdz", &bs_err_dxdz_);
  t->Branch("bs_err_dydz", &bs_err_dydz_);
  t->Branch("bs_err_width", &bs_err_width_);
  t->Branch("pv_x", &pv_x_);
  t->Branch("pv_y", &pv_y_);
  t->Branch("pv_z", &pv_z_);
  t->Branch("pv_chi2dof", &pv_chi2dof_);
  t->Branch("pv_ndof", &pv_ndof_);
  t->Branch("pv_score", &pv_score_);
  t->Branch("pv_cxx", &pv_cxx_);
  t->Branch("pv_cxy", &pv_cxy_);
  t->Branch("pv_cxz", &pv_cxz_);
  t->Branch("pv_cyy", &pv_cyy_);
  t->Branch("pv_cyz", &pv_cyz_);
  t->Branch("pv_czz", &pv_czz_);
  t->Branch("tk_qpt", &tk_qpt_);
  t->Branch("tk_eta", &tk_eta_);
  t->Branch("tk_phi", &tk_phi_);
  t->Branch("tk_dxybs", &tk_dxybs_);
  t->Branch("tk_dxypv", &tk_dxypv_);
  t->Branch("tk_dzpv", &tk_dzpv_);
  t->Branch("tk_vx", &tk_vx_);
  t->Branch("tk_vy", &tk_vy_);
  t->Branch("tk_vz", &tk_vz_);
  t->Branch("tk_err_pt", &tk_err_pt_);
  t->Branch("tk_err_eta", &tk_err_eta_);
  t->Branch("tk_err_phi", &tk_err_phi_);
  t->Branch("tk_err_dxy", &tk_err_dxy_);
  t->Branch("tk_err_dz", &tk_err_dz_);
  t->Branch("tk_chi2dof", &tk_chi2dof_);
  t->Branch("tk_hp", &tk_hp_);
  t->Branch("tk_minhit", &tk_minhit_);
  t->Branch("tk_maxhit", &tk_maxhit_);
  t->Branch("tk_maxpxhit", &tk_maxpxhit_);
}

void TrackingTree::read_from_tree(TTree* t) {
  t->SetBranchAddress("run", &run_);
  t->SetBranchAddress("lumi", &lumi_);
  t->SetBranchAddress("event", &event_);
  t->SetBranchAddress("npu", &npu_);
  t->SetBranchAddress("bs_x", &bs_x_);
  t->SetBranchAddress("bs_y", &bs_y_);
  t->SetBranchAddress("bs_z", &bs_z_);
  t->SetBranchAddress("bs_sigmaz", &bs_sigmaz_);
  t->SetBranchAddress("bs_dxdz", &bs_dxdz_);
  t->SetBranchAddress("bs_dydz", &bs_dydz_);
  t->SetBranchAddress("bs_width", &bs_width_);
  t->SetBranchAddress("bs_err_x", &bs_err_x_);
  t->SetBranchAddress("bs_err_y", &bs_err_y_);
  t->SetBranchAddress("bs_err_z", &bs_err_z_);
  t->SetBranchAddress("bs_err_sigmaz", &bs_err_sigmaz_);
  t->SetBranchAddress("bs_err_dxdz", &bs_err_dxdz_);
  t->SetBranchAddress("bs_err_dydz", &bs_err_dydz_);
  t->SetBranchAddress("bs_err_width", &bs_err_width_);
  t->SetBranchAddress("pv_x", &p_pv_x_);
  t->SetBranchAddress("pv_y", &p_pv_y_);
  t->SetBranchAddress("pv_z", &p_pv_z_);
  t->SetBranchAddress("pv_chi2dof", &p_pv_chi2dof_);
  t->SetBranchAddress("pv_ndof", &p_pv_ndof_);
  t->SetBranchAddress("pv_score", &p_pv_score_);
  t->SetBranchAddress("pv_cxx", &p_pv_cxx_);
  t->SetBranchAddress("pv_cxy", &p_pv_cxy_);
  t->SetBranchAddress("pv_cxz", &p_pv_cxz_);
  t->SetBranchAddress("pv_cyy", &p_pv_cyy_);
  t->SetBranchAddress("pv_cyz", &p_pv_cyz_);
  t->SetBranchAddress("pv_czz", &p_pv_czz_);
  t->SetBranchAddress("tk_qpt", &p_tk_qpt_);
  t->SetBranchAddress("tk_eta", &p_tk_eta_);
  t->SetBranchAddress("tk_phi", &p_tk_phi_);
  t->SetBranchAddress("tk_dxybs", &p_tk_dxybs_);
  t->SetBranchAddress("tk_dxypv", &p_tk_dxypv_);
  t->SetBranchAddress("tk_dzpv", &p_tk_dzpv_);
  t->SetBranchAddress("tk_vx", &p_tk_vx_);
  t->SetBranchAddress("tk_vy", &p_tk_vy_);
  t->SetBranchAddress("tk_vz", &p_tk_vz_);
  t->SetBranchAddress("tk_err_pt", &p_tk_err_pt_);
  t->SetBranchAddress("tk_err_eta", &p_tk_err_eta_);
  t->SetBranchAddress("tk_err_phi", &p_tk_err_phi_);
  t->SetBranchAddress("tk_err_dxy", &p_tk_err_dxy_);
  t->SetBranchAddress("tk_err_dz", &p_tk_err_dz_);
  t->SetBranchAddress("tk_chi2dof", &p_tk_chi2dof_);
  t->SetBranchAddress("tk_hp", &p_tk_hp_);
  t->SetBranchAddress("tk_minhit", &p_tk_minhit_);
  t->SetBranchAddress("tk_maxhit", &p_tk_maxhit_);
  t->SetBranchAddress("tk_maxpxhit", &p_tk_maxpxhit_);
}
