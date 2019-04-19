from JMTucker.Tools.ROOTTools import *

ntk = 3
f = ROOT.TFile('output_btags_vs_bquarks_MiniTreeV23m_ntk%s/background.root' % ntk)

def njets(hname):
  h = f.Get(hname)
  return h.Integral() * h.GetMean()

def btag_eff_per_jet(nvtx, jet_flavor, bdisc):
  num = njets('h_%dv_n%sjets_%s_btag' % (nvtx, jet_flavor, bdisc))
  den = njets('h_%dv_n%sjets' % (nvtx, jet_flavor))
  return num/den

def scale_factor(nvtx, jet_flavor, bdisc):
  h = f.Get('h_%dv_scalefactor_%s_%s_btag' % (nvtx, jet_flavor, bdisc))
  return h.GetMean()

def btag_eff_per_event_from_btag_eff_per_jet(nvtx, event_flavor, effb, effc, effl):
  num = 0
  den = 0
  for nc in range(0, 40):
    h_nbnl = f.Get('h_%dv_%dcjets_nljets_vs_nbjets' % (nvtx, nc))
    for nb in range(0, 1) if event_flavor == 'nobjets' else range(1, h_nbnl.GetNbinsX()):
      for nl in range(0, h_nbnl.GetNbinsY()):
        n_nbnl = h_nbnl.GetBinContent(nb+1, nl+1)
        num += n_nbnl * (1 - (1-effb)**nb * (1-effc)**nc * (1-effl)**nl)
        den += n_nbnl
  return num/den

def btag_eff_per_event(nvtx, event_flavor, bdisc):
  h = f.Get('h_%s_%dv_1%s_btag_flavor_code' % (event_flavor, nvtx, bdisc))
  return h.GetBinContent(2) / h.Integral()

print '%10s%90s%26s%26s%26s' % ('', 'per-jet', 'per-event from per-jet', 'per-event from per-jet*SF', 'per-event')
fmt = '%10s' + '%10s'*9 + '%13s'*6
print fmt % ('', 'effb', 'SFb', 'effb*SFb', 'effc', 'SFc', 'effc*SFc', 'effl', 'SFl', 'effl*SFl', 'btag eff', 'fake rate', 'btag eff', 'fake rate', 'btag eff', 'fake rate')
for nvtx in [1, 2]:
  print '%d-track %d-vertex' % (ntk, nvtx)
  for bdisc in ['loose', 'medium', 'tight']:
    effb, sfb = btag_eff_per_jet(nvtx, 'b', bdisc), scale_factor(nvtx, 'b', bdisc)
    effc, sfc = btag_eff_per_jet(nvtx, 'c', bdisc), scale_factor(nvtx, 'c', bdisc)
    effl, sfl = btag_eff_per_jet(nvtx, 'l', bdisc), scale_factor(nvtx, 'l', bdisc)
    print fmt % (bdisc, '%.3f' % effb, '%.3f' % sfb, '%.3f' % (effb*sfb),
                        '%.3f' % effc, '%.3f' % sfc, '%.3f' % (effc*sfc),
                        '%.3f' % effl, '%.3f' % sfl, '%.3f' % (effl*sfl),
                        '%.3f' % btag_eff_per_event_from_btag_eff_per_jet(nvtx, 'bjets', effb, effc, effl),
                        '%.3f' % btag_eff_per_event_from_btag_eff_per_jet(nvtx, 'nobjets', effb, effc, effl),
                        '%.3f' % btag_eff_per_event_from_btag_eff_per_jet(nvtx, 'bjets', effb*sfb, effc*sfc, effl*sfl),
                        '%.3f' % btag_eff_per_event_from_btag_eff_per_jet(nvtx, 'nobjets', effb*sfb, effc*sfc, effl*sfl),
                        '%.3f' % btag_eff_per_event(nvtx, 'bquarks', bdisc),
                        '%.3f' % btag_eff_per_event(nvtx, 'nobquarks', bdisc))