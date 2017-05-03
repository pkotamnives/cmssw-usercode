#!/usr/bin/env python

# why doesn't comparehists work on 2v_from_jets output?

import sys, os
from JMTucker.Tools.ROOTTools import *
set_style()

output_dir = sys.argv[1]
fns = [fn for fn in sys.argv[2:] if os.path.isfile(fn)]

if not os.path.isdir(os.path.dirname(output_dir)):
    output_dir = plot_dir(output_dir)

ps = plot_saver(output_dir, size=(600,600))

fs = [ROOT.TFile(fn) for fn in fns]
hs = [f.Get('h_c1v_dvv') for f in fs]
colors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kMagenta]
while len(colors) < len(fns):
    colors.append(colors[-1]+1)

for i,h in enumerate(hs):
    h.SetLineColor(colors[i])
    h.Scale(1/h.Integral(0,100000000))
    h.SetStats(0)
    h.SetLineWidth(2)

    if i == 0:
        h.Draw('hist e')
    else:
        h.Draw('hist e same')

    print '        ', '%18s' % '0-400 um', '%18s' % '400-700 um', '%18s' % '>700um'
    for fn, h in zip(fns, hs):
        name = os.path.basename(fn).replace('.root', '')
        print '%8s:' % name, '%7.4f +- %7.4f' % get_integral(hdef, 0, 0.04), '%7.4f +- %7.4f' % get_integral(hdef, 0.04, 0.07), '%7.4f +- %7.4f' % get_integral(hdef, 0.07, 4)
    ps.save('ntracks%i' % ntk)
