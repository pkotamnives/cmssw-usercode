from JMTucker.Tools.ROOTTools import *
from array import array
ROOT.TH1.AddDirectory(0)

year = '2016'

mode = 'vary_eff'
#mode = 'vary_dphi'

set_style()
ROOT.gStyle.SetOptFit(0)
ps = plot_saver('../plots/bkgest/v14/compare_dvvc_%s_%s' % (mode, year), size=(700,700), root=False, log=False)

fn1 = ['2v_from_jets_%s_3track_default_v14.root' % year, '2v_from_jets_%s_3track_%s_v14.root' % (year, mode)]
fn2 = ['2v_from_jets_%s_4track_default_v14.root' % year, '2v_from_jets_%s_4track_%s_v14.root' % (year, mode)]
fn3 = ['2v_from_jets_%s_5track_default_v14.root' % year, '2v_from_jets_%s_5track_%s_v14.root' % (year, mode)]

if mode == 'vary_eff':
    ls = ['vertex pair survival efficiency', 'vertex survival efficiency']

if mode == 'vary_dphi':
    ls = ['|#Delta#phi| from 3-track #Delta#phi_{JJ}', 'uniform |#Delta#phi|']

fns = [fn1, fn2, fn3]
ntk = ['3-track', '4-track', '5-track']

n2v = [934., 7., 1.]
ebin1 = [0.0026, 0.0062, 0.0124]
ebin2 = [0.0023, 0.0078, 0.0328]
ebin3 = [0.0060, 0.0231, 0.1019]

if year == '2015':
    n2v = [44., 1., 1.]
    ebin1 = [0.0118, 0.0291, 0.0567]
    ebin2 = [0.0106, 0.0371, 0.1513]
    ebin3 = [0.0271, 0.1100, 0.4761]

if year == '2015p6':
    n2v = [978., 7., 1.]
    ebin1 = [0.0025, 0.0062, 0.0123]
    ebin2 = [0.0023, 0.0076, 0.0323]
    ebin3 = [0.0058, 0.0229, 0.1019]

colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange, ROOT.kViolet, ROOT.kPink+1]

x = []
ex = []
y1 = []
ey1 = []
y2 = []
ey2 = []
y3 = []
ey3 = []
for i in range(3):
    print ntk[i]
    l1 = ROOT.TLegend(0.50,0.70,0.85,0.85)
    hh = ROOT.TFile(fns[i][0]).Get('h_2v_dvv')
    hh.SetTitle('%s;d_{VV} (cm);events' % ntk[i])
    hh.SetStats(0)
    hh.SetLineColor(ROOT.kBlack)
    hh.SetLineWidth(5)
    if hh.Integral() > 0:
        hh.Scale(n2v[i]/hh.Integral())
    hh.Draw()
    l1.AddEntry(hh, 'simulated events')

    hs = []
    for j in range(len(ls)):
        h = ROOT.TFile(fns[i][j]).Get('h_c1v_dvv')
        h.SetStats(0)
        h.SetLineColor(colors[j])
        h.SetLineWidth(2)
        h.Scale(n2v[i]/h.Integral())
        h.Draw('hist e sames')
        hs.append(h)
        l1.AddEntry(h, ls[j])

    l1.SetFillColor(0)
    l1.Draw()
    ps.save(ntk[i])

    if i == 2:
        hs = []
        l1 = ROOT.TLegend(0.50,0.75,0.85,0.85)
        for j in range(len(ls)):
            h = ROOT.TFile(fns[i][j]).Get('h_c1v_dvv')
            h.SetStats(0)
            h.SetLineColor(colors[j])
            h.SetLineWidth(2)
            h.Scale(n2v[i]/h.Integral())
            if j == 0:
                h.SetTitle(';d_{VV}^{C} (cm);Events')
                h.Draw('hist e')
            elif j == 1:
                h.Draw('hist e sames')
            hs.append(h)
            if j > 1:
                continue
            l1.AddEntry(h, ls[j])
        l1.SetFillColor(0)
        l1.Draw()
        ps.save('compare_dvvc_%s' % mode)

    l2 = ROOT.TLegend(0.15,0.70,0.50,0.85)
    h2 = ROOT.TFile(fns[i][0]).Get('h_2v_absdphivv')
    h2.SetTitle('%s;|#Delta#phi_{VV}|;events' % ntk[i])
    h2.SetStats(0)
    h2.SetLineColor(ROOT.kBlack)
    h2.SetLineWidth(5)
    if h2.Integral() > 0:
        h2.Scale(n2v[i]/h2.Integral())
    h2.SetMinimum(0)
    h2.Draw()
    l2.AddEntry(h2, 'simulated events')

    h2s = []
    for j in range(len(ls)):
        h = ROOT.TFile(fns[i][j]).Get('h_c1v_absdphivv')
        h.SetStats(0)
        h.SetLineColor(colors[j])
        h.SetLineWidth(2)
        h.Scale(n2v[i]/h.Integral())
        h.Draw('hist e sames')
        h2s.append(h)
        l2.AddEntry(h, ls[j])

        chi2 = 0
        for k in range(1,h.GetNbinsX()+1):
            if (h2.GetBinError(k) > 0):
                chi2 += (h.GetBinContent(k)-h2.GetBinContent(k))**2 / h2.GetBinError(k)**2
        print '%35s: deltaphi chi2/ndf = %f' % (ls[j], chi2/(h.GetNbinsX()-1))

    l2.SetFillColor(0)
    l2.Draw()
    ps.save('%s_dphi'%ntk[i])

    if i == 2:
        h2s = []
        l2 = ROOT.TLegend(0.15,0.75,0.50,0.85)
        for j in range(len(ls)):
            h = ROOT.TFile(fns[i][j]).Get('h_c1v_absdphivv')
            h.SetStats(0)
            h.SetLineColor(colors[j])
            h.SetLineWidth(2)
            h.Scale(n2v[i]/h.Integral())
            if j == 0:
                h.SetTitle(';|#Delta#phi_{VV}|;Events')
                h.SetMinimum(0)
                h.Draw('hist e')
            elif j == 1:
                h.Draw('hist e sames')
            h2s.append(h)
            if j > 1:
                continue
            l2.AddEntry(h, ls[j])
        l2.SetFillColor(0)
        l2.Draw()
        ps.save('compare_dphi_%s' % mode)

    es1 = ROOT.Double(0)
    s1 = hh.IntegralAndError(1,4,es1)
    es2 = ROOT.Double(0)
    s2 = hh.IntegralAndError(5,7,es2)
    es3 = ROOT.Double(0)
    s3 = hh.IntegralAndError(8,40,es3)

    c1 = hs[0].Integral(1,4)
    ec1 = ebin1[i] * c1
    c2 = hs[0].Integral(5,7)
    ec2 = ebin2[i] * c2
    c3 = hs[0].Integral(8,40)
    ec3 = ebin3[i] * c3

    v1 = hs[1].Integral(1,4)
    ev1 = ebin1[i] * v1
    v2 = hs[1].Integral(5,7)
    ev2 = ebin2[i] * v2
    v3 = hs[1].Integral(8,40)
    ev3 = ebin3[i] * v3

    r1 = v1/c1
    er1 = (v1/c1) * ((ev1/v1)**2 + (ec1/c1)**2)**0.5
    r2 = v2/c2
    er2 = (v2/c2) * ((ev2/v2)**2 + (ec2/c2)**2)**0.5
    r3 = v3/c3
    er3 = (v3/c3) * ((ev3/v3)**2 + (ec3/c3)**2)**0.5

    if mode == 'vary_dphi':
        if year == '2015':
            if i == 0:
                er1 *= 0.122017
                er2 *= 0.217034
                er3 *= 0.178617
            if i == 1:
                er1 *= 0.109714
                er2 *= 0.219117
                er3 *= 0.189308
            if i == 2:
                er1 *= 0.076194
                er2 *= 0.228316
                er3 *= 0.142351
        if year == '2016':
            if i == 0:
                er1 *= 0.114920
                er2 *= 0.213661
                er3 *= 0.170409
            if i == 1:
                er1 *= 0.105188
                er2 *= 0.220908
                er3 *= 0.187441
            if i == 2:
                er1 *= 0.073343
                er2 *= 0.219650
                er3 *= 0.160530
        if year == '2015p6':
            if i == 0:
                er1 *= 0.114961
                er2 *= 0.213431
                er3 *= 0.170297
            if i == 1:
                er1 *= 0.104845
                er2 *= 0.219375
                er3 *= 0.185960
            if i == 2:
                er1 *= 0.071913
                er2 *= 0.215773
                er3 *= 0.157685
    if mode == 'vary_eff':
        er2 *= (abs(r2-1))**0.5 / (1+r2)**0.5
        er3 *= (abs(r3-1))**0.5 / (1+r3)**0.5

    print
    print '    simulated events: 0-400 um: %6.2f +/- %5.2f, 400-700 um: %6.2f +/- %5.2f, 700-40000 um: %6.2f +/- %5.2f' % (s1, es1, s2, es2, s3, es3)
    print 'default construction: 0-400 um: %6.2f +/- %5.2f, 400-700 um: %6.2f +/- %5.2f, 700-40000 um: %6.2f +/- %5.2f' % (c1, ec1, c2, ec2, c3, ec3)
    print '           variation: 0-400 um: %6.2f +/- %5.2f, 400-700 um: %6.2f +/- %5.2f, 700-40000 um: %6.2f +/- %5.2f' % (v1, ev1, v2, ev2, v3, ev3)
    print ' variation / default: 0-400 um: %6.2f +/- %5.2f, 400-700 um: %6.2f +/- %5.2f, 700-40000 um: %6.2f +/- %5.2f' % (r1, er1, r2, er2, r3, er3)
    print
    print

    x.append(i-2)
    ex.append(0)
    y1.append(r1)
    ey1.append(er1)
    y2.append(r2)
    ey2.append(er2)
    y3.append(r3)
    ey3.append(er3)

bins = ['bin1', 'bin2', 'bin3']
dvvc = ['d_{VV}^{C} < 400 #mum', '400 #mum < d_{VV}^{C} < 700 #mum', 'd_{VV}^{C} > 700 #mum']
ys = [y1, y2, y3]
eys = [ey1, ey2, ey3]
for i in range(3):
    g = ROOT.TGraphErrors(len(x), array('d',x), array('d',ys[i]), array('d',ex), array('d',eys[i]))
    g.SetMarkerStyle(21)
    g.SetTitle('variation / default (%s);3-track%12s4-track%12s5-or-more-track%2s' % (dvvc[i], '','',''))
    g.GetXaxis().SetLimits(-3,1)
    g.GetXaxis().SetLabelSize(0)
    g.GetXaxis().SetTitleOffset(0.5)
    g.GetYaxis().SetRangeUser(0,2)
    g.Draw('AP')

    line = ROOT.TLine(-3,1,1,1)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()

    if mode == 'vary_dphi':
        t5 = ROOT.TLatex()
        t5.SetTextFont(42)
        t5.SetTextSize(0.04)
        t5.DrawLatex(0, ys[i][2] - eys[i][2] - t5.GetTextSize(), '%.2f #pm %.2f' % (ys[i][2], eys[i][2]))

        t = ROOT.TLatex()
        t.SetTextFont(42)
        t.SetTextSize(0.04)
        t.DrawLatex(-2.5, 0.2, '#splitline{difference of 5-or-more-track ratio from 1:}{%.2f #pm %.2f}' % (abs(ys[i][2] - 1), eys[i][2]))

    else:
        r = g.Fit('pol1','S','',-2,-1)
        g5 = ROOT.TGraphErrors(1, array('d',[-0.1]), array('d',[r.Value(0)]), array('d',[0]), array('d',[r.ParError(0)]))
        g5.SetLineColor(ROOT.kRed)
        g5.SetMarkerColor(ROOT.kRed)
        g5.SetMarkerStyle(21)
        g5.Draw('P')

        t5 = ROOT.TLatex()
        t5.SetTextFont(42)
        t5.SetTextSize(0.04)
        t5.SetTextColor(ROOT.kRed)
        t5.DrawLatex(-0.1, r.Value(0) - r.ParError(0) - t5.GetTextSize(), '%.2f #pm %.2f' % (r.Value(0), r.ParError(0)))

        t3 = ROOT.TLatex()
        t3.SetTextFont(42)
        t3.SetTextSize(0.04)
        t3.DrawLatex(-2, ys[i][0] - eys[i][0] - t3.GetTextSize(), '%.2f #pm %.2f' % (ys[i][0], eys[i][0]))

        t = ROOT.TLatex()
        t.SetTextFont(42)
        t.SetTextSize(0.04)
        t.DrawLatex(-2.5, 0.2, '#splitline{difference from 3-track to 5-or-more-track:}{%.2f #pm %.2f}' % (abs(r.Value(0) - ys[i][0]), (r.ParError(0)**2 + eys[i][0]**2)**0.5))

    ps.save('ratio_%s_%s' % (bins[i], mode))
