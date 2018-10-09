from JMTucker.Tools.ROOTTools import *
ROOT.TH1.AddDirectory(0)

is_mc = True
only_10pc = False
year = '2017'

set_style()
ps = plot_saver(plot_dir('closure%s%s_%s' % ('' if is_mc else '_data', '_10pc' if only_10pc else '', year)), size=(700,700), root=False, log=False)

fns = ['2v_from_jets%s_%s_3track_default_v20mp2.root' % ('' if is_mc else '_data', year), '2v_from_jets%s_%s_7track_default_v20mp2.root' % ('' if is_mc else '_data', year), '2v_from_jets%s_%s_4track_default_v20mp2.root' % ('' if is_mc else '_data', year), '2v_from_jets%s_%s_5track_default_v20mp2.root' % ('' if is_mc else '_data', year)]
ntk = ['3-track', '4-track-3-track', '4-track', '5-track']
names = ['3-track x 3-track', '4-track x 3-track', '4-track x 4-track', '#geq 5-track x #geq 5-track']

n2v = [773., 9., 5., 1.] if year == '2017' else [991., 213., 8., 1.]

if is_mc:
    ebin1 = [0.0032, 0.0091, 0.0249, 0.0091] if year == '2017' else [1, 1, 1, 1]
    ebin2 = [0.0030, 0.0113, 0.0648, 0.0113] if year == '2017' else [1, 1, 1, 1]
    ebin3 = [0.0078, 0.0348, 0.2124, 0.0348] if year == '2017' else [1, 1, 1, 1]
elif only_10pc:
    ebin1 = [1, 1, 1, 1] if year == '2017' else [1, 1, 1, 1]
    ebin2 = [1, 1, 1, 1] if year == '2017' else [1, 1, 1, 1]
    ebin3 = [1, 1, 1, 1] if year == '2017' else [1, 1, 1, 1]
else:
    ebin1 = [1, 1, 1, 1] if year == '2017' else [1, 1, 1, 1]
    ebin2 = [1, 1, 1, 1] if year == '2017' else [1, 1, 1, 1]
    ebin3 = [1, 1, 1, 1] if year == '2017' else [1, 1, 1, 1]

def write(font, size, x, y, text):
    w = ROOT.TLatex()
    w.SetNDC()
    w.SetTextFont(font)
    w.SetTextSize(size)
    w.DrawLatex(x, y, text)
    return w

for i in range(4):
    if not is_mc and i > 2:
        h = ROOT.TFile(fns[i]).Get('h_c1v_dvv')
        h.SetTitle(';d_{VV}^{C} (cm);')
        h.SetStats(0)
        h.SetLineColor(ROOT.kRed)
        h.SetLineWidth(2)
        h.Scale(1./h.Integral())
        h.Draw('hist e')
        ps.save('%s_dvvc' % ntk[i])

        ec = ROOT.Double(0)
        c = h.IntegralAndError(1,40,ec)
        c1 = h.Integral(1,4)
        ec1 = ebin1[i] * c1
        c2 = h.Integral(5,7)
        ec2 = ebin2[i] * c2
        c3 = h.Integral(8,40)
        ec3 = ebin3[i] * c3

        print ntk[i]
        print ' constructed events: %7.2f +/- %5.2f, 0-400 um: %7.2f +/- %5.2f, 400-700 um: %6.2f +/- %5.2f, 700-40000 um: %6.2f +/- %5.2f' % (c, ec, c1, ec1, c2, ec2, c3, ec3)
        print '    dVVC normalized: %7.3f +/- %5.3f, 0-400 um: %7.3f +/- %5.3f, 400-700 um: %6.3f +/- %5.3f, 700-40000 um: %6.3f +/- %5.3f' % (c/c, ec/c, c1/c, ec1/c, c2/c, ec2/c, c3/c, ec3/c)

    hh = ROOT.TFile(fns[i]).Get('h_2v_dvv')
    hh.SetTitle(';d_{VV} (cm);Events')
    hh.SetStats(0)
    hh.SetLineColor(ROOT.kBlue)
    hh.SetLineWidth(2)
    if is_mc:
        if hh.Integral() > 0:
            hh.Scale(n2v[i]/hh.Integral())
        else:
            hh.SetMaximum(0.4)
    hh.SetMinimum(0)
    hh.Draw()

    h = ROOT.TFile(fns[i]).Get('h_c1v_dvv')
    h.SetStats(0)
    h.SetLineColor(ROOT.kRed)
    h.SetLineWidth(2)
    if is_mc:
        h.Scale(n2v[i]/h.Integral())
    else:
        if hh.Integral() > 0:
            h.Scale(hh.Integral()/h.Integral())
        else:
            h.Scale(1./h.Integral())
    h.Draw('hist e sames')

    l1 = ROOT.TLegend(0.35, 0.75, 0.85, 0.85)
    l1.AddEntry(hh, 'Simulated events' if is_mc else 'Data')
    l1.AddEntry(h, 'd_{VV}^{C}')
    l1.SetFillColor(0)
    l1.Draw()
    ps.save(ntk[i])

    es = ROOT.Double(0)
    s = hh.IntegralAndError(1,40,es)
    if s == 0:
        s = 1
    es1 = ROOT.Double(0)
    s1 = hh.IntegralAndError(1,4,es1)
    if s1 == 0:
        s1 = 1
    es2 = ROOT.Double(0)
    s2 = hh.IntegralAndError(5,7,es2)
    if s2 == 0:
        s2 = 1
    es3 = ROOT.Double(0)
    s3 = hh.IntegralAndError(8,40,es3)
    if s3 == 0:
        s3 = 1

    ec = ROOT.Double(0)
    c = h.IntegralAndError(1,40,ec)
    c1 = h.Integral(1,4)
    ec1 = ebin1[i] * c1
    c2 = h.Integral(5,7)
    ec2 = ebin2[i] * c2
    c3 = h.Integral(8,40)
    ec3 = ebin3[i] * c3

    r = c/s
    er = (c/s) * ((ec/c)**2 + (es/s)**2)**0.5
    r1 = c1/s1
    er1 = (c1/s1) * ((ec1/c1)**2 + (es1/s1)**2)**0.5
    r2 = c2/s2
    er2 = (c2/s2) * ((ec2/c2)**2 + (es2/s2)**2)**0.5
    r3 = c3/s3
    er3 = (c3/s3) * ((ec3/c3)**2 + (es3/s3)**2)**0.5

    print ntk[i]
    print '   simulated events: %7.2f +/- %5.2f, 0-400 um: %7.2f +/- %5.2f, 400-700 um: %6.2f +/- %5.2f, 700-40000 um: %6.2f +/- %5.2f' % (s, es, s1, es1, s2, es2, s3, es3)
    print ' constructed events: %7.2f +/- %5.2f, 0-400 um: %7.2f +/- %5.2f, 400-700 um: %6.2f +/- %5.2f, 700-40000 um: %6.2f +/- %5.2f' % (c, ec, c1, ec1, c2, ec2, c3, ec3)
    print '     dVV normalized: %7.3f +/- %5.3f, 0-400 um: %7.3f +/- %5.3f, 400-700 um: %6.3f +/- %5.3f, 700-40000 um: %6.3f +/- %5.3f' % (s/s, es/s, s1/s, es1/s, s2/s, es2/s, s3/s, es3/s)
    print '    dVVC normalized: %7.3f +/- %5.3f, 0-400 um: %7.3f +/- %5.3f, 400-700 um: %6.3f +/- %5.3f, 700-40000 um: %6.3f +/- %5.3f' % (c/c, ec/c, c1/c, ec1/c, c2/c, ec2/c, c3/c, ec3/c)
    print '   ratio dVVC / dVV: %7.2f +/- %5.2f, 0-400 um: %7.2f +/- %5.2f, 400-700 um: %6.2f +/- %5.2f, 700-40000 um: %6.2f +/- %5.2f' % (r, er, r1, er1, r2, er2, r3, er3)

    hh = ROOT.TFile(fns[i]).Get('h_2v_absdphivv')
    hh.SetTitle(';|#Delta#phi_{VV}|;Events')
    hh.SetStats(0)
    hh.SetLineColor(ROOT.kBlue)
    hh.SetLineWidth(2)
    if is_mc:
        if hh.Integral() > 0:
            hh.Scale(n2v[i]/hh.Integral())
        else:
            hh.SetMaximum(0.4)
    hh.SetMinimum(0)
    hh.Draw()

    h = ROOT.TFile(fns[i]).Get('h_c1v_absdphivv')
    h.SetStats(0)
    h.SetLineColor(ROOT.kRed)
    h.SetLineWidth(2)
    if is_mc:
        h.Scale(n2v[i]/h.Integral())
    else:
        if hh.Integral() > 0:
            h.Scale(hh.Integral()/h.Integral())
        else:
            h.Scale(1./h.Integral())
    h.Draw('hist e sames')

    l1 = ROOT.TLegend(0.15, 0.75, 0.65, 0.85)
    l1.AddEntry(hh, 'Simulated events' if is_mc else 'Data')
    l1.AddEntry(h, 'Construction')
    l1.SetFillColor(0)
    l1.Draw()

    if not is_mc and not only_10pc and year == '2015p6':
        write(42, 0.040, 0.150, 0.700, names[i])
        write(61, 0.050, 0.098, 0.913, 'CMS')
        write(52, 0.035, 0.200, 0.913, 'Preliminary')
        write(42, 0.050, 0.560, 0.913, '38.5 fb^{-1} (13 TeV)')

    ps.save('%s_dphi' % ntk[i])
