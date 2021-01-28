#! /usr/bin/env python

import ROOT
import sys,os
import math
import numpy as np

#from DataFormats.FWLite import Events, Handle

from JMTucker.Tools.ROOTTools import *
cmssw_setup()

# FIXME you can replace this with the usual stuff for putting plots into our publicweb areas and generating the html
outputdir = "/eos/user/p/pekotamn/www/Evernote_results_err_of_sharedjet_studies
outputdir += "/" # in case we forget it...
os.system("mkdir -p "+outputdir)






# Create histograms, etc.
ROOT.gROOT.SetBatch() # don't pop up canvases

# create an output root file
outfile = ROOT.TFile(outputdir+"out.root", "RECREATE")

# Define histograms here (obviously this one doesn't matter for you, but I stole it from some other code of mine)

n = 4
x = [1,2,3,4]
y = [39.5,69.3,79.0,88.1]
ex = [0,0,0,0]
ey = [1.2,0.9,0.9,0.4]

x = np.asaray(x)
y = np.asaray(x)
ex = np.asaray(ex)
ey = np.asaray(ey)

mg = ROOT.TMultiGraph()

g_blue = ROOT.TGraphErrors(4,x,y,ex,ey)
g_blue.SetLineColor( ROOT.kBlue )
g_blue.SetLineWidth( 3 )
mg.Add(g_blue)

mg.SetTitle("Dijets: Validation of pT variables by varying shared-track ratios; # major vtx's shared tracks/ # of minor vtx's shared tracks; signal efficiency (%)")


# make a canvas, draw, and save it


c0 = ROOT.TCanvas()
mg.Draw("ACP")
c0.Print (outputdir+"plot_dd.png")
c0.Print (outputdir+"plot_dd.root")

