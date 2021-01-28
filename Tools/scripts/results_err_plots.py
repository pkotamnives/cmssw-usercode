#! /usr/bin/env python

import ROOT
import sys,os
import math
import numpy as np

#from DataFormats.FWLite import Events, Handle

from JMTucker.Tools.ROOTTools import *
cmssw_setup()

# FIXME you can replace this with the usual stuff for putting plots into our publicweb areas and generating the html
outputdir = "/eos/user/p/pekotamn/www/Evernote_results_err_of_sharedjet_studies"
outputdir += "/" # in case we forget it...
os.system("mkdir -p "+outputdir)






# Create histograms, etc.
ROOT.gROOT.SetBatch() # don't pop up canvases

# create an output root file
outfile = ROOT.TFile(outputdir+"out.root", "RECREATE")

# Define histograms here (obviously this one doesn't matter for you, but I stole it from some other code of mine)

n = 3
x = [300.0,1000.0,10000.0]
y1 = [91.16,96.98,99.23]
ex = [0.0,0.0,0.0]
ey1 = [0.23,0.08,0.03]

x = np.array(x)
y1 = np.array(y1)
ex = np.array(ex)
ey1 = np.array(ey1)


y2 = [88.75,96.04,99.01]
ey2 = [0.28,0.1,0.04]
y2 = np.array(y2)
ey2 = np.array(ey2)


y3 = [93.73,98.07,99.61]
ey3 = [0.2,0.07,0.02]
y3 = np.array(y3)
ey3 = np.array(ey3)

y4 = [91.66,97.16,99.4]
ey4 = [0.25,0.09,0.03]
y4 = np.array(y4)
ey4 = np.array(ey4)


mg = ROOT.TMultiGraph()

g1 = ROOT.TGraphErrors(4,x,y1,ex,ey1)
g1.SetLineColor( ROOT.kRed )
g1.SetTitle( "dijets w/ sum shared-track pT" )
g1.SetLineWidth( 2 )
#g1.SetMarkerStyle(21)
g2 = ROOT.TGraphErrors(4,x,y2,ex,ey2)
g2.SetLineColor( ROOT.kBlue )
g2.SetTitle( "multijets w/ sum shared-track pT" )
g2.SetLineWidth( 2 )
#g2.SetMarkerStyle(21)

g3 = ROOT.TGraphErrors(4,x,y3,ex,ey3)
g3.SetLineColor( ROOT.kPink+9 )
g3.SetTitle( "dijets w/ ratio sum shared-track pT" )
g3.SetLineWidth( 2 )
#g3.SetMarkerStyle(21)
g4 = ROOT.TGraphErrors(4,x,y4,ex,ey4)
g4.SetLineColor( ROOT.kAzure+1 )
g4.SetTitle( "multijets w/ ratio sum shared-track pT" )
g4.SetLineWidth( 2 )

#g4.SetMarkerStyle(21)
mg.Add(g1)
mg.Add(g2)
mg.Add(g3)
mg.Add(g4)


mg.SetTitle("1.6TeV-mass signal efficiency after removing shared tracks and apply >=5trks/vtx; lifetimes (microns); signal efficiency from full-selection events (%)")


# make a canvas, draw, and save it


c0 = ROOT.TCanvas()
mg.Draw("ALP")
c0.BuildLegend(0.45,0.21,0.75,0.51)
c0.Print (outputdir+"plot_lifetimes.png")
c0.Print (outputdir+"plot_lifetimes.root")

