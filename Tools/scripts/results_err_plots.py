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

n = 4
x = [3.0,5.0,7.0,9.0]
y1 = [90.88,89.85,88.87,93.53]
ex = [0.0,0.0,0.0,0.0]
ey1 = [0.26,0.27,0.28,0.22]

x = np.array(x)
y1 = np.array(y1)
ex = np.array(ex)
ey1 = np.array(ey1)


y2 = [90.9,90.81,89.41,93.42]
ey2 = [0.23,0.23,0.24,0.20]
y2 = np.array(y2)
ey2 = np.array(ey2)

"""
y3 = [44.7,65.1,71.4,79.4]
ey3 = [1.3,0.9,1.0,0.5]
y3 = np.array(y3)
ey3 = np.array(ey3)

y4 = [31.5,76.8,83.6,93.2]
ey4 = [1.2,0.8,0.8,0.3]
y4 = np.array(y4)
ey4 = np.array(ey4)
"""

mg = ROOT.TMultiGraph()

g1 = ROOT.TGraphErrors(4,x,y1,ex,ey1)
g1.SetLineColor( ROOT.kBlue )
g1.SetTitle( "multijets" )
g1.SetLineWidth( 2 )
#g1.SetMarkerStyle(21)
g2 = ROOT.TGraphErrors(4,x,y2,ex,ey2)
g2.SetLineColor( ROOT.kRed )
g2.SetTitle( "dijets" )
g2.SetLineWidth( 2 )
#g2.SetMarkerStyle(21)
"""
g3 = ROOT.TGraphErrors(4,x,y3,ex,ey3)
g3.SetLineColor( ROOT.kOrange )
g3.SetTitle( "ratio average shared-track pT" )
g3.SetLineWidth( 2 )
#g3.SetMarkerStyle(21)
g4 = ROOT.TGraphErrors(4,x,y4,ex,ey4)
g4.SetLineColor( ROOT.kGreen )
g4.SetTitle( "ratio sum shared-track pT" )
g4.SetLineWidth( 2 )
"""
#g4.SetMarkerStyle(21)
mg.Add(g1)
mg.Add(g2)
#mg.Add(g3)
#mg.Add(g4)


mg.SetTitle("Signal efficiency after removing shared tracks and apply >=5trks/vtx; Method #; signal efficiency by pT variables (%)")


# make a canvas, draw, and save it


c0 = ROOT.TCanvas()
mg.Draw("ALP")
c0.BuildLegend(0.45,0.21,0.75,0.51)
c0.Print (outputdir+"plot_eff.png")
c0.Print (outputdir+"plot_eff.root")

