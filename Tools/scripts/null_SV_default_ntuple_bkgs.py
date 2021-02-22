#! /usr/bin/env python

import ROOT
import sys,os
import math
import numpy as np

from DataFormats.FWLite import Events, Handle

from JMTucker.Tools.ROOTTools import *
cmssw_setup()

# FIXME you can replace this with the usual stuff for putting plots into our publicweb areas and generating the html
outputdir = "/eos/user/p/pekotamn/www/Null_SV_tbs_1mm_1.6TeV_name_1000n_ntuple"
outputdir += "/" # in case we forget it...
os.system("mkdir -p "+outputdir)


events_ntuple1 = Events (sys.argv[1])

# create handle and labels outside of loop
vertices_handle1  = Handle ("std::vector<MFVVertexAux>")
vertices_label1 = ("mfvVerticesAux")

# FIXME should be able to do the same for our mfv events, just didn't do anything with it yet
# event_handle1  = Handle ("MFVEvent")
# event_label1 = ("mfvEvent")

# Create histograms, etc.
ROOT.gROOT.SetBatch() # don't pop up canvases

# create an output root file
outfile = ROOT.TFile(outputdir+"out.root", "RECREATE")

# Define histograms here (obviously this one doesn't matter for you, but I stole it from some other code of mine)
h_qual_nsv_event = ROOT.TH1F ("h_qual_nsv_event", ";# of >=5trk-SVs/event", 50, 0, 50)
h_dBV =  ROOT.TH1F("h_dBV", ";dist2d(beamspot, >=5trk-SV in a hemisphere) (cm);arb. units", 100, 0, 0.02)


nevents_processed = 0
nevents_fiducial_cuts = 0
nevents_nsv01_fiducial_cuts = 0

for event1 in events_ntuple1 :

    if nevents_processed % 1000 == 0 :
        print "Processing event #%s" % (nevents_processed)

    run_lumi_event_number1 = (event1.eventAuxiliary().id().run(), event1.eventAuxiliary().id().luminosityBlock(), event1.eventAuxiliary().id().event())

    # associate the handle with the label and the event
    event1.getByLabel (vertices_label1, vertices_handle1)

    # get the product: these are the objects that we'll access!
    vertices_from_ntuple1 = vertices_handle1.product()

    mfv_event1_handle1 = Handle ("MFVEvent")
    mfv_event1_label1 = ("mfvEvent")
    event1.getByLabel (mfv_event1_label1, mfv_event1_handle1)
    mevent = mfv_event1_handle1.product()

   

    nevents_processed += 1
    if True :
         
         beamspot1 = Handle ("reco::BeamSpot")
         beamspot_label1 = ("offlineBeamSpot")
         event1.getByLabel (beamspot_label1, beamspot1)
         beamspot = beamspot1.product() 
         bsx = beamspot.position().x()
         bsy = beamspot.position().y()
         bsz = beamspot.position().z()
         print "the beam spot's r is #%s" % ( np.linalg.norm(np.array([bsx,bsy]))   )
         
              
         if True: # apply fiducial cuts
        #if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0], mevent.gen_lsp_phi[1])) > 2.7: # no apply fiducial cuts
            
               nevents_fiducial_cuts += 1
               n_vertex_seed_tracks = mevent.n_vertex_seed_tracks()
               qual_nsv = 0

               

               for vtx_ntuple1 in vertices_from_ntuple1 : # loop raw verices

                   dBV_vtx_ntuple1 = np.array([vtx_ntuple1.x - bsx,vtx_ntuple1.y - bsy])
                   
                   if vtx_ntuple1.ntracks()>=5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09 and np.linalg.norm(dBV_vtx_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                          qual_nsv += 1
                          h_dBV.Fill(dBV_vtx_ntuple1) 

               
                   
               if qual_nsv < 2:
                   nevents_nsv01_fiducial_cuts += 1

               h_qual_nsv_event.Fill(qual_nsv)
                   

    else:
        break



print "Total processed event #%s" % (nevents_processed-1)
print "Total processed event (fiducial) #%s" % (nevents_fiducial_cuts)
print "Total nsv<2 events (ficudial) #%s" % (nevents_nsv01_fiducial_cuts)

# make a canvas, draw, and save it



c01 = ROOT.TCanvas()                                                      
h_qual_nsv_event.Draw("colz")
c01.Print (outputdir+"h_qual_nsv_event.png")
c01.Print (outputdir+"h_qual_nsv_event.root")

c1 = ROOT.TCanvas()                                                      
h_dBV.Draw("colz")
c1.Print (outputdir+"h_dBV.png")
c1.Print (outputdir+"h_dBV.root")





