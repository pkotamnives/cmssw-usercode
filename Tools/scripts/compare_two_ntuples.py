#! /usr/bin/env python

import ROOT
import sys,os
import math
import numpy as np

from DataFormats.FWLite import Events, Handle

from JMTucker.Tools.ROOTTools import *
cmssw_setup()

# FIXME you can replace this with the usual stuff for putting plots into our publicweb areas and generating the html
outputdir = "/eos/user/p/pekotamn/www/Compare_10mm_3.0TeV_4sigma_default_ntuples"
outputdir += "/" # in case we forget it...
os.system("mkdir -p "+outputdir)

events_ntuple1 = Events (sys.argv[1])
events_ntuple2 = Events (sys.argv[2])

# create handle and labels outside of loop
vertices_handle1  = Handle ("std::vector<MFVVertexAux>")
vertices_label1 = ("mfvVerticesAux")

vertices_handle2  = Handle ("std::vector<MFVVertexAux>")
vertices_label2 = ("mfvVerticesAux")

# FIXME should be able to do the same for our mfv events, just didn't do anything with it yet
# event_handle1  = Handle ("MFVEvent")
# event_label1 = ("mfvEvent")

# Create histograms, etc.
ROOT.gROOT.SetBatch() # don't pop up canvases

# create an output root file
outfile = ROOT.TFile(outputdir+"out.root", "RECREATE")

# Define histograms here (obviously this one doesn't matter for you, but I stole it from some other code of mine)
h_2D_nsv = ROOT.TH2F ("h_2D_nsv", ";# of SVs (4sigma);# of SVs in (default)", 15, 0, 15, 15, 0, 15)

nevents_processed = 0
nevents_presel_nsv01_ntuple1 = 0
nevents_presel_nsv01_ntuple2 = 0

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
    mevent1 = mfv_event1_handle1.product()

 
    for event2 in events_ntuple2 :
        
        run_lumi_event_number2 = (event2.eventAuxiliary().id().run(), event2.eventAuxiliary().id().luminosityBlock(), event2.eventAuxiliary().id().event())

        # FIXME don't forget, you'll need to think about how you want to handle cases where an event is present in one ntuple but not the other
        if run_lumi_event_number1 != run_lumi_event_number2 : continue

        # Okay! we have the same event in both ntuples here now
        mfv_event2_handle2 = Handle ("MFVEvent")
        mfv_event2_label2 = ("mfvEvent")
        event2.getByLabel (mfv_event2_label2, mfv_event2_handle2)
        mevent = mfv_event2_handle2.product()

        event2.getByLabel (vertices_label2, vertices_handle2)
        vertices_from_ntuple2 = vertices_handle2.product()

        #print "lsp2d ntuple1 %s" % mevent1.lspdist2d()
        #print "lsp2d ntuple2 %s" % mevent.lspdist2d()
        
        if  0.0150 < mevent.lspdist2d() < 2 :   # apply fiducial cuts
            if  len(vertices_from_ntuple1) < 2 :
                nevents_presel_nsv01_ntuple1 += 1
            if  len(vertices_from_ntuple2) < 2 :
                nevents_presel_nsv01_ntuple2 += 1


            
            #if math.fabs(ROOT.reco.deltaPhi(sv0_ntuple1_phi,sv1_ntuple1_phi)) > 1.57 and math.fabs(ROOT.reco.deltaPhi(sv0_ntuple2_phi,sv1_ntuple2_phi)) > 1.57:
            nsv_ntuple1 = 0
            for vtx_ntuple1 in vertices_from_ntuple1 :
                 vtx_ntuple1_phi = math.atan2(vtx_ntuple1.y - mevent.bsy_at_z(vtx_ntuple1.z), vtx_ntuple1.x - mevent.bsx_at_z(vtx_ntuple1.z))
                 if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0],vtx_ntuple1_phi)) < 1.57:
                    #if exclude_beampipe and !jmt.Geometry.inside_beampipe(is_mc, vtx_ntuple1.x, vtx_ntuple1.y): 
                    
                    dBV_ntuple1 = np.array([vtx_ntuple1.x - mevent.bsx_at_z(vtx_ntuple1.z),vtx_ntuple1.y - mevent.bsy_at_z(vtx_ntuple1.z)])
                   
                    if vtx_ntuple1.ntracks()>=5 and np.linalg.norm(dBV_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                          nsv_ntuple1 += 1
                        
            nsv_ntuple2 = 0
            for vtx_ntuple2 in vertices_from_ntuple2 :
                 vtx_ntuple2_phi = math.atan2(vtx_ntuple2.y - mevent.bsy_at_z(vtx_ntuple2.z), vtx_ntuple2.x - mevent.bsx_at_z(vtx_ntuple2.z))
                 if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0],vtx_ntuple2_phi)) < 1.57:
                    #if exclude_beampipe and !jmt.Geometry.inside_beampipe(is_mc, vtx_ntuple2.x, vtx_ntuple2.y):
                    
                    dBV_ntuple2 = np.array([vtx_ntuple2.x - mevent.bsx_at_z(vtx_ntuple2.z),vtx_ntuple2.y - mevent.bsy_at_z(vtx_ntuple2.z)])

                    if vtx_ntuple2.ntracks()>=5 and np.linalg.norm(dBV_ntuple2) > 0.0100 and vtx_ntuple2.rescale_bs2derr < 0.0025:
                          nsv_ntuple2 += 1
                            
            h_2D_nsv.Fill(nsv_ntuple1,nsv_ntuple2)
            break
                


                
            #####

            # reset event2
        events_ntuple2.toBegin()

        nevents_processed += 1


print "Total processed event (fiducial) #%s" % (nevents_processed)
print "Total nsv<2 events in ntuple1 (ficudial) #%s" % (nevents_presel_nsv01_ntuple1)
print "Total nsv<2 events in ntuple2 (ficudial) #%s" % (nevents_presel_nsv01_ntuple2)
# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
h_2D_nsv.Draw("colz")
c1.Print (outputdir+"h_2D_nsv.png")
c1.Print (outputdir+"h_2D_nsv.root")
