#! /usr/bin/env python

import ROOT
import sys,os
import math
import numpy as np

from DataFormats.FWLite import Events, Handle

from JMTucker.Tools.ROOTTools import *
cmssw_setup()

# FIXME you can replace this with the usual stuff for putting plots into our publicweb areas and generating the html
outputdir = "/eos/user/p/pekotamn/www/Null_SV_tbs_10mm_1.6TeV_default_ntuple"
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
h_nsv = ROOT.TH1F ("h_nsv", ";# of <5trk-SVs in LSP0 or LSP1's hemisphere", 15, 0, 15)
h_dist3d_to_lsp_SV0 = ROOT.TH1F ("h_dist3d_to_lsp_SV0", ";dist3d(<5trk-SV0, closest gen vtx) (cm)", 200, 0, 0.2)
h_significance_dist3d_to_lsp_SV0 = ROOT.TH1F ("h_significance_dist3d_to_lsp_SV0", ";N#sigma(dist3d(<5trk-SV0, closest gen vtx))", 200, 0, 100)
h_ntracks_SV0 = ROOT.TH1F ("h_ntracks_SV0", ";# of tracks/<5trk-SV0", 10, 0, 10)
h_mass_SV0 = ROOT.TH1F ("h_mass_SV0", ";SV0 tracks-plus-jets-by-ntracks mass (GeV)", 100, 0, 5000)
h_dist3d_to_lsp_SV1 = ROOT.TH1F ("h_dist3d_to_lsp_SV1", ";dist3d(<5trk-SV1, closest gen vtx) (cm)", 200, 0, 0.2)
h_significance_dist3d_to_lsp_SV1 = ROOT.TH1F ("h_significance_dist3d_to_lsp_SV1", ";N#sigma(dist3d(<5trk-SV1, closest gen vtx))", 200, 0, 100)
h_ntracks_SV1 = ROOT.TH1F ("h_ntracks_SV1", ";# of tracks/<5trk-SV1", 10, 0, 10)
h_mass_SV1 = ROOT.TH1F ("h_mass_SV1", ";SV1 tracks-plus-jets-by-ntracks mass (GeV)", 100, 0, 5000)
h_dist3d_to_lsp_SV2 = ROOT.TH1F ("h_dist3d_to_lsp_SV2", ";dist3d(<5trk-SV2, closest gen vtx) (cm)", 200, 0, 0.2)
h_significance_dist3d_to_lsp_SV2 = ROOT.TH1F ("h_significance_dist3d_to_lsp_SV2", ";N#sigma(dist3d(<5trk-SV2, closest gen vtx))", 200, 0, 100)
h_ntracks_SV2 = ROOT.TH1F ("h_ntracks_SV2", ";# of tracks/<5trk-SV2", 10, 0, 10)
h_mass_SV2 = ROOT.TH1F ("h_mass_SV2", ";SV2 tracks-plus-jets-by-ntracks mass (GeV)", 100, 0, 5000)

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
    if nevents_processed <= 5000 :
        if  0.0150 < math.sqrt((mevent.gen_lsp_decay[0])**2 + (mevent.gen_lsp_decay[1])**2) < 2 and  0.0150 < math.sqrt((mevent.gen_lsp_decay[3])**2 + (mevent.gen_lsp_decay[4])**2) < 2 and math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0], mevent.gen_lsp_phi[1])) > 2.7: # apply fiducial cuts
               nevents_fiducial_cuts += 1

               qual_nsv = 0

               ls_of_qual_nsv_lsp0 = []
               ls_of_unqual_nsv_lsp0 = []
               ls_of_qual_nsv_lsp1 = []
               ls_of_unqual_nsv_lsp1 = []

               for vtx_ntuple1 in vertices_from_ntuple1 : # loop raw verices

                   dBV_vtx_ntuple1 = np.array([vtx_ntuple1.x - mevent.bsx_at_z(vtx_ntuple1.z),vtx_ntuple1.y - mevent.bsy_at_z(vtx_ntuple1.z)])
                   
                   if vtx_ntuple1.ntracks()>=5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09 and np.linalg.norm(dBV_vtx_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                          qual_nsv += 1

                   vtx_ntuple1_phi = math.atan2(vtx_ntuple1.y - mevent.bsy_at_z(vtx_ntuple1.z), vtx_ntuple1.x - mevent.bsx_at_z(vtx_ntuple1.z))
                  
                   if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0],vtx_ntuple1_phi)) < 1.57:
                        if vtx_ntuple1.ntracks()>=5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09 and np.linalg.norm(dBV_vtx_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                           ls_of_qual_nsv_lsp0.append(vtx_ntuple1)
                        #if vtx_ntuple1.ntracks()<5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09 and np.linalg.norm(dBV_vtx_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                        if vtx_ntuple1.ntracks()<5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09:
                           ls_of_unqual_nsv_lsp0.append(vtx_ntuple1)

                   if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[1],vtx_ntuple1_phi)) < 1.57:
                        if vtx_ntuple1.ntracks()>=5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09 and np.linalg.norm(dBV_vtx_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                           ls_of_qual_nsv_lsp1.append(vtx_ntuple1)
                        #if vtx_ntuple1.ntracks()<5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09 and np.linalg.norm(dBV_vtx_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                        if vtx_ntuple1.ntracks()<5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09: 
                          ls_of_unqual_nsv_lsp1.append(vtx_ntuple1)
               

               if qual_nsv < 2:
                   nevents_nsv01_fiducial_cuts += 1

               if len(ls_of_qual_nsv_lsp0) == 0:
                   h_nsv.Fill(len(ls_of_unqual_nsv_lsp0))
                   unqual_nsv_lsp0 = 0
                   for unqual_vtx in ls_of_unqual_nsv_lsp0:
                           unqual_nsv_lsp0 += 1
                           if unqual_nsv_lsp0 == 1 :
                               h_dist3d_to_lsp_SV0.Fill(unqual_vtx.gen3ddist)
                               h_significance_dist3d_to_lsp_SV0.Fill(unqual_vtx.gen3dsig())
                               h_ntracks_SV0.Fill(unqual_vtx.ntracks())
                               h_mass_SV0.Fill(unqual_vtx.mass[ROOT.mfv.PTracksPlusJetsByNtracks])
                           
                           if unqual_nsv_lsp0 == 2 :
                               h_dist3d_to_lsp_SV1.Fill(unqual_vtx.gen3ddist)
                               h_significance_dist3d_to_lsp_SV1.Fill(unqual_vtx.gen3dsig())
                               h_ntracks_SV1.Fill(unqual_vtx.ntracks())
                               h_mass_SV1.Fill(unqual_vtx.mass[ROOT.mfv.PTracksPlusJetsByNtracks])

                           if unqual_nsv_lsp0 == 3 :
                               h_dist3d_to_lsp_SV2.Fill(unqual_vtx.gen3ddist)
                               h_significance_dist3d_to_lsp_SV2.Fill(unqual_vtx.gen3dsig())
                               h_ntracks_SV2.Fill(unqual_vtx.ntracks())
                               h_mass_SV2.Fill(unqual_vtx.mass[ROOT.mfv.PTracksPlusJetsByNtracks])

               if len(ls_of_qual_nsv_lsp1) == 0:
                   h_nsv.Fill(len(ls_of_unqual_nsv_lsp1))
                   unqual_nsv_lsp1 = 0
                   for unqual_vtx in ls_of_unqual_nsv_lsp1:
                           unqual_nsv_lsp1 += 1
                           if unqual_nsv_lsp1 == 1 :
                               h_dist3d_to_lsp_SV0.Fill(unqual_vtx.gen3ddist)
                               h_significance_dist3d_to_lsp_SV0.Fill(unqual_vtx.gen3dsig())
                               h_ntracks_SV0.Fill(unqual_vtx.ntracks())
                               h_mass_SV0.Fill(unqual_vtx.mass[ROOT.mfv.PTracksPlusJetsByNtracks])
                           
                           if unqual_nsv_lsp1 == 2 :
                               h_dist3d_to_lsp_SV1.Fill(unqual_vtx.gen3ddist)
                               h_significance_dist3d_to_lsp_SV1.Fill(unqual_vtx.gen3dsig())
                               h_ntracks_SV1.Fill(unqual_vtx.ntracks())
                               h_mass_SV1.Fill(unqual_vtx.mass[ROOT.mfv.PTracksPlusJetsByNtracks])

                           if unqual_nsv_lsp1 == 3 :
                               h_dist3d_to_lsp_SV2.Fill(unqual_vtx.gen3ddist)
                               h_significance_dist3d_to_lsp_SV2.Fill(unqual_vtx.gen3dsig())
                               h_ntracks_SV2.Fill(unqual_vtx.ntracks())
                               h_mass_SV2.Fill(unqual_vtx.mass[ROOT.mfv.PTracksPlusJetsByNtracks])
                   

    else:
        break



print "Total processed event (5000 events) #%s" % (nevents_processed)
print "Total processed event (fiducial) #%s" % (nevents_fiducial_cuts)
print "Total nsv<2 events in ntuple1 (ficudial) #%s" % (nevents_nsv01_fiducial_cuts)

# make a canvas, draw, and save it
c1 = ROOT.TCanvas()
h_nsv.Draw("colz")
c1.Print (outputdir+"h_nsv.png")
c1.Print (outputdir+"h_nsv.root")

c2 = ROOT.TCanvas()
h_dist3d_to_lsp_SV0.Draw("colz")
c2.Print (outputdir+"h_dist3d_to_lsp_SV0.png")
c2.Print (outputdir+"h_dist3d_to_lsp_SV0.root")

c3 = ROOT.TCanvas()
h_dist3d_to_lsp_SV1.Draw("colz")
c3.Print (outputdir+"h_dist3d_to_lsp_SV1.png")
c3.Print (outputdir+"h_dist3d_to_lsp_SV1.root")

c4 = ROOT.TCanvas()
h_dist3d_to_lsp_SV2.Draw("colz")
c4.Print (outputdir+"h_dist3d_to_lsp_SV2.png")
c4.Print (outputdir+"h_dist3d_to_lsp_SV2.root")

c5 = ROOT.TCanvas()
h_significance_dist3d_to_lsp_SV0.Draw("colz")
c5.Print (outputdir+"h_significance_dist3d_to_lsp_SV0.png")
c5.Print (outputdir+"h_significance_dist3d_to_lsp_SV0.root")

c6 = ROOT.TCanvas()
h_significance_dist3d_to_lsp_SV1.Draw("colz")
c6.Print (outputdir+"h_significance_dist3d_to_lsp_SV1.png")
c6.Print (outputdir+"h_significance_dist3d_to_lsp_SV1.root")

c7 = ROOT.TCanvas()
h_significance_dist3d_to_lsp_SV2.Draw("colz")
c7.Print (outputdir+"h_significance_dist3d_to_lsp_SV2.png")
c7.Print (outputdir+"h_significance_dist3d_to_lsp_SV2.root")

c8 = ROOT.TCanvas()
h_ntracks_SV0.Draw("colz")
c8.Print (outputdir+"h_ntracks_SV0.png")
c8.Print (outputdir+"h_ntracks_SV0.root")

c9 = ROOT.TCanvas()
h_ntracks_SV1.Draw("colz")
c9.Print (outputdir+"h_ntracks_SV1.png")
c9.Print (outputdir+"h_ntracks_SV1.root")

c10 = ROOT.TCanvas()
h_ntracks_SV2.Draw("colz")
c10.Print (outputdir+"h_ntracks_SV2.png")
c10.Print (outputdir+"h_ntracks_SV2.root")

c11 = ROOT.TCanvas()
h_mass_SV0.Draw("colz")
c11.Print (outputdir+"h_mass_SV0.png")
c11.Print (outputdir+"h_mass_SV0.root")

c12 = ROOT.TCanvas()
h_mass_SV1.Draw("colz")
c12.Print (outputdir+"h_mass_SV1.png")
c12.Print (outputdir+"h_mass_SV1.root")

c13 = ROOT.TCanvas()
h_mass_SV2.Draw("colz")
c13.Print (outputdir+"h_mass_SV2.png")
c13.Print (outputdir+"h_mass_SV2.root")