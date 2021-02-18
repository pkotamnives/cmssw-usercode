#! /usr/bin/env python

import ROOT
import sys,os
import math
import numpy as np

from DataFormats.FWLite import Events, Handle

from JMTucker.Tools.ROOTTools import *
cmssw_setup()

# FIXME you can replace this with the usual stuff for putting plots into our publicweb areas and generating the html
outputdir = "/eos/user/p/pekotamn/www/Null_SV_tbs_1mm_1.6TeV_name_100n_ntuple"
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
h_unqual_nsv = ROOT.TH1F ("h_unqual_nsv", ";Categories of unqualified SVs/LSP", 7, 0, 7)
h_qual_nsv = ROOT.TH1F ("h_qual_nsv", ";# of >=5trk-SVs/LSP", 15, 0, 15)
h_nsv = ROOT.TH1F ("h_nsv", ";# of unqualified SVs/LSP", 15, 0, 15)
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

h_non_qual_nsv_lsp_r =  ROOT.TH1F ("h_non_qual_nsv_lsp_r", ";LSP r w.r.t bsp (cm)", 100, 0, 0.5)
h_non_qual_nsv_lsp_z =  ROOT.TH1F ("h_non_qual_nsv_lsp_z", ";LSP z (cm)", 100, -25, 25)
h_non_qual_nsv_lsp_seed_tracks =  ROOT.TH1F ("h_non_qual_nsv_lsp_seed_tracks", ";# of seed tracks/LSP", 50, 0, 50)
h_qual_nsv_lsp_r =  ROOT.TH1F ("h_qual_nsv_lsp_r", ";LSP r w.r.t bsp  (cm)", 100, 0, 0.5)
h_qual_nsv_lsp_z =  ROOT.TH1F ("h_qual_nsv_lsp_z", ";LSP z (cm)", 100, -25, 25)
h_qual_nsv_lsp_seed_tracks =  ROOT.TH1F ("h_qual_nsv_lsp_seed_tracks", ";# of seed tracks/LSP", 50, 0, 50)

h_qual_nsv_distr_sv_lsp =  ROOT.TH1F ("h_qual_nsv_distr_sv_lsp", ";>=5trk-SV's r - closest gen LSP's r w.r.t bsp (cm)", 200, -0.08, 0.08)
h_qual_nsv_dphi_sv_lsp =  ROOT.TH1F ("h_qual_nsv_dphi_sv_lsp", ";deltaPhi(>=5trk-SV, closest gen LSP)", 314, 0, 3.14)
h_qual_nsv_dist3d_sv_lsp = ROOT.TH1F ("h_qual_nsv_dist3d_sv_lsp", ";dist3d(>=5trk-SV, closest gen LSP) (cm)", 200, 0, 0.2)
h_qual_nsv_significance_dist3d_sv_lsp = ROOT.TH1F ("h_qual_nsv_significance_dist3d_sv_lsp", ";N#sigma(dist3d(>=5trk-SV, closest gen LSP))", 200, 0, 100)

h_unqual_type4_dBV =  ROOT.TH1F("h_unqual_type4_dBV", ";dist2d(beamspot, unqual >=5trk-SV in a hemisphere) (cm);arb. units", 100, 0, 0.02); 
h_unqual_type4_r =  ROOT.TH1F("h_unqual_type4_r", "; unqual >=5trk-SV's r (cm);arb. units", 100, 0, 2);
h_unqual_type4_LSP_r =  ROOT.TH1F("h_unqual_type4_LSP_r", "; LSP's r  w.r.t bsp (cm);arb. units", 100, 0, 2);
h_unqual_type4_LSP_r_ctr =  ROOT.TH1F("h_unqual_type4_LSP_r_ctr", "; LSP's r  w.r.t origin (cm);arb. units", 100, 0, 2);
h_unqual_type4_beamspot_r =  ROOT.TH1F("h_unqual_type4_beamspot_r", "; beamspot's r (cm);arb. units", 100, 0, 2);
h_unqual_type4_bs2derr =  ROOT.TH1F("h_unqual_type4_bs2derr", "; unqual >=5trk-SV's bs2derr (cm);arb. units", 100, 0, 0.01);
h_unqual_type4_distr_sv_lsp =  ROOT.TH1F("h_unqual_type4_distr_sv_lsp", "; unqual >=5trk-SV's r - closest gen LSP's r  w.r.t bsp  (cm);arb. units", 200, -0.08, 0.08);
h_unqual_type4_dist3d_to_lsp  = ROOT.TH1F ("h_unqual_type4_dist3d_to_lsp", ";dist3d(unqual >=5trk-SV, closest gen vtx) (cm)", 200, 0, 0.2)
h_unqual_type4_dist2d_to_lsp  = ROOT.TH1F ("h_unqual_type4_dist2d_to_lsp", ";dist2d(unqual >=5trk-SV, closest gen vtx) (cm)", 200, 0, 0.2)
h_non_qual_nsv_distr_sv_lsp =  ROOT.TH1F ("h_non_qual_nsv_distr_sv_lsp", ";unqualified-SV's r - closest gen LSP's r  w.r.t bsp (cm)", 200, -0.08, 0.08)


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
    if nevents_processed <= 1000 :
         vtx_bsp = vertices_from_ntuple1[0]    # change vertex when a beamspot is not constant
         print "the beam spot's r is #%s" % ( np.linalg.norm(np.array([mevent.bsx_at_z(vtx_bsp.z),mevent.bsy_at_z(vtx_bsp.z)]))   )
         dBLSP0_vtx_ntuple1 = np.array([mevent.gen_lsp_decay[0] - mevent.bsx_at_z(vtx_bsp.z),mevent.gen_lsp_decay[1] - mevent.bsy_at_z(vtx_bsp.z)])     
         dBLSP1_vtx_ntuple1 = np.array([mevent.gen_lsp_decay[3] - mevent.bsx_at_z(vtx_bsp.z),mevent.gen_lsp_decay[4] - mevent.bsy_at_z(vtx_bsp.z)]) 
               
         if 0.0100 < np.linalg.norm(dBLSP0_vtx_ntuple1) and math.sqrt((mevent.gen_lsp_decay[0])**2 + (mevent.gen_lsp_decay[1])**2) < 2.09 and  0.0100 < np.linalg.norm(dBLSP1_vtx_ntuple1) and math.sqrt((mevent.gen_lsp_decay[3])**2 + (mevent.gen_lsp_decay[4])**2) < 2.09 and math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0], mevent.gen_lsp_phi[1])) > 2.7: # apply fiducial cuts
        #if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0], mevent.gen_lsp_phi[1])) > 2.7: # no apply fiducial cuts
            
               nevents_fiducial_cuts += 1
               n_vertex_seed_tracks = mevent.n_vertex_seed_tracks()
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
                        #if vtx_ntuple1.ntracks()<5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09:
                        else:
                           ls_of_unqual_nsv_lsp0.append(vtx_ntuple1)

                   if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[1],vtx_ntuple1_phi)) < 1.57:
                        if vtx_ntuple1.ntracks()>=5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09 and np.linalg.norm(dBV_vtx_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                           ls_of_qual_nsv_lsp1.append(vtx_ntuple1)
                        #if vtx_ntuple1.ntracks()<5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09 and np.linalg.norm(dBV_vtx_ntuple1) > 0.0100 and vtx_ntuple1.rescale_bs2derr < 0.0025:
                        #if vtx_ntuple1.ntracks()<5 and math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2) < 2.09: 
                        else:
                          ls_of_unqual_nsv_lsp1.append(vtx_ntuple1)
               
                   
               if qual_nsv < 2:
                   nevents_nsv01_fiducial_cuts += 1

               h_qual_nsv.Fill(len(ls_of_qual_nsv_lsp0))
               h_qual_nsv.Fill(len(ls_of_qual_nsv_lsp1))
               h_qual_nsv_event.Fill(len(ls_of_qual_nsv_lsp0)+len(ls_of_qual_nsv_lsp1))

               if len(ls_of_qual_nsv_lsp0) == 0:
                   h_nsv.Fill(len(ls_of_unqual_nsv_lsp0))
                   unqual_nsv_lsp0 = 0
                   for unqual_vtx in ls_of_unqual_nsv_lsp0:
                           unqual_nsv_lsp0 += 1
                           dBV_vtx_ntuple1 = np.array([unqual_vtx.x - mevent.bsx_at_z(unqual_vtx.z),unqual_vtx.y - mevent.bsy_at_z(unqual_vtx.z)])
                           h_non_qual_nsv_distr_sv_lsp.Fill(np.linalg.norm(dBV_vtx_ntuple1)- np.linalg.norm(dBLSP0_vtx_ntuple1))
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

                   h_non_qual_nsv_lsp_r.Fill(np.linalg.norm(dBLSP0_vtx_ntuple1))
                   h_non_qual_nsv_lsp_z.Fill(mevent.gen_lsp_decay[2])

                   
                   count_seed_tracks = 0
                   for i in range(n_vertex_seed_tracks):
                       if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0],mevent.vertex_seed_track_phi[i])) < 1.57 :
                           count_seed_tracks += 1
                   h_non_qual_nsv_lsp_seed_tracks.Fill(count_seed_tracks)

                   if len(ls_of_unqual_nsv_lsp0) == 0:
                        h_unqual_nsv.Fill(1)
                   if len(ls_of_unqual_nsv_lsp0) == 1:
                        if ls_of_unqual_nsv_lsp0[0].ntracks()==3:
                            h_unqual_nsv.Fill(2)
                        if ls_of_unqual_nsv_lsp0[0].ntracks()==4:
                            h_unqual_nsv.Fill(3)
                        if ls_of_unqual_nsv_lsp0[0].ntracks()>=5:
                            h_unqual_nsv.Fill(4)
                            vtx_ntuple1 = ls_of_unqual_nsv_lsp0[0]
                            dBV_vtx_ntuple1 = np.array([vtx_ntuple1.x - mevent.bsx_at_z(vtx_ntuple1.z),vtx_ntuple1.y - mevent.bsy_at_z(vtx_ntuple1.z)])
                            h_unqual_type4_dBV.Fill(np.linalg.norm(dBV_vtx_ntuple1))
                            h_unqual_type4_r.Fill(math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2))
                            dBLSP_vtx_ntuple1 = np.array([mevent.gen_lsp_decay[0] - mevent.bsx_at_z(vtx_ntuple1.z),mevent.gen_lsp_decay[1] - mevent.bsy_at_z(vtx_ntuple1.z)])     # change vertex when a beamspot is not constant 
                            h_unqual_type4_LSP_r.Fill(np.linalg.norm(dBLSP_vtx_ntuple1))
                            h_unqual_type4_LSP_r_ctr.Fill(math.sqrt(mevent.gen_lsp_decay[0]**2 + mevent.gen_lsp_decay[1]**2))
                            h_unqual_type4_beamspot_r.Fill(math.sqrt(mevent.bsx_at_z(vtx_ntuple1.z)**2 + mevent.bsy_at_z(vtx_ntuple1.z)))
                            h_unqual_type4_bs2derr.Fill(vtx_ntuple1.rescale_bs2derr)
                            h_unqual_type4_distr_sv_lsp.Fill(np.linalg.norm(dBV_vtx_ntuple1)-np.linalg.norm(dBLSP_vtx_ntuple1))
                            h_unqual_type4_dist3d_to_lsp.Fill(vtx_ntuple1.gen3ddist)
                            h_unqual_type4_dist2d_to_lsp.Fill(vtx_ntuple1.gen2ddist)

                   if len(ls_of_unqual_nsv_lsp0) == 2:
                        if ls_of_unqual_nsv_lsp0[0].ntracks()<5 and ls_of_unqual_nsv_lsp0[1].ntracks()<5 :
                            h_unqual_nsv.Fill(5)
                        if ls_of_unqual_nsv_lsp0[0].ntracks()>=5 or ls_of_unqual_nsv_lsp0[1].ntracks()>=5 :
                            h_unqual_nsv.Fill(6)

               else:

                   dBLSP_vtx_ntuple1 = np.array([mevent.gen_lsp_decay[0] - mevent.bsx_at_z(vtx_ntuple1.z),mevent.gen_lsp_decay[1] - mevent.bsy_at_z(vtx_ntuple1.z)])          # change vertex when a beamspot is not constant 
                   h_qual_nsv_lsp_r.Fill(np.linalg.norm(dBLSP_vtx_ntuple1))
                   h_qual_nsv_lsp_z.Fill(mevent.gen_lsp_decay[2])
                   
                   for i in range(len(ls_of_qual_nsv_lsp0)):
                       sv_phi = math.atan2(ls_of_qual_nsv_lsp0[i].y - mevent.bsy_at_z(ls_of_qual_nsv_lsp0[i].z), ls_of_qual_nsv_lsp0[i].x - mevent.bsx_at_z(ls_of_qual_nsv_lsp0[i].z))
                       dBV_vtx_ntuple1 = np.array([ls_of_qual_nsv_lsp0[i].x - mevent.bsx_at_z(ls_of_qual_nsv_lsp0[i].z),ls_of_qual_nsv_lsp0[i].y - mevent.bsy_at_z(ls_of_qual_nsv_lsp0[i].z)])
                       h_qual_nsv_distr_sv_lsp.Fill(np.linalg.norm(dBV_vtx_ntuple1)-np.linalg.norm(dBLSP0_vtx_ntuple1))
                       h_qual_nsv_dphi_sv_lsp.Fill(math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0],sv_phi)))
                       h_qual_nsv_dist3d_sv_lsp.Fill(ls_of_qual_nsv_lsp0[i].gen3ddist)
                       h_qual_nsv_significance_dist3d_sv_lsp.Fill(ls_of_qual_nsv_lsp0[i].gen3dsig())
                   
                   count_seed_tracks = 0
                   for i in range(n_vertex_seed_tracks):
                       if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[0],mevent.vertex_seed_track_phi[i])) < 1.57 :
                           count_seed_tracks += 1
                   h_qual_nsv_lsp_seed_tracks.Fill(count_seed_tracks)

               if len(ls_of_qual_nsv_lsp1) == 0:
                   h_nsv.Fill(len(ls_of_unqual_nsv_lsp1))
                   unqual_nsv_lsp1 = 0
                   for unqual_vtx in ls_of_unqual_nsv_lsp1:
                           unqual_nsv_lsp1 += 1
                           dBV_vtx_ntuple1 = np.array([unqual_vtx.x - mevent.bsx_at_z(unqual_vtx.z),unqual_vtx.y - mevent.bsy_at_z(unqual_vtx.z)])
                           h_non_qual_nsv_distr_sv_lsp.Fill(np.linalg.norm(dBV_vtx_ntuple1)- np.linalg.norm(dBLSP1_vtx_ntuple1))
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

                   h_non_qual_nsv_lsp_r.Fill(np.linalg.norm(dBLSP1_vtx_ntuple1))
                   h_non_qual_nsv_lsp_z.Fill(mevent.gen_lsp_decay[5])
                   
                   count_seed_tracks = 0
                   for i in range(n_vertex_seed_tracks):
                       if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[1],mevent.vertex_seed_track_phi[i])) < 1.57 :
                           count_seed_tracks += 1
                   h_non_qual_nsv_lsp_seed_tracks.Fill(count_seed_tracks)

                   if len(ls_of_unqual_nsv_lsp1) == 0:
                        h_unqual_nsv.Fill(1)
                   if len(ls_of_unqual_nsv_lsp1) == 1:
                        if ls_of_unqual_nsv_lsp1[0].ntracks()==3:
                            h_unqual_nsv.Fill(2)
                        if ls_of_unqual_nsv_lsp1[0].ntracks()==4:
                            h_unqual_nsv.Fill(3)
                        if ls_of_unqual_nsv_lsp1[0].ntracks()>=5:
                            h_unqual_nsv.Fill(4)
                            vtx_ntuple1 = ls_of_unqual_nsv_lsp1[0]
                            dBV_vtx_ntuple1 = np.array([vtx_ntuple1.x - mevent.bsx_at_z(vtx_ntuple1.z),vtx_ntuple1.y - mevent.bsy_at_z(vtx_ntuple1.z)])
                            h_unqual_type4_dBV.Fill(np.linalg.norm(dBV_vtx_ntuple1))
                            h_unqual_type4_r.Fill(math.sqrt(vtx_ntuple1.x**2 + vtx_ntuple1.y**2))
                            dBLSP_vtx_ntuple1 = np.array([mevent.gen_lsp_decay[3] - mevent.bsx_at_z(vtx_ntuple1.z),mevent.gen_lsp_decay[4] - mevent.bsy_at_z(vtx_ntuple1.z)])     # change vertex when a beamspot is not constant 
                            h_unqual_type4_LSP_r.Fill(np.linalg.norm(dBLSP_vtx_ntuple1))
                            h_unqual_type4_LSP_r_ctr.Fill(math.sqrt(mevent.gen_lsp_decay[3]**2 + mevent.gen_lsp_decay[4]**2))
                            h_unqual_type4_beamspot_r.Fill(math.sqrt(mevent.bsx_at_z(vtx_ntuple1.z)**2 + mevent.bsy_at_z(vtx_ntuple1.z)))
                            h_unqual_type4_bs2derr.Fill(vtx_ntuple1.rescale_bs2derr)
                            h_unqual_type4_distr_sv_lsp.Fill(np.linalg.norm(dBV_vtx_ntuple1)-np.linalg.norm(dBLSP_vtx_ntuple1))
                            h_unqual_type4_dist3d_to_lsp.Fill(vtx_ntuple1.gen3ddist)
                            h_unqual_type4_dist2d_to_lsp.Fill(vtx_ntuple1.gen2ddist)

                   if len(ls_of_unqual_nsv_lsp1) == 2:
                        if ls_of_unqual_nsv_lsp1[0].ntracks()<5 and ls_of_unqual_nsv_lsp1[1].ntracks()<5 :
                            h_unqual_nsv.Fill(5)
                        if ls_of_unqual_nsv_lsp1[0].ntracks()>=5 or ls_of_unqual_nsv_lsp1[1].ntracks()>=5 :
                            h_unqual_nsv.Fill(6)

               else:

                   dBLSP_vtx_ntuple1 = np.array([mevent.gen_lsp_decay[3] - mevent.bsx_at_z(vtx_ntuple1.z),mevent.gen_lsp_decay[4] - mevent.bsy_at_z(vtx_ntuple1.z)])          # change vertex when a beamspot is not constant 
                   h_qual_nsv_lsp_r.Fill(np.linalg.norm(dBLSP_vtx_ntuple1))
                   h_qual_nsv_lsp_z.Fill(mevent.gen_lsp_decay[5])

                   for i in range(len(ls_of_qual_nsv_lsp1)):
                       sv_phi = math.atan2(ls_of_qual_nsv_lsp1[i].y - mevent.bsy_at_z(ls_of_qual_nsv_lsp1[i].z), ls_of_qual_nsv_lsp1[i].x - mevent.bsx_at_z(ls_of_qual_nsv_lsp1[i].z))
                       dBV_vtx_ntuple1 = np.array([ls_of_qual_nsv_lsp1[i].x - mevent.bsx_at_z(ls_of_qual_nsv_lsp1[i].z),ls_of_qual_nsv_lsp1[i].y - mevent.bsy_at_z(ls_of_qual_nsv_lsp1[i].z)])
                       h_qual_nsv_distr_sv_lsp.Fill(np.linalg.norm(dBV_vtx_ntuple1)-np.linalg.norm(dBLSP1_vtx_ntuple1))
                       h_qual_nsv_dphi_sv_lsp.Fill(math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[1],sv_phi)))
                       h_qual_nsv_dist3d_sv_lsp.Fill(ls_of_qual_nsv_lsp1[i].gen3ddist)
                       h_qual_nsv_significance_dist3d_sv_lsp.Fill(ls_of_qual_nsv_lsp1[i].gen3dsig())
                   
                   count_seed_tracks = 0
                   for i in range(n_vertex_seed_tracks):
                       if math.fabs(ROOT.reco.deltaPhi(mevent.gen_lsp_phi[1],mevent.vertex_seed_track_phi[i])) < 1.57 :
                           count_seed_tracks += 1
                   h_qual_nsv_lsp_seed_tracks.Fill(count_seed_tracks)
                   

    else:
        break



print "Total processed event #%s" % (nevents_processed-1)
print "Total processed event (fiducial) #%s" % (nevents_fiducial_cuts)
print "Total nsv<2 events (ficudial) #%s" % (nevents_nsv01_fiducial_cuts)

# make a canvas, draw, and save it

c00 = ROOT.TCanvas()
h_unqual_nsv.Draw("colz")
c00.Print (outputdir+"h_unqual_nsv.png")
c00.Print (outputdir+"h_unqual_nsv.root")

c01 = ROOT.TCanvas()                                                      
h_qual_nsv_event.Draw("colz")
c01.Print (outputdir+"h_qual_nsv_event.png")
c01.Print (outputdir+"h_qual_nsv_event.root")

c0 = ROOT.TCanvas()
h_qual_nsv.Draw("colz")
c0.Print (outputdir+"h_qual_nsv.png")
c0.Print (outputdir+"h_qual_nsv.root")

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

c14 = ROOT.TCanvas()
h_non_qual_nsv_lsp_r.Draw("colz")
c14.Print (outputdir+"h_non_qual_nsv_lsp_r.png")
c14.Print (outputdir+"h_non_qual_nsv_lsp_r.root")

c15 = ROOT.TCanvas()
h_non_qual_nsv_lsp_z.Draw("colz")
c15.Print (outputdir+"h_non_qual_nsv_lsp_z.png")
c15.Print (outputdir+"h_non_qual_nsv_lsp_z.root")

c16 = ROOT.TCanvas()
h_non_qual_nsv_lsp_seed_tracks.Draw("colz")
c16.Print (outputdir+"h_non_qual_nsv_lsp_seed_tracks.png")
c16.Print (outputdir+"h_non_qual_nsv_lsp_seed_tracks.root")

c17 = ROOT.TCanvas()
h_qual_nsv_lsp_r.Draw("colz")
c17.Print (outputdir+"h_qual_nsv_lsp_r.png")
c17.Print (outputdir+"h_qual_nsv_lsp_r.root")

c18 = ROOT.TCanvas()
h_qual_nsv_lsp_z.Draw("colz")
c18.Print (outputdir+"h_qual_nsv_lsp_z.png")
c18.Print (outputdir+"h_qual_nsv_lsp_z.root")

c19 = ROOT.TCanvas()
h_qual_nsv_lsp_seed_tracks.Draw("colz")
c19.Print (outputdir+"h_qual_nsv_lsp_seed_tracks.png")
c19.Print (outputdir+"h_qual_nsv_lsp_seed_tracks.root")

c20 = ROOT.TCanvas()
h_qual_nsv_distr_sv_lsp.Draw("colz")
c20.Print (outputdir+"h_qual_nsv_distr_sv_lsp.png")
c20.Print (outputdir+"h_qual_nsv_distr_sv_lsp.root")

c21 = ROOT.TCanvas()
h_qual_nsv_dphi_sv_lsp.Draw("colz")
c21.Print (outputdir+"h_qual_nsv_dphi_sv_lsp.png")
c21.Print (outputdir+"h_qual_nsv_dphi_sv_lsp.root")

c22 = ROOT.TCanvas()
h_qual_nsv_dist3d_sv_lsp.Draw("colz")
c22.Print (outputdir+"h_qual_nsv_dist3d_sv_lsp.png")
c22.Print (outputdir+"h_qual_nsv_dist3d_sv_lsp.root")

c23 = ROOT.TCanvas()
h_qual_nsv_significance_dist3d_sv_lsp.Draw("colz")
c23.Print (outputdir+"h_qual_nsv_significance_dist3d_sv_lsp.png")
c23.Print (outputdir+"h_qual_nsv_significance_dist3d_sv_lsp.root")

c24 = ROOT.TCanvas()
h_unqual_type4_dBV.Draw("colz")
c24.Print (outputdir+"h_unqual_type4_dBV.png")
c24.Print (outputdir+"h_unqual_type4_dBV.root")

c25 = ROOT.TCanvas()
h_unqual_type4_r.Draw("colz")
c25.Print (outputdir+"h_unqual_type4_r.png")
c25.Print (outputdir+"h_unqual_type4_r.root")

c26 = ROOT.TCanvas()
h_unqual_type4_bs2derr.Draw("colz")
c26.Print (outputdir+"h_unqual_type4_bs2derr.png")
c26.Print (outputdir+"h_unqual_type4_bs2derr.root")

c27 = ROOT.TCanvas()
h_unqual_type4_distr_sv_lsp.Draw("colz")
c27.Print (outputdir+"h_unqual_type4_distr_sv_lsp.png")
c27.Print (outputdir+"h_unqual_type4_distr_sv_lsp.root")

c28 = ROOT.TCanvas()
h_non_qual_nsv_distr_sv_lsp.Draw("colz")
c28.Print (outputdir+"h_non_qual_nsv_distr_sv_lsp.png")
c28.Print (outputdir+"h_non_qual_nsv_distr_sv_lsp.root")

c29 = ROOT.TCanvas()
h_unqual_type4_beamspot_r.Draw("colz")
c29.Print (outputdir+"h_unqual_type4_beamspot_r.png")
c29.Print (outputdir+"h_unqual_type4_beamspot_r.root")

c30 = ROOT.TCanvas()
h_unqual_type4_LSP_r.Draw("colz")
c30.Print (outputdir+"h_unqual_type4_LSP_r.png")
c30.Print (outputdir+"h_unqual_type4_LSP_r.root")

c31 = ROOT.TCanvas()
h_unqual_type4_dist3d_to_lsp.Draw("colz")
c31.Print (outputdir+"h_unqual_type4_dist3d_to_lsp.png")
c31.Print (outputdir+"h_unqual_type4_dist3d_to_lsp.root")

c32 = ROOT.TCanvas()
h_unqual_type4_dist2d_to_lsp.Draw("colz")
c32.Print (outputdir+"h_unqual_type4_dist2d_to_lsp.png")
c32.Print (outputdir+"h_unqual_type4_dist2d_to_lsp.root")

c33 = ROOT.TCanvas()
h_unqual_type4_LSP_r_ctr.Draw("colz")
c33.Print (outputdir+"h_unqual_type4_LSP_r_ctr.png")
c33.Print (outputdir+"h_unqual_type4_LSP_r_ctr.root")



