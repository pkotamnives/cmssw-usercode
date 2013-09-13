import FWCore.ParameterSet.Config as cms

mfvVertices = cms.EDProducer('MFVVertexer',
                             kvr_params = cms.PSet(
                                 maxDistance = cms.double(0.01),
                                 maxNbrOfIterations = cms.int32(10),
                                 doSmoothing = cms.bool(True),
                                 ),
                             avr_params = cms.PSet(
                                 finder = cms.string('avr'),
                                 primcut = cms.double(1.0),
                                 seccut = cms.double(3),
                                 smoothing = cms.bool(True)
                                 ),
                             beamspot_src = cms.InputTag('offlineBeamSpot'),
                             use_tracks = cms.bool(True),
                             track_src = cms.InputTag('generalTracks'),
                             use_pf_candidates = cms.bool(False),
                             pf_candidate_src = cms.InputTag('particleFlow'),
                             use_pf_jets = cms.bool(False),
                             pf_jet_src = cms.InputTag('ak5PFJets'),
                             use_pat_jets = cms.bool(False),
                             pat_jet_src = cms.InputTag('selectedPatJetsPF'),
                             min_seed_jet_pt = cms.double(30),
                             min_seed_track_pt = cms.double(1),
                             min_seed_track_dxy = cms.double(0.01),
                             min_seed_track_nhits = cms.int32(8),
                             max_seed_vertex_chi2 = cms.double(5),
                             use_2d_vertex_dist = cms.bool(False),
                             use_2d_track_dist = cms.bool(False),
                             merge_anyway_dist = cms.double(-1),
                             merge_anyway_sig = cms.double(3),
                             merge_shared_dist = cms.double(-1),
                             merge_shared_sig = cms.double(4),
                             max_track_vertex_dist = cms.double(-1),
                             max_track_vertex_sig = cms.double(5),
                             min_track_vertex_sig_to_remove = cms.double(1.5),
                             remove_one_track_at_a_time = cms.bool(True),
                             histos = cms.untracked.bool(False),
                             verbose = cms.untracked.bool(False),
                             )

import JMTucker.MFVNeutralino.GenParticles_cff
mfvGenVertices = JMTucker.MFVNeutralino.GenParticles_cff.mfvGenVertices.clone()

import JMTucker.MFVNeutralino.VertexSelector_cfi
mfvSelectedVertices = JMTucker.MFVNeutralino.VertexSelector_cfi.mfvSelectedVertices.clone()

mfvSelectedJets = cms.EDFilter('PATJetSelector',
                               src = cms.InputTag('selectedPatJetsPF'),
                               cut = cms.string('pt > 20')
                               )

mfvVerticesToJets = cms.EDProducer('MFVJetVertexAssociator',
                                   jet_src = cms.InputTag('mfvSelectedJets'),
                                   vertex_src = cms.InputTag('mfvSelectedVertices'),
                                   tag_info_name = cms.string('secondaryVertexMaxDR6p0'),
                                   min_vertex_track_weight = cms.double(0.5),
                                   min_tracks_shared = cms.int32(1),
                                   min_track_pt = cms.double(3),
                                   min_hits_shared = cms.int32(1),
                                   max_cos_angle_diff = cms.double(2),
                                   max_miss_dist = cms.double(1e6),
                                   max_miss_sig = cms.double(1e6),
                                   histos = cms.untracked.bool(True),
                                   verbose = cms.untracked.bool(False),
                                   )

def use_pf_candidates(vertex_producer, src):
    vertex_producer.use_tracks = False
    vertex_producer.use_pf_jets = False
    vertex_producer.use_pat_jets = False
    vertex_producer.use_pf_candidates = True
    vertex_producer.pf_candidate_src = src

def use_jets(vertex_producer, kind):
    vertex_producer.use_tracks = False
    vertex_producer.use_pf_candidates = False
    if kind == 'pat':
        vertex_producer.use_pf_jets = False
        vertex_producer.use_pat_jets = True
    elif kind == 'pf':
        vertex_producer.use_pf_jets = True
        vertex_producer.use_pat_jets = False
    else:
        raise ValueError("don't know anything about kind = %r" % kind)

mfvVertexSequence = cms.Sequence(mfvVertices * mfvGenVertices * mfvSelectedVertices * mfvSelectedJets * mfvVerticesToJets)

mfvVerticesFromCands = mfvVertices.clone()
use_pf_candidates(mfvVerticesFromCands, 'particleFlow')

mfvVerticesFromNoPUCands = mfvVertices.clone()
use_pf_candidates(mfvVerticesFromNoPUCands, 'pfNoPileUpPF')

mfvVerticesFromNoPUZCands = mfvVertices.clone()
use_pf_candidates(mfvVerticesFromNoPUZCands, 'pfNoPileUpPFClosestZVertex')

mfvVerticesFromJets = mfvVertices.clone()
use_jets(mfvVerticesFromJets, 'pat')

mfvVerticesFromPFJets = mfvVertices.clone()
use_jets(mfvVerticesFromPFJets, 'pf')

mfvPFCandVertexSequence = cms.Sequence(mfvVerticesFromCands + mfvVerticesFromNoPUCands + mfvVerticesFromNoPUZCands)
mfvJetVertexSequence = cms.Sequence(mfvVerticesFromJets + mfvVerticesFromPFJets)
mfvExtraVertexSequence = cms.Sequence(mfvPFCandVertexSequence + mfvJetVertexSequence)

mfvNonPATJetVertexSequence = cms.Sequence(mfvVerticesFromPFJets)
mfvNonPATExtraVertexSequence = cms.Sequence(mfvPFCandVertexSequence + mfvNonPATJetVertexSequence)
