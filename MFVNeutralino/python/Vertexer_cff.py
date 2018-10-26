import FWCore.ParameterSet.Config as cms

from JMTucker.MFVNeutralino.GenParticles_cff import mfvGenParticles
from JMTucker.MFVNeutralino.Vertexer_cfi import mfvVertexTracks, mfvVertices
from JMTucker.MFVNeutralino.VertexAuxProducer_cfi import mfvVerticesAuxTmp, mfvVerticesAux
from JMTucker.MFVNeutralino.VertexSelector_cfi import mfvSelectedVertices
from JMTucker.MFVNeutralino.JetVertexAssociator_cfi import mfvVerticesToJets

mfvSelectedVerticesTmp = mfvSelectedVertices.clone(vertex_aux_src = 'mfvVerticesAuxTmp',
                                                   produce_refs = True,
                                                   min_ntracks = 3)

mfvVerticesAuxPresel = mfvVerticesAux
mfvVerticesAux = mfvSelectedVertices.clone(vertex_aux_src = 'mfvVerticesAuxPresel', min_ntracks = 3)

mfvVertexSequenceBare = cms.Sequence(
    mfvVertexTracks *
    mfvVertices
    )

mfvVertexSequence = cms.Sequence(
    mfvVertexSequenceBare *
    mfvGenParticles *
    mfvVerticesAuxTmp *
    mfvSelectedVerticesTmp *
    mfvVerticesToJets *
    mfvVerticesAuxPresel *
    mfvVerticesAux
    )

def modifiedVertexSequence(process, name, **kwargs):
    mfvVertexTracksNew = process.mfvVertexTracks.clone(**kwargs)
    mfvVerticesNew = process.mfvVertices.clone(seed_tracks_src = 'mfvVertexTracks%s' % name, **kwargs)
    mfvVerticesAuxTmpNew = process.mfvVerticesAuxTmp.clone(vertex_src = 'mfvVertices%s' % name)
    mfvSelectedVerticesTmpNew = process.mfvSelectedVerticesTmp.clone(vertex_src = 'mfvVertices%s' % name,
                                                                     vertex_aux_src = 'mfvVerticesAuxTmp%s' % name)
    mfvVerticesToJetsNew = process.mfvVerticesToJets.clone(vertex_src = 'mfvSelectedVerticesTmp%s' % name)
    mfvVerticesAuxPreselNew = process.mfvVerticesAuxPresel.clone(vertex_src = 'mfvVertices%s' % name,
                                                                 sv_to_jets_src = 'mfvVerticesToJets%s' % name)
    mfvVerticesAuxNew = process.mfvVerticesAux.clone(vertex_aux_src = 'mfvVerticesAuxPresel%s' % name)

    setattr(process, 'mfvVertexTracks%s'        % name, mfvVertexTracksNew)
    setattr(process, 'mfvVertices%s'            % name, mfvVerticesNew)
    setattr(process, 'mfvVerticesAuxTmp%s'      % name, mfvVerticesAuxTmpNew)
    setattr(process, 'mfvSelectedVerticesTmp%s' % name, mfvSelectedVerticesTmpNew)
    setattr(process, 'mfvVerticesToJets%s'      % name, mfvVerticesToJetsNew)
    setattr(process, 'mfvVerticesAuxPresel%s'   % name, mfvVerticesAuxPreselNew)
    setattr(process, 'mfvVerticesAux%s'         % name, mfvVerticesAuxNew)

    seq = cms.Sequence(
        mfvVertexTracksNew *
        mfvVerticesNew *
        mfvVerticesAuxTmpNew *
        mfvSelectedVerticesTmpNew *
        mfvVerticesToJetsNew *
        mfvVerticesAuxPreselNew *
        mfvVerticesAuxNew)
    return seq
