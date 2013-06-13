import os, sys
from JMTucker.Tools.BasicAnalyzer_cfg import cms, process, add_analyzer

process.maxEvents.input = 1000
process.source.fileNames = ['file:/uscms/home/tucker/nobackup/fromt3/mfv_neutralino_tau1000um_M0400_jtuple_v6_547d3313903142038335071634b26604_pat_1_1_Dpa.root']
process.TFileService.fileName = 'track_play.root'

#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'START53_V21::All'
#process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

ana = process.mfvTrackingPlay = cms.EDAnalyzer('MFVTrackPlay',
                                               tracks_src = cms.InputTag('generalTracks'),
                                               dxy_cut = cms.double(0),
                                               quality_cut = cms.string(''),
                                               )
process.p = cms.Path(process.mfvTrackingPlay)

cuts = [
    ('Qpt10eta2p4',        'pt > 10 && abs(eta) < 2.4'),
    ('Qpt50eta2p4',        'pt > 50 && abs(eta) < 2.4'),
    ('Qpt10eta2p4dxy0p1',  'pt > 10 && abs(eta) < 2.4'),
    ('Qpt10eta2p4dxy0p5',  'pt > 10 && abs(eta) < 2.4'),
    ('Qpt10eta2p4highPur', 'pt > 10 && abs(eta) < 2.4'),
    ]

for cut_name, cut in cuts:
    filt = cms.EDFilter('TrackSelector',
                        src = ana.tracks_src,
                        filter = cms.bool(False),
                        cut = cms.string(cut),
                        )
    setattr(process, cut_name, filt)
    cut_ana = ana.clone(tracks_src = cut_name)
    if 'dxy0p1' in cut_name:
        cut_ana.dxy_cut = 0.1
    if 'dxy0p5' in cut_name:
        cut_ana.dxy_cut = 0.5
    if 'highPur' in cut_name:
        cut_ana.quality_cut = 'highPurity'
    setattr(process, 'mfvTrackingPlay' + cut_name, cut_ana)
    process.p *= filt * cut_ana

#process.add_(cms.Service('SimpleMemoryCheck'))

if __name__ == '__main__' and hasattr(sys, 'argv') and 'submit' in sys.argv:
    if 'debug' in sys.argv:
        raise RuntimeError('refusing to submit jobs in debug (verbose print out) mode')

    crab_cfg = '''
[CRAB]
jobtype = cmssw
scheduler = %(scheduler)s 

[CMSSW]
%(dbs_url)s
datasetpath = %(dataset)s
pset = track_play_crab.py
total_number_of_events = 100000
events_per_job = 20000

[USER]
ui_working_dir = crab/TrackPlay/crab_mfv_trackplay_%(name)s
jmt_skip_input_files = src/EGamma/EGammaAnalysisTools/data/*
return_data = 1
'''

    os.system('mkdir -p crab/TrackPlay')
    
    testing = 'testing' in sys.argv
    from JMTucker.Tools.Samples import mfv_neutralino_tau0000um_M0400, mfv_neutralino_tau1000um_M0400, mfv_neutralino_tau9900um_M0400, ttbarincl, TupleOnlyMCSample
    samples = [mfv_neutralino_tau0000um_M0400, mfv_neutralino_tau1000um_M0400, mfv_neutralino_tau9900um_M0400, ttbarincl]

    mfv_neutralino_tau9900um_M0400.dataset = '/crabfake_mfv_neutralino_tau9900um_M0400_jtuple_v6_547d3313903142038335071634b26604/tucker-crabfake_mfv_neutralino_tau9900um_M0400_jtuple_v6_547d3313903142038335071634b26604-5bdce5833f35b995ab0c308220e77250/USER'
    mfv_neutralino_tau1000um_M0400.dataset = '/crabfake_mfv_neutralino_tau1000um_M0400_jtuple_v6_547d3313903142038335071634b26604/tucker-crabfake_mfv_neutralino_tau1000um_M0400_jtuple_v6_547d3313903142038335071634b26604-5bdce5833f35b995ab0c308220e77250/USER'
    mfv_neutralino_tau0000um_M0400.dataset = '/crabfake_mfv_neutralino_tau0000um_M0400_jtuple_v6_547d3313903142038335071634b26604/tucker-crabfake_mfv_neutralino_tau0000um_M0400_jtuple_v6_547d3313903142038335071634b26604-5bdce5833f35b995ab0c308220e77250/USER'
    
    for sample in samples:
        sample.scheduler_name = 'condor'
        open('crab.cfg', 'wt').write(crab_cfg % sample)
        new_py = open('track_play.py').read()
        open('track_play_crab.py', 'wt').write(new_py)
        if not testing:
            os.system('crab -create -submit')
            os.system('rm crab.cfg track_play_crab.py track_play_crab.pyc')
