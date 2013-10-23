import FWCore.ParameterSet.Config as cms

mfvGenParticleFilter = cms.EDFilter('MFVGenParticleFilter',
                                    gen_src = cms.InputTag('genParticles'),
                                    print_info = cms.bool(False),
                                    cut_invalid = cms.bool(True),
                                    required_num_leptonic = cms.int32(-1),
                                    allowed_decay_types = cms.vint32(),
                                    min_lepton_pt = cms.double(0),
                                    max_lepton_eta = cms.double(1e99),
                                    min_rho0 = cms.double(-1),
                                    max_rho0 = cms.double(-1),
                                    min_rho1 = cms.double(-1),
                                    max_rho1 = cms.double(-1),
                                    min_r0 = cms.double(-1),
                                    max_r0 = cms.double(-1),
                                    min_r1 = cms.double(-1),
                                    max_r1 = cms.double(-1),
                                    min_rhobigger = cms.double(-1),
                                    max_rhobigger = cms.double(-1),
                                    min_rhosmaller = cms.double(-1),
                                    max_rhosmaller = cms.double(-1),
                                    min_rbigger = cms.double(-1),
                                    max_rbigger = cms.double(-1),
                                    min_rsmaller = cms.double(-1),
                                    max_rsmaller = cms.double(-1),
                                    )
