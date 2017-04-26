import FWCore.ParameterSet.Config as cms
from Year import year

pileup_weights = {
    2015: [ 0.543103, 0.683627, 1.06283, 1.36101, 1.6323, 1.9294, 1.45325, 1.28716, 1.38616, 1.36714, 1.26216, 1.16094, 1.05804, 0.904997, 0.701966, 0.495774, 0.334415, 0.244172, 0.213622, 0.185556, 0.120686, 0.054814, 0.0187897, 0.00535259, 0.00152566, 0.000508218, 0.000208025, 9.61443e-05, 4.83211e-05, 2.51288e-05, 1.35832e-05, 7.12796e-06, 3.80237e-06, 1.96791e-06, 1.01934e-06, 5.04546e-07, 2.465e-07, 1.11323e-07, 4.79865e-08, 2.1649e-08, 7.19711e-09, 2.9701e-09, 8.03544e-10, 3.4155e-10, 9.65737e-11, 2.12712e-11, 7.34873e-12, 1.09032e-12, ],
    2016: [ 0.344207, 0.918094, 1.19356, 0.948006, 1.12332, 1.16726, 0.785461, 0.493116, 0.742028, 0.879304, 0.963246, 1.07105, 1.12495, 1.17388, 1.20124, 1.20699, 1.19942, 1.18401, 1.1452, 1.09663, 1.06599, 1.05092, 1.0514, 1.05105, 1.0493, 1.05901, 1.07264, 1.08203, 1.09413, 1.10867, 1.09382, 1.08406, 1.04061, 0.985874, 0.909698, 0.820672, 0.716504, 0.610538, 0.503291, 0.404671, 0.310418, 0.228336, 0.164019, 0.113181, 0.0773672, 0.050999, 0.0318103, 0.0201421, 0.012261, 0.00741128, 0.00438758, 0.002621, 0.00156723, 0.000973837, 0.000726163, 0.000674079, 0.000723403, 0.000945946, 0.00134649, 0.00183786, 0.00317501, 0.00387694, 0.00471, 0.00507489, 0.00587621, 0.00511519, 0.0055894, 0.00421601, 0.00409846, 0.00367874, 0.00299944, 0.00276595, 0.00233498, 0.00191472, 0.00173505, ],
    '15MCto16MC': [ 5.76204, 15.3144, 6.27469, 4.24925, 5.0273, 7.53141, 18.794, 29.2987, 44.1202, 42.5113, 28.3124, 17.1247, 10.0655, 5.91514, 3.30473, 1.725, 0.868268, 0.432714, 0.220083, 0.124471, 0.0845176, 0.066368, 0.0560245, 0.0487893, 0.0420591, 0.0367177, 0.0316873, 0.0272032, 0.0228623, 0.0189862, 0.0147349, 0.0116224, 0.00848947, 0.00611264, 0.00413508, 0.00276414, 0.00173674, 0.00112434, 0.000706043, 0.000405201, 0.000290442, 0.000159825, 0.000125765, 6.16102e-05, 4.3736e-05, 3.79775e-05, 1.97114e-05, 2.35089e-05, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08, 4.46766e-08 ],
    '2016BCD': [ 0.0084964, 0.495008, 0.911958, 1.06983, 1.38199, 1.45624, 1.12292, 0.91547, 1.73947, 2.203, 2.09853, 1.93158, 1.85658, 1.82924, 1.74761, 1.64912, 1.58912, 1.55964, 1.51132, 1.44009, 1.3649, 1.28189, 1.20741, 1.13193, 1.05054, 0.96824, 0.873805, 0.763084, 0.64812, 0.535762, 0.419438, 0.320844, 0.230899, 0.159118, 0.103615, 0.0640486, 0.0372565, 0.0206316, 0.0108421, 0.00549721, 0.00266001, 0.00125454, 0.000603023, 0.000304827, 0.000180231, 0.000130545, 0.000114322, 0.000121891, 0.000139973, 0.000171536, 0.000216452, 0.00028497, 0.000380288, 0.000512733, 0.000750859, 0.0011413, 0.00163102, 0.0024407, 0.00367205, 0.00511987, 0.00891572, 0.0109194, 0.0132809, 0.0143162, 0.0165798, 0.0144337, 0.0157723, 0.011897, 0.0115654, 0.0103811, 0.00846416, 0.00780529, 0.00658914, 0.00540318, 0.00489617, ],
    '2016EF': [ 0.0103252, 0.739282, 1.19343, 0.845502, 0.967378, 1.04694, 0.692875, 0.202823, 0.1346, 0.174675, 0.420345, 0.590907, 0.589285, 0.677052, 0.869995, 1.0262, 1.06891, 1.05181, 1.01596, 0.977932, 0.950848, 0.943689, 0.97363, 1.02756, 1.09454, 1.17747, 1.26312, 1.33795, 1.40683, 1.46735, 1.47471, 1.47231, 1.40595, 1.30614, 1.16262, 0.993467, 0.805627, 0.625108, 0.460575, 0.325663, 0.216804, 0.137017, 0.0839585, 0.0491912, 0.0284753, 0.015879, 0.00837953, 0.0044919, 0.00231596, 0.0011847, 0.000591002, 0.000293888, 0.000141946, 6.65161e-05, 3.25877e-05, 1.59993e-05, 7.15925e-06, 3.26838e-06, 1.4703e-06, 6.04431e-07, 3.07641e-07, 1.09555e-07, 3.85861e-08, 1.19806e-08, 3.9631e-09, 9.73442e-10, 2.95299e-10, 6.06246e-11, 1.56818e-11, 3.65536e-12, 7.51407e-13, 1.70221e-13, 3.39418e-14, 6.26074e-15, 1.42718e-15, ],
    '2016G': [ 0.081847, 0.850802, 0.932273, 0.856234, 0.930757, 0.999119, 0.588744, 0.226755, 0.150877, 0.16996, 0.340903, 0.702185, 0.967549, 1.06707, 1.03605, 1.01435, 1.00903, 0.992279, 0.94431, 0.894769, 0.882885, 0.909783, 0.96819, 1.0352, 1.10282, 1.18208, 1.26949, 1.35802, 1.44978, 1.53163, 1.54507, 1.5291, 1.42877, 1.28383, 1.09534, 0.891981, 0.688059, 0.509047, 0.359955, 0.246823, 0.161569, 0.102094, 0.0637469, 0.0388437, 0.0238807, 0.0144377, 0.00842299, 0.00508199, 0.00299714, 0.00177923, 0.00104345, 0.000617124, 0.0003582, 0.000203542, 0.000121799, 7.33852e-05, 4.03517e-05, 2.25665e-05, 1.23327e-05, 6.07475e-06, 3.63651e-06, 1.49026e-06, 5.90394e-07, 2.01768e-07, 7.20814e-08, 1.88249e-08, 5.99766e-09, 1.28107e-09, 3.42262e-10, 8.18886e-11, 1.725e-11, 3.98134e-12, 8.14704e-13, 1.56493e-13, 3.13868e-14, ],
    '2016H': [ 1.34222, 1.74825, 1.83651, 0.932058, 1.03678, 0.985467, 0.534111, 0.339354, 0.281775, 0.119824, 0.273535, 0.515269, 0.621425, 0.707113, 0.810069, 0.870013, 0.896572, 0.904524, 0.885266, 0.86234, 0.878275, 0.920587, 0.957206, 0.964682, 0.963803, 0.988797, 1.03874, 1.10295, 1.18683, 1.29253, 1.38469, 1.50545, 1.59892, 1.68474, 1.73095, 1.73405, 1.67104, 1.55786, 1.39012, 1.19619, 0.971124, 0.748378, 0.558233, 0.397027, 0.27799, 0.186737, 0.118186, 0.0756491, 0.0463797, 0.0281113, 0.0165719, 0.00972361, 0.00554147, 0.00306947, 0.00178463, 0.00104672, 0.000565044, 0.000315409, 0.000176558, 9.2329e-05, 6.1435e-05, 2.95937e-05, 1.47153e-05, 6.79693e-06, 3.55164e-06, 1.46724e-06, 7.94372e-07, 3.0651e-07, 1.55525e-07, 7.36064e-08, 3.17311e-08, 1.5437e-08, 6.83873e-09, 2.92186e-09, 1.3682e-09, ],
    }

mfvWeight = cms.EDProducer('MFVWeightProducer',
                           throw_if_no_mcstat = cms.bool(True),
                           mevent_src = cms.InputTag('mfvEvent'),
                           enable = cms.bool(True),
                           prints = cms.untracked.bool(False),
                           histos = cms.untracked.bool(True),
                           weight_gen = cms.bool(False),
                           weight_gen_sign_only = cms.bool(False),
                           weight_pileup = cms.bool(True),
                           pileup_weights = cms.vdouble(*pileup_weights[year]),
                           weight_npv = cms.bool(False),
                           npv_weights = cms.vdouble(),
                           )
