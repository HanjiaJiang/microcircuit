import numpy as np

rat_dict = {
    # rat dimension data: (Meyer, 2010, Cerebral Cortex)
    # area ~= 0.1359 mm^2
    'radius': np.sqrt(1.359*1e5 / np.pi),    # ~= 207.98 um
    'thickness': np.array([
        [1276.0, 1656.0],
        [1276.0, 1656.0],
        [1276.0, 1656.0],
        [1276.0, 1656.0],
        [1013.0, 1276.0],
        [1013.0, 1276.0],
        [1013.0, 1276.0],
        [466.0, 1013],
        [466.0, 1013],
        [466.0, 1013],
        [0.0, 466.0],
        [0.0, 466.0],
        [0.0, 466.0]
    ])
}

mouse_dict = {
    # area ~= 0.06 mm2 (Petersen, 2019, Nat. Rev. Neurosci.)
    'radius': np.sqrt((0.6*1e5) / np.pi),    # ~= 138.19 um.
    # thickness: (Lefort, Petersen, 2009, Neuron)
    'thickness': np.array([
        [736.0, 1026.0],
        [736.0, 1026.0],
        [736.0, 1026.0],
        [736.0, 1026.0],
        [566.0, 736.0],
        [566.0, 736.0],
        [566.0, 736.0],
        [264.0, 566],
        [264.0, 566],
        [264.0, 566],
        [0.0, 264.0],
        [0.0, 264.0],
        [0.0, 264.0]
    ])
}

# Blue Brain Project data (Rat S1 cortex, morphological algorithmic estimation)
bbp = np.array(
    [[0.056 , 0.0744, 0.15  , 0.0392, 0.0056, 0.0152, 0.034 , 0.0009, 0.0012, 0.012 , 0.    , 0.    , 0.    ],
     [0.0566, 0.0693, 0.1238, 0.0605, 0.0033, 0.0123, 0.0236, 0.0004, 0.0011, 0.0114, 0.0001, 0.0001, 0.0049],
     [0.028 , 0.0697, 0.11  , 0.029 , 0.003 , 0.0088, 0.033 , 0.0003, 0.0012, 0.017 , 0.    , 0.    , 0.    ],
     [0.0197, 0.0431, 0.0497, 0.0239, 0.0019, 0.0077, 0.0128, 0.0002, 0.0008, 0.0048, 0.    , 0.0001, 0.0045],
     [0.0639, 0.0264, 0.0445, 0.0525, 0.0684, 0.0678, 0.0641, 0.0124, 0.0209, 0.0389, 0.0025, 0.0018, 0.0203],
     [0.0299, 0.012 , 0.0171, 0.0398, 0.0474, 0.0697, 0.0572, 0.0071, 0.0208, 0.0354, 0.0011, 0.0023, 0.0153],
     [0.03  , 0.0136, 0.018 , 0.0363, 0.0396, 0.0656, 0.034 , 0.0079, 0.024 , 0.035 , 0.0015, 0.0018, 0.019 ],
     [0.0783, 0.0343, 0.0768, 0.0503, 0.1028, 0.0522, 0.0566, 0.0815, 0.0738, 0.0848, 0.0148, 0.0134, 0.0383],
     [0.015 , 0.0038, 0.0087, 0.0159, 0.0287, 0.0228, 0.0171, 0.0533, 0.0298, 0.0406, 0.0038, 0.0098, 0.0195],
     [0.014 , 0.0034, 0.0084, 0.0157, 0.0312, 0.0226, 0.019 , 0.0518, 0.0286, 0.038 , 0.0049, 0.0112, 0.024 ],
     [0.0127, 0.0016, 0.0028, 0.0068, 0.0235, 0.0127, 0.0082, 0.0385, 0.0148, 0.0126, 0.0899, 0.0651, 0.0606],
     [0.0025, 0.0001, 0.0002, 0.0018, 0.0045, 0.0019, 0.0004, 0.0137, 0.    , 0.0026, 0.0431, 0.038 , 0.0611],
     [0.002 , 0.    , 0.0002, 0.001 , 0.0041, 0.0017, 0.0003, 0.0129, 0.    , 0.0017, 0.0672, 0.0448, 0.051 ]]
)

# raw experimental data
# 200423
# L6 intra-layer connectivity is the average across L2/3~L5
exp = np.array([
    [0.1182, 0.4300, 0.6250, 0.0000, 0.1316, 0.0000, 0.0000, 0.0264, 0.1014, 0.2100, 0.0000, 0.0000, 0.0000],
    [0.5100, 0.4680, 0.2903, 0.0000, 0.0000, 0.0000, 0.0000, 0.0289, 0.2182, 0.2520, 0.0000, 0.0000, 0.0000],
    [0.3100, 0.5714, 0.0357, 0.3548, 0.0000, 0.0000, 0.0000, 0.0000, 0.0250, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],

    [0.0159, 0.0000, 0.0000, 0.0000, 0.2428, 0.6300, 0.3800, 0.0073, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.1290, 0.4800, 0.5600, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.1250, 0.6100, 0.0385, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],

    [0.0947, 0.0595, 0.0811, 0.0000, 0.1044, 0.0000, 0.0000, 0.1279, 0.2500, 0.2105, 0.0115, 0.0000, 0.0000],
    [0.0794, 0.1364, 0.0500, 0.0000, 0.0000, 0.0000, 0.0000, 0.1130, 0.4776, 0.3467, 0.0000, 0.0000, 0.0000],
    [0.1124, 0.0081, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0784, 0.1579, 0.0498, 0.0000, 0.0000, 0.0000],

    [0.0000, 0.0000, 0.0000, 0.0000, 0.0323, 0.0000, 0.0000, 0.0465, 0.0000, 0.0000, 0.0282, 0.4367, 0.4052],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.2507, 0.4752, 0.3990],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1711, 0.4464, 0.0413]
])

# not using
allen = np.array([
    [0.0000, 0.3469, 0.2075, 0.0208, 0.0000, 0.4375, 0.2727, 0.0000, 0.0000, 0.0667, 0.0000, 0.0000, 0.0000],
    [0.3800, 0.3711, 0.1563, 0.0213, 0.2500, 0.2105, 0.0000, 0.0357, 0.0588, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.2712, 0.0735, 0.0357, 0.1481, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.1765, 0.0000, 0.1481, 0.0000, 0.0455, 0.0000, 0.1667, 0.0000, 0.0000, 0.1667, 0.0000, 0.0000, 0.0000],

    [0.0000, 0.1250, 0.0435, 0.0000, 0.0000, 0.2059, 0.3077, 0.0000, 0.2500, 0.0526, 0.0000, 0.0000, 0.0000],
    [0.3125, 0.1765, 0.0000, 0.0000, 0.1290, 0.4727, 0.0000, 0.0000, 0.2400, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0909, 0.0000, 0.0000, 0.1667, 0.0000, 0.0000, 0.0000, 0.0769, 0.0000, 0.0303, 0.0000, 0.0000, 0.0000],

    [0.0000, 0.0357, 0.0811, 0.0000, 0.0000, 0.1818, 0.0000, 0.0000, 0.1389, 0.1928, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.1250, 0.0500, 0.0909, 0.2308, 0.1786, 0.0000, 0.0946, 0.2143, 0.1071, 0.0000, 0.1000, 0.1538],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0294, 0.1158, 0.0862, 0.0498, 0.0000, 0.0357, 0.0000],

    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.2319, 0.0441],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0968, 0.1034, 0.1447, 0.3475, 0.1235],
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0408, 0.1299, 0.0122, 0.0455]
])

# Flags for integration (previously diameters)
dia = np.array([
    [1, 1, 1, 0,  1, 0, 0,  1, 1, 1,  1, 0, 0],
    [1, 1, 1, 0,  0, 0, 0,  1, 1, 1,  0, 0, 0],
    [1, 1, 1, 1,  0, 0, 0,  0, 1, 0,  0, 0, 0],
    [0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0],

    [1, 0, 0, 0,  1, 1, 1,  1, 0, 0,  1, 0, 0],
    [0, 0, 0, 0,  1, 1, 1,  0, 0, 0,  0, 0, 0],
    [0, 0, 0, 0,  1, 1, 1,  0, 0, 0,  0, 0, 0],

    [1, 1, 1, 0,  1, 0, 0,  1, 1, 1,  1, 0, 0],
    [1, 1, 1, 0,  0, 0, 0,  1, 1, 1,  0, 0, 0],
    [1, 1, 0, 0,  0, 0, 0,  1, 1, 1,  0, 0, 0],

    [1, 0, 0, 0,  1, 0, 0,  1, 0, 0,  1, 1, 1],
    [0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  1, 1, 1],
    [0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  1, 1, 1]
])

dia_allen = np.array([
    [1, 1, 1, 2,  1, 2, 0,  1, 1, 2,  1, 2, 0],
    [1, 2, 1, 2,  2, 2, 0,  2, 2, 1,  2, 0, 0],
    [1, 1, 2, 2,  0, 0, 0,  0, 1, 2,  0, 0, 0],
    [2, 2, 2, 2,  2, 2, 2,  2, 2, 2,  0, 0, 0],

    [1, 2, 0, 2,  1, 1, 1,  1, 2, 2,  1, 0, 0],
    [2, 2, 0, 2,  1, 1, 1,  2, 2, 0,  0, 0, 0],
    [2, 0, 0, 2,  1, 1, 1,  2, 0, 2,  0, 0, 0],

    [1, 2, 1, 2,  1, 2, 2,  1, 2, 2,  1, 0, 0],
    [1, 2, 1, 2,  2, 2, 0,  2, 2, 2,  2, 2, 2],
    [1, 1, 0, 2,  2, 0, 2,  2, 2, 2,  0, 2, 2],

    [1, 0, 0, 0,  1, 0, 0,  1, 2, 0,  1, 2, 2],
    [0, 0, 0, 0,  0, 0, 0,  0, 2, 2,  2, 2, 2],
    [0, 0, 0, 0,  0, 0, 0,  0, 2, 2,  2, 2, 2]
])
