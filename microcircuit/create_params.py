import os
import sys
import numpy as np
import pickle
from microcircuit.network_params import net_dict
from microcircuit.sim_params import sim_dict
from microcircuit.stimulus_params import stim_dict
from microcircuit.functions import special_dict
from stp.stp_dicts import allen_stp, doiron_stp, doiron_stp_weak

def set_constant():
    sim_dict['t_sim'] = 2000.0
    net_dict['g'] = 4.0
    net_dict['bg_rate'] = 4.0
    net_dict['animal'] = 'mouse'
    net_dict['renew_conn'] = False
    stim_dict['thalamic_input'] = False
    stim_dict['th_start'] = np.arange(1500.0, sim_dict['t_sim'], 500.0)
    special_dict['orient_tuning'] = False
    special_dict['stp_dict'] = doiron_stp_weak


def params_single(path):
    set_constant()
    sim_dict['master_seed'] = 55
    net_dict['K_ext'] = np.array([2000, 2000, 1500, 600,
                                  2000, 2000, 1500,
                                  2000, 2000, 1500,
                                  2000, 2000, 1500])
    # net_dict['K_ext'] = np.array([2000, 2000, 1000, 500,
    #                               2000, 2000, 1000,
    #                               2000, 2000, 1000,
    #                               2000, 2000, 1000])
    net_dict['conn_probs'] = \
        np.array([[0.0872, 0.3173, 0.4612, 0.0443, 0.1056, 0.4011, 0.0374, 0.0234, 0.09  , 0.1864, 0.    , 0.    , 0.    ],
           [0.3763, 0.3453, 0.2142, 0.0683, 0.0802, 0.0135, 0.026 , 0.0257, 0.1937, 0.2237, 0.0001, 0.0001, 0.0051],
           [0.2288, 0.4216, 0.0263, 0.2618, 0.0033, 0.0097, 0.0363, 0.0003, 0.0222, 0.018 , 0.    , 0.    , 0.    ],
           [0.0222, 0.0487, 0.0561, 0.027 , 0.0021, 0.0085, 0.0141, 0.0002, 0.0008, 0.0051, 0.    , 0.0001, 0.0047],
           [0.0128, 0.0668, 0.049 , 0.0578, 0.1764, 0.4577, 0.2761, 0.0059, 0.0229, 0.0427, 0.    , 0.0019, 0.0212],
           [0.0329, 0.0132, 0.0188, 0.0438, 0.0937, 0.3487, 0.4068, 0.0078, 0.0228, 0.0389, 0.0011, 0.0024, 0.016 ],
           [0.033 , 0.015 , 0.0198, 0.2618, 0.2906, 0.4432, 0.028 , 0.0087, 0.0263, 0.0384, 0.0016, 0.0019, 0.0198],
           [0.0841, 0.0528, 0.072 , 0.0534, 0.0844, 0.0573, 0.0621, 0.0957, 0.1871, 0.1575, 0.0094, 0.0146, 0.0418],
           [0.0705, 0.1211, 0.0444, 0.0169, 0.0315, 0.025 , 0.0188, 0.0846, 0.3574, 0.2594, 0.0041, 0.0107, 0.0213],
           [0.0998, 0.0072, 0.0089, 0.2618, 0.0343, 0.0248, 0.0209, 0.0587, 0.1182, 0.0373, 0.0054, 0.0122, 0.0262],
           [0.    , 0.0017, 0.0029, 0.007 , 0.0297, 0.0133, 0.0086, 0.0381, 0.0162, 0.0138, 0.021 , 0.3249, 0.3014],
           [0.0026, 0.0001, 0.0002, 0.0019, 0.0047, 0.002 , 0.0004, 0.015 , 0.    , 0.0028, 0.1865, 0.3535, 0.2968],
           [0.0021, 0.    , 0.0002, 0.2618, 0.0043, 0.0018, 0.0003, 0.0141, 0.    , 0.0019, 0.1955, 0.3321, 0.0307]])
    # # Allen map
    # net_dict['conn_probs'] = \
    #     np.array([[0.1014, 0.3688, 0.536 , 0.0178, 0.1174, 0.3903, 0.0429, 0.0247, 0.0948, 0.0623, 0.    , 0.    , 0.    ],
    #        [0.4374, 0.3183, 0.249 , 0.0183, 0.223 , 0.1878, 0.0298, 0.0334, 0.055 , 0.2355, 0.    , 0.0001, 0.0054],
    #        [0.2659, 0.4901, 0.0306, 0.127 , 0.0038, 0.0111, 0.0417, 0.0004, 0.0234, 0.    , 0.    , 0.    , 0.    ],
    #        [0.1514, 0.    , 0.127 , 0.    , 0.0406, 0.    , 0.1487, 0.    , 0.    , 0.1558, 0.    , 0.0001, 0.005 ],
    #        [0.0142, 0.1115, 0.0562, 0.    , 0.2058, 0.5339, 0.322 , 0.0065, 0.2232, 0.047 , 0.    , 0.0021, 0.0232],
    #        [0.2788, 0.1575, 0.0216, 0.    , 0.1093, 0.4068, 0.4746, 0.    , 0.2143, 0.0446, 0.0013, 0.0026, 0.0175],
    #        [0.0811, 0.0172, 0.0227, 0.1487, 0.339 , 0.5169, 0.0326, 0.0687, 0.0303, 0.0271, 0.0017, 0.0021, 0.0217],
    #        [0.0885, 0.0334, 0.0758, 0.    , 0.0932, 0.1623, 0.    , 0.1098, 0.1192, 0.1655, 0.0103, 0.0167, 0.0477],
    #        [0.0742, 0.1168, 0.0467, 0.085 , 0.2061, 0.1595, 0.0216, 0.0812, 0.184 , 0.0919, 0.    , 0.09  , 0.1384],
    #        [0.1051, 0.0076, 0.0098, 0.    , 0.    , 0.0285, 0.0262, 0.0994, 0.074 , 0.0428, 0.0061, 0.0321, 0.    ],
    #        [0.    , 0.0018, 0.0031, 0.0075, 0.0306, 0.0145, 0.0094, 0.0418, 0.    , 0.0157, 0.0241, 0.1985, 0.0377],
    #        [0.0028, 0.0001, 0.0002, 0.002 , 0.0052, 0.0022, 0.0005, 0.0171, 0.0871, 0.093 , 0.1238, 0.2974, 0.1057],
    #        [0.0022, 0.    , 0.0002, 0.0011, 0.0047, 0.0019, 0.0003, 0.0161, 0.    , 0.0367, 0.1112, 0.0104, 0.0389]])
    net_dict['neuron_params']['tau_syn_ex'] = 2.1
    net_dict['neuron_params']['tau_syn_inh'] = 3.2
    # net_dict['renew_conn'] = True
    stim_dict['orientation'] = 0.0
    para_dict = {
        'net_dict': net_dict,
        'sim_dict': sim_dict,
        'stim_dict': stim_dict,
        'special_dict': special_dict
    }
    with open(path, 'wb') as h:
        pickle.dump(para_dict, h)


if __name__ == "__main__":
    output_list = sys.argv[1:]
    set_constant()

    # conn_probs
    cwd = os.getcwd()
    conn_probs_list = []
    conn_folder = os.path.join(cwd, 'conn_probs')
    for file in os.listdir(conn_folder):
        if file.endswith(".npy") and 'conn_probs' in file:
            tmp = np.load(os.path.join(conn_folder, file))
            if tmp.shape == net_dict['conn_probs'].shape:
                conn_probs_list.append(tmp)
    # net_dict['conn_probs'] = np.load('conn_probs.npy')

    # lists of parameters
    som_len = vip_len = np.floor(np.sqrt(len(output_list)/len(conn_probs_list)))
    som_list = np.linspace(100.0, 2000.0, som_len)
    vip_list = np.linspace(100.0, 2000.0, vip_len)

    # to output
    idx = 0
    scan_folder = os.path.dirname(output_list[idx])
    for a, conn_probs in enumerate(conn_probs_list):
        for b, som in enumerate(som_list):
            for c, vip in enumerate(vip_list):
                sim_dict['data_path'] = os.path.join(scan_folder, 'data{}'.format(idx))
                net_dict['K_ext'] = np.array([2000, 2000, som, vip,
                                              2000, 2000, som,
                                              2000, 2000, som,
                                              2000, 2000, som])
                net_dict['conn_probs'] = conn_probs
                para_dict = {
                    'net_dict': net_dict,
                    'sim_dict': sim_dict,
                    'stim_dict': stim_dict,
                    'special_dict': special_dict
                }
                print('no. {}: K_ext(som) = {}, K_ext(som) = {}'.format(idx, som, vip))
                with open(output_list[idx], 'wb') as handle:
                    pickle.dump(para_dict, handle)

                idx += 1



