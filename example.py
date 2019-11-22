import numpy as np
import time
import network
from network_params import net_dict
from sim_params import sim_dict
from stimulus_params import stim_dict
from functions import special_dict
from shutil import copy
import multiprocessing as mp
from conn import *

cwd = os.getcwd()

# settings
run_sim = True
on_server = False
copy_py_file = True
use_conn_file = False

sim_dict['t_sim'] = 2000.0
sim_dict['master_seed'] = 55

net_dict['K_ext'] = np.array([2000, 2000, 1000, 300,
                              2000, 2000, 1000,
                              2000, 2000, 1000,
                              2000, 2000, 1000])
net_dict['g'] = 4.0
net_dict['bg_rate'] = 4.0
net_dict['animal'] = 'rat'
net_dict['renew_conn'] = False
net_dict['conn_probs'] = np.array([[0.0872, 0.3173, 0.4612, 0.0443, 0.1056, 0.4011, 0.0374, 0.0234, 0.09  , 0.1864, 0.    , 0.    , 0.    ],
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
# net_dict['N_full'] = np.array([5096, 400, 136, 136, 4088, 224, 112, 3264, 360, 320, 4424, 224, 184])
# net_dict['rec_dev'] = ['spike_detector', 'voltmeter']

stim_dict['thalamic_input'] = True
stim_dict['th_start'] = np.arange(1500.0, sim_dict['t_sim'], 500.0)
stim_dict['orientation'] = 0.0
# stim_dict['PSP_th'] = 0.15
# stim_dict['PSP_sd'] = 0.1

special_dict['orient_tuning'] = False
special_dict['som_fac'] = True
special_dict['pv_dep'] = True
special_dict['pv2all_dep'] = True
special_dict['weak_dep'] = True

if on_server:
    cpu_ratio = 1
else:
    cpu_ratio = 0.5
sim_dict['local_num_threads'] = int(mp.cpu_count() * cpu_ratio)

if use_conn_file is True:
    if os.path.isfile(cwd + '/conn_probs.npy') is True:
        net_dict['conn_probs'] = np.load(cwd + '/conn_probs.npy')
    else:
        print('no conn_probs.npy file; using original map')

# Initialize the network and pass parameters to it.
tic = time.time()
net = network.Network(sim_dict, net_dict, stim_dict, special_dict)
toc = time.time() - tic
print("Time to initialize the network: %.2f s" % toc)

# Connect all nodes.
tic = time.time()
if run_sim:
    net.setup()
    toc = time.time() - tic
    print("Time to create the connections: %.2f s" % toc)
    # Simulate.
    tic = time.time()
    net.simulate()
    toc = time.time() - tic
    print("Time to simulate: %.2f s" % toc)

# evaluation
raster_plot_time_idx = np.array(
    # [500.0, 1500.0]
    [stim_dict['th_start'][0] - 100.0, stim_dict['th_start'][0] + 100.0]
)
fire_rate_time_idx = np.array([500.0, stim_dict['th_start'][0]])
net.evaluate(raster_plot_time_idx, fire_rate_time_idx)

if copy_py_file is True:
    file_source = os.listdir(os.getcwd())
    file_dest = sim_dict['data_path']
    for file in file_source:
        if os.path.isfile(file) and '.py' in file:
            copy(file, file_dest)