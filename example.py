import time
import numpy as np
import network
from network_params import net_dict
from sim_params import sim_dict
from stimulus_params import stim_dict

import os
from shutil import copy
import multiprocessing as mp
from scan_params import *
from conn import *

# settings
run_sim = True
on_server = False
custom_conn = \
    {
        # use conn_probs.npy
        'file': False,
        # use conn_barrel_integrate()
        'integrate': False
    }
copy_py_file = True

# set parameters according to given data
cwd = os.getcwd()
if on_server:
    cpu_ratio = 1
else:
    cpu_ratio = 0.5
if custom_conn['file'] is True:
    if custom_conn['integrate'] is True:
        conn_barrel_integrate()
    if os.path.isfile(cwd + '/conn_probs.npy') is True:
        net_dict['conn_probs'] = np.load(cwd + '/conn_probs.npy')
    else:
        print('no conn_probs.npy file; using original map')
sim_dict['local_num_threads'] = int(mp.cpu_count() * cpu_ratio)
sim_dict['t_sim'] = 2000.0
# net_dict['K_ext'] = np.array([3000, 2600, 1200, 500,
#                                   2700, 2400, 2800,
#                                   1900, 2600, 1300,
#                                   2400, 2400, 2100])
net_dict['K_ext'] = np.array([2000, 1500, 1000, 550,
                              2000, 1500, 1000,
                              2000, 1500, 1000,
                              2000, 1500, 1000])
net_dict['g'] = g_scan
net_dict['bg_rate'] = bg_scan
stim_dict['thalamic_input'] = True
stim_dict['th_start'] = np.arange(1500.0, sim_dict['t_sim'], 250.0)
stim_dict['th_duration'] = stim_duration
stim_dict['th_rate'] = stim_rate
stim_dict['orientation'] = 0.0

# Initialize the network and pass parameters to it.
tic = time.time()
net = network.Network(sim_dict, net_dict, stim_dict)
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
    [stim_dict['th_start'][0] - 100.0, stim_dict['th_start'][0] + 100.0]
)
fire_rate_time_idx = np.array([500.0, 1500.0])
net.evaluate(raster_plot_time_idx, fire_rate_time_idx)

if copy_py_file is True:
    file_source = os.listdir(os.getcwd())
    file_dest = sim_dict['data_path']
    for file in file_source:
        if os.path.isfile(file) and '.py' in file:
            copy(file, file_dest)
