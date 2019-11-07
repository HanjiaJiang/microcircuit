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
run_sim = False
on_server = False
copy_py_file = True
custom_conn = {'file': False, 'integrate': False}

sim_dict['t_sim'] = 2000.0

net_dict['K_ext'] = np.array([2000, 2000, 1500, 500,
                              2000, 2000, 1500,
                              2000, 2000, 1500,
                              2000, 2000, 1500])
# net_dict['K_ext'] = np.array([2000, 1500, 1000, 500,
#                               2000, 1500, 1000,
#                               2000, 1500, 1000,
#                               2000, 1500, 1000])
net_dict['g'] = 4.0
net_dict['bg_rate'] = 4.0

stim_dict['thalamic_input'] = True
stim_dict['th_start'] = np.arange(1500.0, sim_dict['t_sim'], 250.0)
stim_dict['orientation'] = 0.0
# stim_dict['PSP_th'] = 0.15
# stim_dict['PSP_sd'] = 0.1

special_dict['orient_tuning'] = False
special_dict['som_fac'] = True
special_dict['pv_dep'] = True
special_dict['pv2all_dep'] = True
special_dict['weak_dep'] = False

if on_server:
    cpu_ratio = 1
else:
    cpu_ratio = 0.5
sim_dict['local_num_threads'] = int(mp.cpu_count() * cpu_ratio)

if custom_conn['file'] is True:
    if custom_conn['integrate'] is True:
        conn_barrel_integrate()
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
    [500.0, 1500.0]
    # [stim_dict['th_start'][0] - 20.0, stim_dict['th_start'][0] + 80.0]
)
fire_rate_time_idx = np.array([500.0, stim_dict['th_start'][0]])
net.evaluate(raster_plot_time_idx, fire_rate_time_idx)

if copy_py_file is True:
    file_source = os.listdir(os.getcwd())
    file_dest = sim_dict['data_path']
    for file in file_source:
        if os.path.isfile(file) and '.py' in file:
            copy(file, file_dest)
