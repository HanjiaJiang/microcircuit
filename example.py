# -*- coding: utf-8 -*-
#
# example.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

'''
pynest microcircuit example
---------------------------

Example file to run the microcircuit.

Hendrik Rothe, Hannah Bos, Sacha van Albada; May 2016

This example uses the function GetNodes, which is deprecated. A deprecation
warning is therefore issued. For details about deprecated functions, see
documentation.
'''

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

# HJ
run_sim = True
on_server = False
conn_probs_from_file = False
copy_file = True
cwd = os.getcwd()

# for parameter scan
if on_server:
    cpu_ratio = 1
else:
    cpu_ratio = 0.5
if conn_probs_from_file is True \
        and os.path.isfile(cwd+'/conn_probs.npy') is True:
    net_dict['conn_probs'] = np.load(cwd+'/conn_probs.npy')
sim_dict['local_num_threads'] = int(mp.cpu_count()*cpu_ratio)
sim_dict['t_sim'] = 2000.0

net_dict['K_ext'] = np.array([3000, 2600, 1200, 500,
                              2700, 2400, 2800,
                              1900, 2600, 1300,
                              2400, 2400, 2100])
net_dict['g'] = g_scan
net_dict['bg_rate'] = bg_scan
stim_dict['thalamic_input'] = True
stim_dict['th_start'] = np.arange(1000.0, sim_dict['t_sim'], 1000.0)
stim_dict['th_duration'] = stim_duration
stim_dict['th_rate'] = stim_rate
stim_dict['orientation'] = stim_orient

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

raster_plot_time_idx = np.array(
     [stim_dict['th_start'][0]-100.0, stim_dict['th_start'][0]+100.0]
    )
fire_rate_time_idx = np.array([500.0, 1000.0])
net.evaluate(raster_plot_time_idx, fire_rate_time_idx)

if copy_file is True:
    file_source = os.listdir(os.getcwd())
    file_dest = sim_dict['data_path']
    for file in file_source:
        if os.path.isfile(file) and '.py' in file:
            copy(file, file_dest)
