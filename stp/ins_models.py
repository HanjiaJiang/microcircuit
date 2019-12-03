import numpy as np
import nest
import copy
from microcircuit.network_params import net_dict
import matplotlib.pyplot as plt
from stp_dicts import cell_types, doiron_stp, test_stp
np.set_printoptions(precision=3, suppress=True)

# define
sim_time = 200.0
cell_colors = ['b', 'r', 'orange', 'g']
syn_w = 1000.0
mm_resol = 0.1
order = 1   # which input cell it is
input_type = cell_types[order]
spk_w = 50000.0
plot_Vms = True

# spike generator
def create_spks(start):
    isi = 20.0
    spike_times = np.arange(start, start + sim_time, isi)
    spk = nest.Create('spike_generator')
    nest.SetStatus(spk, {'spike_times': spike_times})
    return spk, spike_times


# multimeter
def create_mm():
    mm = nest.Create('multimeter')
    nest.SetStatus(mm, {"withtime": True, "record_from": ["V_m"], 'interval': mm_resol})
    return mm


# stp experimental data
# columns: pre-synaptic, rows: post-synaptic
allen_1to8 = np.array([[-0.273996, -0.455301, -0.137356, -0.008266],
                       [-0.185856, -0.365362, -0.130458, -0.066462],
                       [0.165965, -0.423736, -0.147765, -0.096198],
                       [-0.007874, -0.327010, 0.132345, -0.128896]])


def set_neuron_status(neuron, subtype, n_dict):
    neuron_defaults = nest.GetDefaults(n_dict['neuron_model'])
    default_keys = [key for key, value in neuron_defaults.items()]
    set_dict = {}
    for key, value in n_dict['neuron_params'].items():
        # print(type(value), value)
        if key in default_keys:
            if type(value) is dict:
                set_dict[key] = value[subtype]
            else:
                set_dict[key] = value
    set_dict['V_m'] = n_dict['neuron_params']['E_L'][subtype]
    set_dict['V_reset'] = n_dict['neuron_params']['E_L'][subtype]
    print(set_dict)
    nest.SetStatus(neuron, set_dict)


def reshape_mm(Vms, ts, cell_n, resolution):
    freq = int(1.0/resolution)
    Vms = np.reshape(Vms, (int(len(Vms) / (cell_n*freq)), cell_n, freq))
    ts = np.reshape(ts, (int(len(ts) / (cell_n*freq)), cell_n, freq))
    Vms_new = []
    ts_new = []
    for i in range(4):
        tmp_ts = tmp_Vms = None
        for j, k in zip(ts, Vms):
            if tmp_ts is None and tmp_Vms is None:
                tmp_ts = j[i]
                tmp_Vms = k[i]
            else:
                tmp_ts = np.concatenate((tmp_ts, j[i]))
                tmp_Vms = np.concatenate((tmp_Vms, k[i]))
        ts_new.append(tmp_ts)
        Vms_new.append(tmp_Vms)
    ts_new = np.array(ts_new).T
    Vms_new = np.array(Vms_new).T
    return Vms_new, ts_new


def calc_psp_amp(spk_ts, rec_ts, rec_vms):
    psps = np.zeros(len(spk_ts) - 1)
    for i in range(len(spk_ts) - 1):
        vms = rec_vms[(rec_ts >= spk_ts[i]) & (rec_ts < spk_ts[i+1])]
        if len(vms) != 0:
            psps[i] = (np.max(vms) - np.min(vms))
    return psps


def connect_cells(g_pre, g_post, stp_dict=doiron_stp, c_types=cell_types, w=syn_w, disconnect=False):
    # set connections
    for post_type in c_types:
        for pre_type in c_types:
            # static
            syn_dict = {
                'model': 'static_synapse',
                'weight': w
            }
            if pre_type != 'Exc':
                syn_dict['weight'] = -w
            discon_spec = 'static_synapse'

            # stp
            if pre_type in stp_dict:
                if post_type in stp_dict[pre_type]:
                    syn_dict = copy.deepcopy(stp_dict[pre_type][post_type])
                    # if pre_type == 'Exc':
                    #     syn_dict['tau_psc'] = net_dict['neuron_params']['tau_syn_ex']
                    #     syn_dict['weight'] = w
                    # else:
                    #     syn_dict['tau_psc'] = net_dict['neuron_params']['tau_syn_in']
                    #     syn_dict['weight'] = -w
                    discon_spec = 'tsodyks_synapse'
            # print(syn_dict)
            if disconnect is False:
                nest.Connect(g_pre[pre_type], g_post[post_type], syn_spec=syn_dict)
            else:
                nest.Disconnect(g_pre[pre_type], g_post[post_type], conn_spec={'rule': 'one_to_one'}, syn_spec={'model': discon_spec})


def simulate(dev, g_post, s_time=sim_time):
    # simulation
    for key, value in g_post.items():
        nest.Connect(dev, value)
    nest.Simulate(s_time)
    dmm = nest.GetStatus(dev)[0]
    Vms = dmm["events"]["V_m"]
    ts = dmm["events"]["times"]
    return Vms, ts


def test_stp_dicts(g_pre, g_post, s_dict=doiron_stp, post_types=cell_types, pre_type=input_type, resol=mm_resol, w=syn_w):
    dt = 1.0
    min_tau = 1.0
    taus = np.arange(min_tau, 50.0, dt)
    stp_factor = np.full((len(post_types) + 1, len(taus)), np.nan)
    # test with a range of tau
    for j, tau in enumerate(taus):
        # temporary stp for test
        tmp_stp = copy.deepcopy(s_dict)
        # mutate test_stp
        for post_type in post_types:
            # syn_dict = test_stp[pre_type][post_type]
            if tmp_stp[pre_type][post_type]['tau_fac'] > 0.0:
                tmp_stp[pre_type][post_type]['tau_fac'] = tau
            else:
                tmp_stp[pre_type][post_type]['tau_rec'] = tau
            if pre_type == 'Exc':
                tmp_stp[pre_type][post_type]['tau_psc'] = net_dict['neuron_params']['tau_syn_ex']
                tmp_stp[pre_type][post_type]['weight'] = w
            else:
                tmp_stp[pre_type][post_type]['tau_psc'] = net_dict['neuron_params']['tau_syn_in']
                tmp_stp[pre_type][post_type]['weight'] = -w
            print('{}=>{}: {}'.format(pre_type, post_type, tmp_stp[pre_type][post_type]))

        mm = create_mm()
        spks, spike_times = create_spks(j*sim_time + mm_resol)
        nest.Connect(spks, g_pre[pre_type], syn_spec={'weight': spk_w})
        connect_cells(g_pre, g_post, stp_dict=tmp_stp)
        Vms, ts = simulate(mm, g_post)
        Vms, ts = reshape_mm(Vms, ts, len(post_types), resol)
        if plot_Vms is True and j < 10:
            plt.plot(ts, Vms)
            plt.show()
        # calculate psps of 4 cells
        for k in range(4):
            psps = calc_psp_amp(spike_times, ts[:, k], Vms[:, k])
            diff = np.abs(((psps[4] - psps[0]) / psps[0]) - allen_1to8[k, order])
            # stp_factor[k + 1, j] = diff
            if diff < 0.015:
                stp_factor[k + 1, j] = diff
        stp_factor[0, j] = tau
        connect_cells(g_pre, g_post, stp_dict=tmp_stp, disconnect=True)
    return stp_factor


# create cells and connect to mm
cells_pre = {}
cells_post = {}
for i, cell_type in enumerate(cell_types):
    neuron_pre = nest.Create(net_dict['neuron_model'])
    neuron_post = nest.Create(net_dict['neuron_model'])
    cells_pre[cell_type] = neuron_pre
    cells_post[cell_type] = neuron_post
    nest.SetStatus(neuron_pre, {'tau_syn_ex': 0.1, 'tau_syn_in': 0.1})
    set_neuron_status(neuron_post, cell_type, net_dict)
    # print(nest.GetStatus(cells_post[cell_type]))

# Connect 1 cell to spike generator
# nest.Connect(spk, cells_pre[input_type], syn_spec={'weight': spk_w})

result = test_stp_dicts(cells_pre, cells_post, s_dict=test_stp)
print(result.T)

# connect_cells(cells_pre, cells_post)
#
# Vms, ts = simulate(mm)
#
# Vms, ts = reshape_mm(Vms, ts, len(cell_types), mm_resol)
#
# # plotting
# plt.title('presyn. type = ' + input_type)
# for i, cell_type in enumerate(cell_types):
#     x = ts[:, i]
#     y = Vms[:, i]
#     psps = calc_psp_amp(spike_times, x, y)
#     stp_factor = (psps[4]-psps[0])/psps[0]
#     print('{}, 1-to-5th relative diff(psp) = {:f}'.format(cell_type, stp_factor))
#     plt.plot(x, y, color=cell_colors[i], label='{}, stp: {:.2f} ({:.2f})'.format(cell_type, stp_factor, allen_1to8[i][0]))
# plt.vlines(spike_times, -70.0, -40.0, linestyles=':')
# plt.legend()
# plt.xlabel('time (ms)')
# plt.ylabel('V')
# plt.show()