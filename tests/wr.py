import nest
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    wr = nest.Create('weight_recorder', params={"withrport": True})
    nest.CopyModel('tsodyks_synapse', 'tsodyks_synapse_rec', {"weight_recorder": wr[0]})
    a = nest.Create('iaf_psc_exp')
    b = nest.Create('iaf_psc_exp')
    nest.SetStatus(a, {'I_e': 380.})
    nest.Connect(a, b, syn_spec='tsodyks_synapse_rec')
    voltmeter = nest.Create('voltmeter')
    nest.Connect(voltmeter, b)
    nest.Simulate(1000)
    dmm = nest.GetStatus(voltmeter)[0]
    Vms = dmm['events']['V_m']
    ts = dmm['events']['times']
    Vms = Vms + 70.
    wr_events = nest.GetStatus(wr, 'events')[0]
    senders = wr_events['senders']
    targets = wr_events['targets']
    times = wr_events['times']
    weights = wr_events['weights']
    data = list(zip(times, weights, senders, targets))
    plt.plot(ts, Vms/np.max(Vms))
    plt.plot(times, weights/np.max(weights))
    plt.show()
