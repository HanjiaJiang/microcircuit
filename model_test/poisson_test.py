import nest
import numpy as np
import matplotlib.pyplot as plt

def reshape_mm(Vms, ts, cell_n, resolution):
    freq = int(1.0/resolution)
    Vms = np.reshape(Vms, (int(len(Vms) / (cell_n*freq)), cell_n, freq))
    ts = np.reshape(ts, (int(len(ts) / (cell_n*freq)), cell_n, freq))
    Vms_new = []
    ts_new = []
    for i in range(cell_n):
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
    ts_new = np.array(ts_new)
    Vms_new = np.array(Vms_new)
    return Vms_new, ts_new

if __name__ == '__main__':
    # population
    pop_len = 10
    w = 1.0
    pop = nest.Create('iaf_psc_exp', pop_len)
    poisson = nest.Create('poisson_generator')
    nest.SetStatus(poisson, {'rate': 2.0})
    conn_dict = {'rule': 'all_to_all'}
    weight_dict = {
       'distribution': 'lognormal',
       'mu': w,
       'sigma': w
       }
    syn_dict_poisson = {
        'model': 'static_synapse',
         'weight': w
        }
    nest.Connect(
        poisson, pop,
        conn_spec=conn_dict,
        syn_spec=syn_dict_poisson
        )

    # multimeter and connect
    resol = 0.1
    mm = nest.Create('multimeter')
    nest.SetStatus(mm, {"withtime": True, "record_from": ["V_m"], 'interval': resol})
    nest.Connect(mm, pop)

    #
    nest.Simulate(1000)

    #
    dmm = nest.GetStatus(mm)[0]
    Vms = dmm["events"]["V_m"]
    ts = dmm["events"]["times"]
    Vms, ts = reshape_mm(Vms, ts, pop_len, resol)
    for i in range(len(Vms)):
        plt.plot(ts[i], Vms[i])
    plt.show()
