import numpy as np
from scipy import integrate
import os
from multiprocessing import Process
from multiprocessing import Manager
from microcircuit.raw_data import rat_dict, mouse_dict, bbp, exp, dia, allen, dia_allen


def conn_barrel_integrate(animal_dict, arr_0, arr_1, arr_2, dia):
    arr_shape = arr_1.shape
    conn_arr_out = np.zeros(arr_shape)

    r_barrel = animal_dict['radius']
    z_arr = animal_dict['thickness']

    # integration
    data_used = []
    procs = []
    return_dict = Manager().dict()
    arr_len = arr_shape[0]
    for i in range(arr_len):
        print('post group={0}'.format(i))
        for j in range(arr_len):
            pre_base = z_arr[j][0]
            pre_top = z_arr[j][1]
            post_base = z_arr[i][0]
            post_top = z_arr[i][1]
            if dia[i, j] == 0:
                r = 210.0
            else:
                r = 100.0
            pass_flag = False
            barrel_pass_flag = False

            for idx, data in enumerate(data_used):
                if data[0] == [pre_base, post_base]:
                    barrel_pass_flag = True
                    if data[1] == r:
                        pass_flag = True

            data_used.append([[pre_base, post_base], r])

            if pass_flag == 0:
                proc_1 = Process(target=integrate_by_radius,
                                 args=(r, pre_base, pre_top, post_base, post_top, return_dict))
                procs.append(proc_1)
                proc_1.start()
                if barrel_pass_flag == 0:
                    proc_2 = Process(target=integrate_by_radius,
                                    args=(r_barrel, pre_base, pre_top, post_base, post_top, return_dict))
                    procs.append(proc_2)
                    proc_2.start()

    for proc in procs:
        proc.join()

    # assign calculated values to conn_arr_out to return
    for i in range(arr_len):
        for j in range(arr_len):
            pre_base = z_arr[j][0]
            post_base = z_arr[i][0]
            r = 100.0
            if dia[i, j] == 0:
                r = 210.0
                conn = arr_0[i, j]
            elif dia[i, j] == 1:
                conn = arr_1[i, j]
            else:
                conn = arr_2[i, j]
            f = return_dict[keygen(r, pre_base, post_base)]
            f_barrel = return_dict[keygen(r_barrel, pre_base, post_base)]
            print('post{0}, pre{1}: f={2}, f_barrel={3}'.format(i, j, f, f_barrel))
            if f != 0:
                conn_arr_out[i, j] = conn*f_barrel/f

    np.set_printoptions(precision=4, suppress=True, linewidth=200)
    print('raw exp map = \n{0}'.format(repr(arr_1)))
    print('raw exp diameter = \n{0}'.format(repr(dia)))
    print('bbp map = \n{0}'.format(repr(arr_0)))
    print('combined map = \n{0}'.format(repr(conn_arr_out)))
    np.save('conn_probs.npy', conn_arr_out)
    return conn_arr_out


def keygen(r, pre_base, post_base):
    return int(pre_base*1e8 + post_base*1e4 + r)


def integrate_by_radius(r, pre_base, pre_top, post_base, post_top, return_dict):
    proc = os.getpid()
    print('start proc {0}'.format(proc))
    # print(pre_base, post_base, r)

    # integration
    lmbda = 160.0
    f_x_vol = integrate.nquad(
        integrand_connectivity_profile_exp,
        ranges=[
            [0, r],
            [0, r],
            [0, 2 * np.pi],
            [pre_base, pre_top],
            [post_base, post_top]
        ],
        args=[2 * np.pi, lmbda],
        opts={'epsrel': 1e-1}
    )[0]

    vol = integrate.nquad(
        integrand_vol,
        ranges=[
            [0, r],
            [0, r],
            [0, 2 * np.pi],
            [pre_base, pre_top],
            [post_base, post_top]
        ],
        args=[2 * np.pi],
        opts={'epsrel': 1e-1}
    )[0]

    f = f_x_vol / vol

    return_dict[keygen(r, pre_base, post_base)] = f
    print('finished proc {0}'.format(proc))

    # return f

def integrand_connectivity_profile_exp(r1, r2, phi, z1, z2, phi_max, lmbda):
    dist = np.sqrt(r1 ** 2 - 2 * r1 * r2 * np.cos(phi) + r2 ** 2 + (z2 - z1) ** 2)
    return 2 * r1* r2 * (phi_max - phi) * np.exp(-dist/lmbda)


def integrand_vol(r1, r2, phi, z1, z2, phi_max):
    return 2 * r1 * r2 * (phi_max - phi)


