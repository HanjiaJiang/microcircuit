import numpy as np
from scipy import integrate
import os
from multiprocessing import Process
from multiprocessing import Manager

def conn_barrel_integrate():
    # BBP data
    conn_arr_bbp = np.array(
        [[0.056 , 0.0751, 0.15  , 0.0397, 0.0056, 0.0148, 0.034 , 0.0009, 0.0011, 0.012 , 0.    , 0.    , 0.    ],
         [0.0589, 0.069 , 0.1211, 0.0605, 0.0031, 0.0109, 0.0233, 0.0004, 0.0012, 0.0111, 0.0001, 0.0001, 0.006 ],
         [0.028 , 0.0658, 0.11  , 0.0293, 0.003 , 0.0088, 0.033 , 0.0003, 0.0011, 0.017 , 0.    , 0.    , 0.    ],
         [0.0198, 0.0457, 0.0502, 0.0243, 0.0019, 0.0078, 0.0129, 0.0002, 0.0008, 0.0048, 0.    , 0.0001, 0.0046],

         [0.0639, 0.0261, 0.0445, 0.053 , 0.0684, 0.0678, 0.0641, 0.0124, 0.0211, 0.0389, 0.002 , 0.0016, 0.0203],
         [0.0288, 0.011 , 0.018 , 0.0389, 0.047 , 0.0697, 0.0555, 0.0066, 0.021 , 0.0336, 0.0009, 0.0019, 0.015 ],
         [0.03  , 0.0131, 0.018 , 0.0366, 0.0396, 0.0689, 0.034 , 0.0079, 0.0234, 0.035 , 0.001 , 0.0017, 0.019 ],

         [0.0783, 0.0343, 0.0768, 0.0508, 0.1028, 0.0497, 0.0566, 0.0815, 0.0739, 0.0848, 0.0124, 0.0127, 0.0383],
         [0.0149, 0.0035, 0.0088, 0.0155, 0.0287, 0.0205, 0.0167, 0.0534, 0.0305, 0.0405, 0.0027, 0.0093, 0.0194],
         [0.014 , 0.0032, 0.0084, 0.0158, 0.0312, 0.0205, 0.019 , 0.0518, 0.027 , 0.038 , 0.0035, 0.0114, 0.024 ],

         [0.0122, 0.0017, 0.0027, 0.0066, 0.0219, 0.012 , 0.008 , 0.0353, 0.0156, 0.0117, 0.0764, 0.0712, 0.0586],
         [0.0024, 0.0001, 0.0003, 0.0019, 0.0043, 0.0015, 0.0004, 0.0136, 0.    , 0.0028, 0.037 , 0.042 , 0.0593],
         [0.002 , 0.    , 0.0002, 0.001 , 0.0004, 0.0013, 0.0003, 0.0129, 0.    , 0.0017, 0.0642, 0.0517, 0.051 ]]
    )

    # raw experimental data
    # 190707
    conn_arr_exp = np.array([
        [0.1182, 0.4300, 0.6250, 0.0000, 0.1316, 0.5000, 0.0000, 0.0264, 0.1014, 0.2100, 0.0000, 0.0000, 0.0000],
        [0.5100, 0.4680, 0.2903, 0.0000, 0.1000, 0.0000, 0.0000, 0.0289, 0.2182, 0.2520, 0.0000, 0.0000, 0.0000],
        [0.3100, 0.1818, 0.0000, 0.3548, 0.0000, 0.0000, 0.0000, 0.0000, 0.0250, 0.0000, 0.0000, 0.0000, 0.0000],
        [0.0000, 0.0426, 0.4375, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],

        [0.0159, 0.0833, 0.0000, 0.0000, 0.2428, 0.6300, 0.3800, 0.0073, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.1290, 0.4800, 0.5600, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.4000, 0.6100, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],   # vip to L4 som still not sure

        [0.0947, 0.0595, 0.0811, 0.0000, 0.1044, 0.0000, 0.0000, 0.1279, 0.2500, 0.2105, 0.0115, 0.0000, 0.0000],
        [0.0794, 0.1364, 0.0500, 0.0286, 0.0000, 0.0000, 0.0000, 0.1130, 0.4776, 0.3467, 0.0000, 0.0000, 0.0000],
        [0.1124, 0.0081, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0784, 0.1579, 0.0000, 0.0000, 0.0000, 0.0000],

        [0.0000, 0.0000, 0.0000, 0.0000, 0.0323, 0.0000, 0.0000, 0.0465, 0.0000, 0.0000, 0.0282, 0.4367, 0.4052],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.2507, 0.4752, 0.3990],
        [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1428, 0.4464, 0.0510],   # vip to L6 som still not sure
    ])

    # radius of experimental data; also as flags for integration
    dia_arr_exp = np.array([
        [100.0, 141.0,  100.0,  0.0,    100.0,  141.0,  0.0,    100.0,  106.0,  106.0,  100.0,  0.0,    0.0],
        [88.00, 141.0,  100.0,  0.0,    141.0,  0.0,    0.0,    106.0,  106.0,  106.0,  0.0,    0.0,    0.0],
        [88.00, 141.0,  0.0,    141.0,  0.0,    0.0,    0.0,    0.0,    106.0,  0.0,    0.0,    0.0,    0.0],
        [0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],

        [100.0, 100.0,  0.0,    0.0,    100.0,  100.0,  100.0,  100.0,  0.0,    0.0,    100.0,  0.0,    0.0],
        [0.0,   0.0,    0.0,    0.0,    100.0,  100.0,  100.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0],
        [0.0,   0.0,    0.0,    0.0,    100.0,  100.0,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],

        [100.0, 106.0,  106.0,  0.0,    100.0,  0.0,    0.0,    100.0,  141.0,  141.0,  100.0,  0.0,    0.0],
        [106.0, 106.0,  106.0,  0.0,    0.0,    0.0,    0.0,    141.0,  141.0,  141.0,  0.0,    0.0,    0.0],
        [106.0, 106.0,  0.0,    0.0,    0.0,    0.0,    0.0,    141.0,  141.0,  0.0,    0.0,    0.0,    0.0],

        [100.0, 0.0,    0.0,    0.0,    100.0,  0.0,    0.0,    100.0,  0.0,    0.0,    100.0,  100.0,  100.0],
        [0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    100.0,  100.0,  100.0],
        [0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    100.0,  100.0,  100.0],
    ])
    # dia_arr_exp = np.zeros(arr_shape)

    # Number and Laminar Distribution of Neurons in a Thalamocortical Projection Column of Rat Vibrissal Cortex
    z_arr = np.array([
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

    arr_shape = conn_arr_exp.shape
    conn_arr_out = np.zeros(arr_shape)
    r_barrel = np.sqrt(1e5 / np.pi)
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
            r = dia_arr_exp[i, j]/2
            if r == 0:
                r = 210.0
            # 0415 integrate exp data if radius != 0
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

    for i in range(arr_len):
        for j in range(arr_len):
            pre_base = z_arr[j][0]
            post_base = z_arr[i][0]
            r = dia_arr_exp[i, j] / 2
            if r == 0:
                r = 210.0
                conn = conn_arr_bbp[i, j]
                # conn_arr_out[i, j] = conn
            else:
                # 0415
                r = 100.0
                conn = conn_arr_exp[i, j]
            f = return_dict[keygen(r, pre_base, post_base)]
            f_barrel = return_dict[keygen(r_barrel, pre_base, post_base)]
            print('post{0}, pre{1}: f={2}, f_barrel={3}'.format(i, j, f, f_barrel))
            if f != 0:
                conn_arr_out[i, j] = conn*f_barrel/f
            # print('{0}, {1}th f_barrel = {2}'.format(i, j, f_barrel))
            # conn_arr_out[i, j] = conn*f_barrel

    np.set_printoptions(precision=4, suppress=True, linewidth=200)
    print('raw exp map = \n{0}'.format(repr(conn_arr_exp)))
    print('raw exp diameter = \n{0}'.format(repr(dia_arr_exp)))
    print('bbp map = \n{0}'.format(repr(conn_arr_bbp)))
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


