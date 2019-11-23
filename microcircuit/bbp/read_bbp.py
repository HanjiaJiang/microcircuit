import numpy as np
import json
# import easygui
# fn_anatomy = easygui.fileopenbox()
# fn_layer = easygui.fileopenbox()
np.set_printoptions(precision=4, suppress=True)
fn_anatomy = 'pathways_anatomy_factsheets_simplified.json'
fn_layer = 'layer_download.json'

cell_types = {
    'EXC': ['PC', 'SS', 'SP'],
    'PV': ['NBC', 'LBC'],
    # 'PV': ['SBC', 'NBC', 'LBC'],
    'SOM': ['MC'],
    'VIP': ['BP', 'DBC']
}

layers = ['L23', 'L4', 'L5', 'L6']

dict_anatomy = {}
dict_layer = {}

# Read JSON data
if fn_anatomy:
    with open(fn_anatomy, 'r') as f:
        dict_anatomy = json.load(f)

if fn_layer:
    with open(fn_layer, 'r') as f:
        dict_layer = json.load(f)

morph_str = "No. of neurons per morphological types"

# Calculation:
# gtype: 'EXC', 'PV', 'SOM', 'VIP'
# mtype: 'PC', 'SBC', ...
# pre_layer, post_layer: 'L23', 'L4', ...
def conn_sum(pre_layer, post_layer, pre_mtypes, post_mtypes):
    prelayer_gps = dict_layer[pre_layer][morph_str]
    postlayer_gps = dict_layer[post_layer][morph_str]
    max_sum = 0
    real_sum = 0
    for pre_gp, pre_n in prelayer_gps.items():  # loop each group in data ('L23_PC', ...)
        for post_gp, post_n in postlayer_gps.items():
            for pre_mtype in pre_mtypes:    # get mtype str in each group  ('PC', 'SS', ...)
                for post_mtype in post_mtypes:
                    if pre_mtype in pre_gp and post_mtype in post_gp:  # see if the str match the data
                        # exception ...
                        if (pre_mtype == 'BP' and 'BPC' in pre_gp) or (post_mtype == 'BP' and 'BPC' in post_gp):
                            continue
                        max_sum += pre_n * post_n   # max possible connection number = pre_n * post_n
                        print('**{} in {}, {} in {}'.format(pre_mtype, pre_gp, post_mtype, post_gp))
                        for key, value in dict_anatomy.items():
                            gps_str = key.split(':')    # key => e.g. 'L23_PC:L23_PC' (connection)
                            if pre_gp == gps_str[0] and post_gp == gps_str[1]:
                                print('{} = {}, {} = {}'.format(pre_gp, gps_str[0], post_gp, gps_str[1]))
                                # real connection number = pre_n * post_n * probability
                                real_sum += pre_n * post_n * value['connection_probability'] / 100.0
    return max_sum, real_sum

len_mtx = len(layers)*len(cell_types)
conn_arr = np.zeros((len_mtx, len_mtx))
pre_i = post_i = 0
for i, pre_ly in enumerate(layers):
    for pre_gtype, pre_mtypes in cell_types.items():
        for j, post_ly in enumerate(layers):
            for post_gtype, post_mtypes in cell_types.items():
                print('\nConnect: {} {} to {} {}:'.format(pre_ly, pre_gtype, post_ly, post_gtype))
                max_sum, real_sum = conn_sum(pre_ly, post_ly, pre_mtypes, post_mtypes)
                conn_prob = float(real_sum)/max_sum
                print('max conn # = {}'.format(max_sum))
                print('real conn # = {}'.format(real_sum))
                print('conn prob = {}'.format(conn_prob))
                conn_arr[post_i, pre_i] = conn_prob
                post_i += 1
        pre_i += 1
        post_i = 0

# delete vip in L4, L5, L6
conn_arr = np.delete(conn_arr, [7, 11, 15], 0)
conn_arr = np.delete(conn_arr, [7, 11, 15], 1)

# print results
print()
print(repr(conn_arr))
np.save('bbp_conn_new.npy', conn_arr)
conn_arr_old = np.array(
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
np.save('bbp_conn_old.npy', conn_arr_old)
boo_arr = np.abs(conn_arr - conn_arr_old) > 0.0001
print()
print(repr(conn_arr_old))
print()
print(repr(boo_arr))