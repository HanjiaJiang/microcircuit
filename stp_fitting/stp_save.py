import os
import pickle

if __name__ == '__main__':
    stp_dicts = {}
    for f in os.listdir():
        if f.startswith('stp') and f.endswith('.pickle'):
            with open(f, 'rb') as p:
                tmp_dict = pickle.load(p)
                pre_subtype = tmp_dict['pre_subtype'].replace('-', '_')
                post_subtype = tmp_dict['post_subtype'].replace('-', '_')
                syn_dict = tmp_dict['syn_dict']
                if pre_subtype not in stp_dicts.keys():
                    stp_dicts[pre_subtype] = {}
                stp_dicts[pre_subtype][post_subtype] = syn_dict
                p.close()
    with open('stp_fitted.pickle', 'wb') as p:
        pickle.dump(stp_dicts, p)
        p.close()
