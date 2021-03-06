'''
save each connection's fitting result to a common .pickle
'''
import os
import pickle

if __name__ == '__main__':
    stp_dicts = {}
    for f in os.listdir():
        if f.startswith('stp') and f.endswith('.pickle'):
            with open(f, 'rb') as p:
                # the individual stp fitting result
                tmp_dict = pickle.load(p)
                # replace names like L4-Exc to L4_Exc; to be improved
                pre_subtype = tmp_dict['pre_subtype'].replace('-', '_')
                post_subtype = tmp_dict['post_subtype'].replace('-', '_')
                # collect
                syn_dict = tmp_dict['syn_dict']
                if pre_subtype not in stp_dicts.keys():
                    stp_dicts[pre_subtype] = {}
                stp_dicts[pre_subtype][post_subtype] = syn_dict
                p.close()
    with open('stp_fitted.pickle', 'wb') as p:
        pickle.dump(stp_dicts, p)
        p.close()
