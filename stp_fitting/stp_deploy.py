import os
import json

if __name__ == '__main__':
    with open('stp_exp-data.json', 'r') as j:
        exp_dicts = json.load(j)
    with open('Snakefile', 'r') as s:
        snake_text = s.read()
    for k, v in exp_dicts.items():
        os.system('mkdir -p {}'.format(k))
        constants = 'article=\'{}\'\npretype=\'{}\'\nposttype=\'{}\'\nspk_isi=\'{}\'\npprs=\'{}\'\npeaks=\'{}\'\n'.format(v['article'], v['pre_subtype'], v['post_subtype'], v['spk_isi'], v['pprs'], v['peaks'])
        snake_tmp = constants + snake_text
        with open('{}/Snakefile'.format(k), 'w') as s:
            s.write(snake_tmp)
        os.system('cp stp*.py *.json *.yml snake-stp.sh {}/'.format(k))
        # os.system('sbatch {}/snake-stp.sh'.format(k))
