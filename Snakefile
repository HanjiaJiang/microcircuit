F1_COUNT = 2
F2_COUNT = 2

G_START = 6
G_END = 6
G_STEP = 1

BG_START = 4
BG_END = 4
BG_STEP = 1

SOM_START = 1300
SOM_END = 1700
SOM_STEP = 100

VIP_START = 580
VIP_END = 620
VIP_STEP = 10

localrules: all, create

rule all:
    input:
        expand('scans/{a}_{b}_{c}_{d}/lts_distr.npy',
            a=range(F1_COUNT),
            b=range(F2_COUNT),
            c=range(G_START, G_END + G_STEP, G_STEP),
            d=range(BG_START, BG_END + BG_STEP, BG_STEP))
    output:
        expand('{a}_ltc-distr.png', a=range(F1_COUNT))
    shell:
        '''
        python ltc_distr.py {input}
        '''

rule create:
    output:
        expand('scans/{a}_{b}_{c}_{d}.pickle',
            a=range(F1_COUNT),
            b=range(F2_COUNT),
            c=range(G_START, G_END + G_STEP, G_STEP),
            d=range(BG_START, BG_END + BG_STEP, BG_STEP))
    shell:
        '''
        python microcircuit/create_params.py {output}
        mkdir -p scans/microcircuit/
        cp microcircuit/*.py scans/microcircuit/
        mkdir -p scans/stp/
        cp stp/*.py scans/stp/
        mkdir -p scans/conn_probs/
        cp conn_probs/* scans/conn_probs/
        '''

rule simulate:
    input: 'scans/{a}_{b}_{c}_{d}.pickle'
    output: 'scans/{a}_{b}_{c}_{d}/lts_distr.npy'
    shell:
        '''
        python -W ignore run_network.py {input} {output}
        '''
