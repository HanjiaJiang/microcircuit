F1_COUNT = 2
F2_COUNT = 2
G_START = 4
G_END = 10
G_STEP = 1
BG_START = 2
BG_END = 8
BG_STEP = 1
SOM_START = 1200
SOM_END = 1800
SOM_STEP = 100
VIP_START = 500
VIP_END = 620
VIP_STEP = 20

localrules: all, create

rule all:
    input:
        expand('scans/{a}_{b}_{c}_{d}/ai.dat', a=range(F1_COUNT), b=range(F2_COUNT), c=range(SOM_START, SOM_END + SOM_STEP, SOM_STEP), d=range(VIP_START, VIP_END + VIP_STEP, VIP_STEP))
    shell:
        '''
        cp * scans/
        '''

rule create:
    output:
        expand('scans/{a}_{b}_{c}_{d}.pickle', a=range(F1_COUNT), b=range(F2_COUNT), c=range(G_START, G_END + G_STEP, G_STEP), d=range(BG_START, BG_END + BG_STEP, BG_STEP))
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
    output: 'scans/{a}_{b}_{c}_{d}/ai.dat'
    shell:
        '''
        python run_network.py {input} {output}
        '''