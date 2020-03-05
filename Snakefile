F1_COUNT = 2
F2_COUNT = 2
G_START = 4
G_STEP = 1
G_END = 10
BG_START = 2
BG_STEP = 1
BG_END = 8

localrules: all, create

rule all:
    input:
        expand('scans/{a}_{b}_{c}_{d}/ai.dat', a=range(F1_COUNT), b=range(F2_COUNT), c=range(G_START, G_END, G_STEP), d=range(BG_START, BG_END, BG_STEP))
    shell:
        '''
        cp * scans/
        '''

rule create:
    output:
        expand('scans/{a}_{b}_{c}_{d}.pickle', a=range(F1_COUNT), b=range(F2_COUNT), c=range(G_START, G_START + G_COUNT), d=range(BG_START, BG_START + BG_COUNT))
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