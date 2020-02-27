F1_COUNT = 2
F2_COUNT = 2
G_START = 4
G_COUNT = 2
BG_START = 4
BG_COUNT = 2

localrules: all, create

rule all:
    input:
        expand('scans/{f1}_{f2}_{g}_{bg}/ai.dat', f1=range(F1_COUNT), f2=range(F2_COUNT), g=range(G_START, G_START + G_COUNT), bg=range(BG_START, BG_START + BG_COUNT))
    shell:
        '''
        cp * scans/
        '''

rule create:
    output:
        expand('scans/{f1}_{f2}_{g}_{bg}.pickle', f1=range(F1_COUNT), f2=range(F2_COUNT), g=range(G_START, G_START + G_COUNT), bg=range(BG_START, BG_START + BG_COUNT))
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
    input: 'scans/{f1}_{f2}_{g}_{bg}.pickle'
    output: 'scans/{f1}_{f2}_{g}_{bg}/ai.dat'
    shell:
        '''
        python run_network.py {input} {output}
        '''