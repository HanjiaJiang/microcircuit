STP_COUNT = 1
G_START = 4
G_COUNT = 5
BG_START = 2
BG_COUNT = 5

localrules: all, create

rule all:
    input:
        expand('scans/{a}_{b}_{c}/ai.dat', a=range(STP_COUNT), b=range(G_START, G_START + G_COUNT), c=range(BG_START, BG_START + BG_COUNT))
    shell:
        '''
        cp * scans/
        '''

rule create:
    output:
        expand('scans/{a}_{b}_{c}.pickle', a=range(STP_COUNT), b=range(G_START, G_START + G_COUNT), c=range(BG_START, BG_START + BG_COUNT))
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
    input: 'scans/{a}_{b}_{c}.pickle'
    output: 'scans/{a}_{b}_{c}/ai.dat'
    shell:
        '''
        python run_network.py {input} {output}
        '''