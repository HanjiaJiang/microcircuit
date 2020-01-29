STP_COUNT = 1
G_COUNT = 10
BG_COUNT = 10

localrules: all, create

rule all:
    input:
        expand('scans/{a}_{b}_{c}/ai.dat', a=range(STP_COUNT), b=range(G_COUNT), c=range(BG_COUNT))
    shell:
        '''
        cp * scans/
        '''

rule create:
    output:
        expand('scans/{a}_{b}_{c}.pickle', a=range(STP_COUNT), b=range(G_COUNT), c=range(BG_COUNT))
    shell:
        '''
        python microcircuit/create_params.py {output}
        mkdir scans/microcircuit/
        cp microcircuit/*.py scans/microcircuit/
        mkdir scans/stp/
        cp stp/*.py scans/stp/
        mkdir scans/conn_probs/
        cp conn_probs/* scans/conn_probs/
        '''

rule simulate:
    input: 'scans/{a}_{b}_{c}.pickle'
    output: 'scans/{a}_{b}_{c}/ai.dat'
    shell:
        '''
        python run_network.py {input} {output}
        '''