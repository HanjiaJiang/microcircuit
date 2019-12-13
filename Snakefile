CONN_COUNT = 2
STP_COUNT = 3
SOM_COUNT = 2
VIP_COUNT = 2

localrules: all, create

rule all:
    input:
        expand('scans/{a}_{b}_{c}_{d}/ai.dat', a=range(CONN_COUNT), b=range(STP_COUNT), c=range(SOM_COUNT), d=range(VIP_COUNT))

rule create:
    output:
        expand('scans/{a}_{b}_{c}_{d}.pickle', a=range(CONN_COUNT), b=range(STP_COUNT), c=range(SOM_COUNT), d=range(VIP_COUNT))
    shell:
        '''
        python ./microcircuit/create_params.py {output}
        '''

rule simulate:
    input: 'scans/{a}_{b}_{c}_{d}.pickle'
    output: 'scans/{a}_{b}_{c}_{d}/ai.dat'
    shell:
        '''
        python run_network.py {input} {output}
        '''
