CONN_COUNT = 2
SOM_COUNT = 3
VIP_COUNT = 3

localrules: all, create

rule all:
    input:
        expand('scans/conn{a}_som{b}_vip{c}/box_plot.png', a=range(CONN_COUNT), b=range(SOM_COUNT), c=range(VIP_COUNT))

rule create:
    output:
        expand('scans/{a}_{b}_{c}.pickle', a=range(CONN_COUNT), b=range(SOM_COUNT), c=range(VIP_COUNT))
    shell:
        '''
        python ./microcircuit/create_params.py {output}
        '''

rule simulate:
    input: 'scans/{a}_{b}_{c}.pickle'
    output: 'scans/conn{a}_som{b}_vip{c}/box_plot.png'
    shell:
        '''
        python run_network.py {input} {output}
        '''
