SIM_COUNT = 32

localrules: all, create

rule all:
    input:
        expand('scans/data{a}/box_plot.png', a=range(SIM_COUNT))

rule create:
    output:
        expand('scans/para_dict{a}.pickle', a=range(SIM_COUNT))
    shell:
        '''
        python ./microcircuit/create_params.py {output}
        '''

rule simulate:
    input: 'scans/para_dict{a}.pickle'
    output: 'scans/data{a}/box_plot.png'
    shell:
        '''
        python run_network.py {input} {output}
        '''
