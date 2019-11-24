SIM_COUNT = 10

localrules: all, create

rule all:
    input:
        expand('sim{a}/data/box_plot.png', a=range(SIM_COUNT))

rule create:
    output:
        expand('sim{a}/para_dict.pickle', a=range(SIM_COUNT))
    shell:
        '''
        python ./microcircuit/create_params.py {output}
        '''

rule simulate:
    input: 'sim{a}/para_dict.pickle'
    output: 'sim{a}/data/box_plot.png'
    shell:
        '''
        python run_network.py {input} {output}
        '''
