localrules: all, create

rule all:
    input:
        expand('scans/{a}_{b}_{c}_{d}/raster_plot.png',
            a=range(8),
            b=range(10),
            c=6,
            d=4)
    shell:
        '''
        cp * scans/
        '''

rule create:
    output:
        expand('scans/{a}_{b}_{c}_{d}.pickle',
            a=range(8),
            b=range(10),
            c=6,
            d=4)
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
    output: 'scans/{a}_{b}_{c}_{d}/raster_plot.png'
    shell:
        '''
        python -W ignore run_network.py {input} {output}
        '''
