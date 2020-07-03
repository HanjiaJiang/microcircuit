conn = ['6-6']     # connectivity
vip_conn = [1]      # adjust vip-to-som connectivity
lyr_epsp = [1]      # layer-sepcific epsp
lyr_ipsp = [0]      # layer-specific ipsp

localrules: all, create

rule all:
    input:
        expand('done_{a}_{b}_{c}_{d}', a=conn, b=vip_conn, c=lyr_epsp, d=lyr_ipsp)
    shell:
        '''
        rm done*
        '''

rule snakes:
    input:
        '{a}_{b}_{c}_{d}/'
    output:
        'done_{a}_{b}_{c}_{d}'
    shell:
        '''
        cp -r microcircuit/ scans/ snake-gs.sh cluster.json config.yml run_network.py {input}
        cd {input}
        snakemake
        # sbatch snake-gs.sh
        cd ..
        > {output}
        '''

rule create:
    output:
        expand('{a}_{b}_{c}_{d}/', a=conn, b=vip_conn, c=lyr_epsp, d=lyr_ipsp)
    shell:
        '''
        python create_snakes.py {output}
        '''
