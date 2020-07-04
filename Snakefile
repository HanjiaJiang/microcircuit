conn = ['5', '6-6']     # connectivity
lyr_epsp = [0, 1]      # layer-sepcific epsp
u_compen = [0, 1]      # compensate w by U
vip_conn = [1]      # adjust vip-to-som connectivity

localrules: all, create

rule all:
    input:
        expand('done_{a}_{b}_{c}_{d}', a=conn, b=lyr_epsp, c=u_compen, d=vip_conn)
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
        expand('{a}_{b}_{c}_{d}/', a=conn, b=lyr_epsp, c=u_compen, d=vip_conn)
    shell:
        '''
        python create_snakes.py {output}
        '''
