conn = ['5', '6-6']     # connectivity
lyr_epsp = [0, 1]      # layer-sepcific epsp
u_compen = [0, 1]      # compensate w by U
stp = [1, 2]

localrules: all, create

rule all:
    input:
        expand('done_{a}_{b}_{c}_{d}', a=conn, b=lyr_epsp, c=u_compen, d=stp)
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
        sbatch snake-gs.sh
        cd ..
        > {output}
        '''

rule create:
    output:
        directory(expand('{a}_{b}_{c}_{d}/', a=conn, b=lyr_epsp, c=u_compen, d=stp))
    shell:
        '''
        python create_snakes.py {output}
        '''
