conn = ['5', '6-6']
stp = [2]
lyr_epsp = [0, 1]
lyr_ipsp = [0]

localrules: all, create

rule all:
    input:
        expand('done_{a}_{b}_{c}_{d}', a=conn, b=stp, c=lyr_epsp, d=lyr_ipsp)
    shell:
        '''
        > done_all
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
        directory(expand('{a}_{b}_{c}_{d}/', a=conn, b=stp, c=lyr_epsp, d=lyr_ipsp))
    shell:
        '''
        python create_snakes.py {output}
        '''
