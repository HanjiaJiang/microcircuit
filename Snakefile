stp = [0]
som = range(0, 1001, 250)
lyr_epsp = 0
lyr_ipsp = 0

localrules: all, create

rule all:
    input:
        expand('done_{a}_{b}_{c}_{d}', a=stp, b=som, c=lyr_epsp, d=lyr_ipsp)
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
        # snakemake
        sbatch snake-gs.sh
        cd ..
        > {output}
        '''

rule create:
    output:
        directory(expand('{a}_{b}_{c}_{d}/', a=stp, b=som, c=lyr_epsp, d=lyr_ipsp))
    shell:
        '''
        python create_snakes.py {output}
        '''
