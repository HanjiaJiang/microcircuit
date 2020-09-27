stp = [0]
vip2som = [0, 1]
epsp = 0
ipsp = 0

localrules: all, create

rule all:
    input:
        expand('done_{a}_{b}_{c}_{d}', a=stp, b=vip2som, c=epsp, d=ipsp)
    shell:
        '''
        rm done*
        mkdir -p png/
        cp -r *.py *.yml Snakefile* microcircuit/ scans/ png/
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
        directory(expand('{a}_{b}_{c}_{d}/', a=stp, b=vip2som, c=epsp, d=ipsp))
    shell:
        '''
        python create_snakes.py {output}
        '''
