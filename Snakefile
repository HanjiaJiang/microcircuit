stp = [0, 2]
vip2som = [0, 1]
epsp = [0, 1]
ipsp = 0

localrules: all, create

rule all:
    input:
        expand('done_{a}_{b}_{c}_{d}', a=stp, b=vip2som, c=epsp, d=ipsp)
    shell:
        '''
        mkdir -p png/
        cp -r *.py *.yml Snakefile* microcircuit/ scans/ png/
        rm {input}
        '''

rule snakes:
    input:
        '{a}_{b}_{c}_{d}/done'
    output:
        'done_{a}_{b}_{c}_{d}'
    shell:
        '''
        cp -r microcircuit/ scans/ snake-gs.sh cluster.json config.yml run_network.py $(dirname {input})
        cd $(dirname {input})
        sbatch snake-gs.sh
        cd ..
        rm {input}
        > {output}
        '''

rule create:
    output:
        expand('{a}_{b}_{c}_{d}/done', a=stp, b=vip2som, c=epsp, d=ipsp)
    shell:
        '''
        python create_snakes.py {output}
        touch {output}
        '''
