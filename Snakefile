exc = 750
pv = 1250
som = range(500, 1001, 250)
vip = range(750, 1251, 250)

localrules: all, create

rule all:
    input:
        expand('done_{a}_{b}_{c}_{d}', a=exc, b=pv, c=som, d=vip)
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
        directory(expand('{a}_{b}_{c}_{d}/', a=exc, b=pv, c=som, d=vip))
    shell:
        '''
        python create_snakes.py {output}
        '''
