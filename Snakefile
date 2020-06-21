conn = ['6-6', '6-7']
ctsp = [1]
stp = [2]
g = [6, 8]

localrules: all, create

rule all:
    input:
        expand('{a}_{b}_{c}_{d}_done', a=conn, b=ctsp, c=stp, d=g)
    shell:
        '''
        > all-done
        '''

rule snakes:
    input:
        '{a}_{b}_{c}_{d}/'
    output:
        '{a}_{b}_{c}_{d}_done'
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
        directory(expand('{a}_{b}_{c}_{d}/', a=conn, b=ctsp, c=stp, d=g))
    shell:
        '''
        python create_snakes.py {output}
        '''
