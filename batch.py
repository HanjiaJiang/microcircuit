import os
import numpy
import random

simname = 'HJ'

# where your python code for the microcircuit model resides
workingdir = os.getcwd()
# output base directory
output_dir = workingdir + '/output/'

# python file to be executed by the queue
pyf_name = 'example.py'

# job description file
jdf_name = 'HJ_batch.jdf'

bg_list = [3.0, 3.5, 4.0]
g_list = [-4]
PV_ext_list = [1200, 1400, 1600, 1800, 2000]
SOM_ext_list = [1200, 1400, 1600, 1800, 2000]
VIP_ext_list = [400, 450, 500, 550, 600]

for bg in bg_list:
    for g in g_list:
        for PV_ext in PV_ext_list:
            for SOM_ext in SOM_ext_list:
                for VIP_ext in VIP_ext_list:

                    # output directory for this parameter combination
                    this_output_dir = 'bg={0}_g={1}_PV_ext={2}_SOM_ext={3}_VIP_ext={4}'.format(bg, g, PV_ext, SOM_ext, VIP_ext)
                    full_output_dir = output_dir + this_output_dir

                    # create directory if it doesn't exist yet
                    if this_output_dir not in os.listdir(output_dir):
                        os.system('mkdir ' + full_output_dir)
                        os.system('mkdir ' + full_output_dir + '/data')

                    os.chdir(workingdir)

                    # copy all the relevant files to the output directory
                    os.system('cp example.py ' + full_output_dir)
                    os.system('cp network.py ' + full_output_dir)
                    os.system('cp network_params.py ' + full_output_dir)
                    os.system('cp sim_params.py ' + full_output_dir)
                    os.system('cp helpers.py ' + full_output_dir)
                    os.system('cp stimulus_params.py ' + full_output_dir)
                    os.system('cp conn.py ' + full_output_dir)
                    os.system('cp functions.py ' + full_output_dir)

                    os.chdir(full_output_dir)

                    # write custom parameter files, which will be invoked in the simulation after importing the default parameters
                    f = open(full_output_dir + '/scan_params.py', 'w')
                    f.write("bg_scan = " + str(bg) + '\n')
                    f.write("g_scan = " + str(g) + '\n')
                    f.write("PV_ext_scan = " + str(PV_ext) + '\n')
                    f.write("SOM_ext_scan = " + str(SOM_ext) + '\n')
                    f.write('VIP_ext_scan = ' + str(VIP_ext) + '\n')

                    # f.write("seed = " + str(seed) + '\n')
                    f.close()

                    this_pyf_name = full_output_dir + '/' + pyf_name

                    # write job description file
                    f = open(full_output_dir + '/' + jdf_name, 'w')
                    f.write('#!/bin/bash \n')
                    # set job name
                    f.write('#SBATCH --job-name ' + simname + '_{0}_{0}_{0}_{0}_{0}'.format(bg, g, PV_ext, SOM_ext, VIP_ext) + '\n')
                    # output file to send standard output stream to
                    f.write('#SBATCH -o ./data/outfile.txt' + '\n')
                    # send error stream to same file
                    f.write('#SBATCH -e ./data/errorfile.txt' + '\n')
                    # request a total of nodes*processes_per_node for this job
                    f.write('#SBATCH --ntasks=1\n')
                    f.write('#SBATCH --cpus-per-task=48\n')
                    f.write('#SBATCH --ntasks-per-node=1\n')
                    # request processor time
                    f.write('#SBATCH --time=02:00:00\n')
                    f.write('source $HOME/.bashrc\n')
                    f.write('conda activate nest\n')
                    f.write('source $HOME/hanjia/opt/nest-lognormal-nest/bin/nest_vars.sh\n')
                    f.write('python %s\n' % this_pyf_name)
                    f.close()

                    # submit job
                    os.system('sbatch ' + jdf_name)



