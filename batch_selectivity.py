import os
import numpy
import random

simname = 'INs_selectivity'

# where your python code for the microcircuit model resides
workingdir = os.getcwd()
# output base directory
output_dir = workingdir + '/output/'

# python file to be executed by the queue
pyf_name = 'example.py'

# job description file
jdf_name = simname + '_batch.jdf'

stim_duration_list = [10.0, 20.0, 30.0, 40.0, 50.0]
stim_rate_list = [60.0, 90.0, 120.0]
som_facilitate_list = [True, False]
stim_orientation_list = ['0.0', '-np.pi/2.0']

for a, duration in enumerate(stim_duration_list):
    for b, rate in enumerate(stim_rate_list):
        for c, facilitate in enumerate(som_facilitate_list):
            for d, orientation in enumerate(stim_orientation_list):

                # output directory for this parameter combination
                this_output_dir = 'd={0}_r={1}_f={2}_o={3}'.\
                    format(duration, rate, facilitate, d)
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
                f.write("import numpy as np\n")
                f.write("bg_scan = 4.0\n")
                f.write("g_scan = -4\n")
                f.write("PV_ext_scan = 2000\n")
                f.write("SOM_ext_scan = 1500\n")
                f.write('VIP_ext_scan = 500\n')
                f.write("stim_duration = " + str(duration) + '\n')
                f.write("stim_rate = " + str(rate) + '\n')
                f.write("som_facilitate = " + str(facilitate) + '\n')
                f.write("stim_orient = " + orientation + '\n')
                f.close()

                this_pyf_name = full_output_dir + '/' + pyf_name

                # write job description file
                f = open(full_output_dir + '/' + jdf_name, 'w')
                f.write('#!/bin/bash \n')
                # set job name
                f.write('#SBATCH --job-name ' + simname + '{0}_{1}_{2}_{3}'.
                        format(duration, rate, facilitate, d) + '\n')
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



