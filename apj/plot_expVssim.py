#! /usr/bin/env python
import os
import sys
from subprocess import call
#from math import pi, sqrt,cos,sin,atan2,fabs, asin, modf
#from subprocess import check_output
#from string import join
#import numpy as np
#import time
import argparse
import datetime

version=2
#Version 2 to rename some flags, remove some options , define default path for madx and bein able to recover full python path and apj files directory. 1 March 2024. 
#Version 1 to make calls to programs in the apj/ directory not the bin/. February 16/2024
#Version 0 was copied  from /Users/jfcp/bin/plot_expVssim5 on February 13/2024.
print("plot_expVssim Version", version)

parser = argparse.ArgumentParser()
parser.add_argument(
        "-t1", "--b1_tbt_dir",
        help="APJ-beam-1  experimental tbt directory",
        required=True,
        dest="b1_tbt_dir"
    )
parser.add_argument(
        "-t2", "--b2_tbt_dir",
        help="APJ-beam-2 experimental tbt directory",
        required=True,
        dest="b2_tbt_dir"
    )
parser.add_argument(
        "-t", "--titles",
        help="Labels for every plot separated by commas",
        dest="titles"
    )
parser.add_argument(
        "-od", "--out_dir",
        help="Optional directory to store output files. By default, the program selects the directory just below b1input_dir",
        dest="out_dir",
    )
parser.add_argument(
        "-ef", "--IR_errors",
        help="Files with IR magnetic errors separated by commas. By default, the program chooses the file IR_errors_4quads.madx in the main directory",
        dest="IR_errors",
    )
parser.add_argument(
        "-mx", "--mad_program",
        help="Location of madx executable",
        dest="mad_program",
        #default='/Users/jfcp/bin/madx507'
        default='/afs/cern.ch/user/m/mad/madx/releases/5.07.00/madx-linux64-gnu'
    )

args = parser.parse_args()
b1_tbt_dir = args.b1_tbt_dir
b2_tbt_dir = args.b2_tbt_dir
IR_errors = args.IR_errors

python_path = sys.executable
programs_dir =  os.path.dirname(sys.argv[0]) + '/'

#start_dir = os.getcwd()
#os.chdir(b1_tbt_dir)
#full_dir_b1=os.getcwd()
#os.chdir(start_dir)

full_dir_b1=os.path.abspath(b1_tbt_dir)
full_dir_b1_s = full_dir_b1.split('/')
full_dir_b1_s.pop(0)
name_tbt_b1=full_dir_b1_s.pop()
b1_dir = '/'
for i in  full_dir_b1_s:
    b1_dir = b1_dir + i + '/'


#os.chdir(b2_tbt_dir)
#full_dir_b2=os.getcwd()
#os.chdir(start_dir)
full_dir_b2=os.path.abspath(b2_tbt_dir)
full_dir_b2_s = full_dir_b2.split('/')
full_dir_b2_s.pop(0)
name_tbt_b2=full_dir_b2_s.pop()
b2_dir = '/'
for i in  full_dir_b2_s:
    b2_dir = b2_dir + i + '/'


full_dir_b1_s.pop()
main_dir = '/'
for i in  full_dir_b1_s:
    main_dir = main_dir + i + '/'


inp_line_b1 = open(b1_dir + 'input_line.dat','r')
for i in inp_line_b1:
    i_s = i.split()
    if (i_s[0] == 'ip'):
        ip_b1 = '-ip=' + i_s[1]
        ip=i_s[1]
    if (i_s[0] == 'beam'): beam_b1 = '-b=' + i_s[1]
    if (i_s[0] == 'model_dir'):  model_dir_b1 = '-md=' +i_s[1]
    if (i_s[0] == 'bpm_for_avermax'): bpm_for_avermax_b1 ='-ab=' +  i_s[1]
    if (i_s[0] == 'ip_phase_bpm'): ip_phase_bpm_b1 ='-ib=' +  i_s[1]
    if (i_s[0] == 'left_kick_bpm'): left_kick_bpm_b1 ='-lb=' +  i_s[1]
    if (i_s[0] == 'right_kick_bpm'): right_kick_bpm_b1 ='-rb=' +  i_s[1]
    if (i_s[0] == 'left_arc_start'): left_arc_start_b1 ='-ls=' +  i_s[1]
    if (i_s[0] == 'right_arc_start'): right_arc_start_b1 = '-rs=' + i_s[1]
    if (i_s[0] == 'left_arc_end'): left_arc_end_b1 ='-le=' +  i_s[1]
    if (i_s[0] == 'right_arc_end'): right_arc_end_b1 ='-re=' +  i_s[1]
inp_line_b1.close()

inp_line_b2 = open(b2_dir + 'input_line.dat','r')
for i in inp_line_b2:
    i_s = i.split()
    if (i_s[0] == 'ip'): ip_b2 = '-ip=' + i_s[1]
    if (i_s[0] == 'beam'): beam_b2 = '-b=' + i_s[1]
    if (i_s[0] == 'model_dir'):  model_dir_b2 = '-md=' +i_s[1]
    if (i_s[0] == 'bpm_for_avermax'): bpm_for_avermax_b2 ='-ab=' +  i_s[1]
    if (i_s[0] == 'ip_phase_bpm'): ip_phase_bpm_b2 ='-ib=' +  i_s[1]
    if (i_s[0] == 'left_kick_bpm'): left_kick_bpm_b2 ='-lb=' +  i_s[1]
    if (i_s[0] == 'right_kick_bpm'): right_kick_bpm_b2 ='-rb=' +  i_s[1]
    if (i_s[0] == 'left_arc_start'): left_arc_start_b2 ='-ls=' +  i_s[1]
    if (i_s[0] == 'right_arc_start'): right_arc_start_b2 = '-rs=' + i_s[1]
    if (i_s[0] == 'left_arc_end'): left_arc_end_b2 ='-le=' +  i_s[1]
    if (i_s[0] == 'right_arc_end'): right_arc_end_b2 ='-re=' +  i_s[1]
inp_line_b2.close()


i1_opt = '-i1=' + b1_dir 
i2_opt = '-i2=' + b2_dir
call([python_path, programs_dir + 'get_IR_corrs.py',ip_b1,i1_opt,i2_opt])
if (args.IR_errors != None):
    err_files = args.IR_errors
else:
    err_files = main_dir + 'IR'+ip+'_errors_4quads.madx' 
err_files_s= err_files.split(',')

#    titles_s = [ ''  for i in err_files_s]
if (args.titles != None):
    titles=args.titles
    titles_s = titles.split(',')
    if (len(err_files_s) != len(titles_s)):
         print("number of err_files and titles must be equal")
         exit()
else:
    titles_s = ['simu']

    
t_opt_p='-t=exp'
i1_opt_p='-i1=' + full_dir_b1 + '/'
i2_opt_p='-i2=' + full_dir_b2 + '/'
for ii in err_files_s:
    ##############Beam1 #################################
    b1_simu_dir =main_dir + 'b1_simu'+ str(err_files_s.index(ii)) + '/' 
    b2_simu_dir =main_dir + 'b2_simu'+ str(err_files_s.index(ii)) + '/' 
    td_opt_b1= '-td=' + full_dir_b1 +'/'
    od_opt_b1= '-od='+ b1_simu_dir
    af_opt_b1 = '-af='+ full_dir_b1 +'/' + 'avermax.sdds.new'
    ef_opt = '-ef=' + ii
    mx_opt = '-mx=' +  args.mad_program
    print("")
    print("")
    print("Simulating Beam1 TBT with file errors", ii, '...')
    call([python_path, programs_dir + 'tbt_simulator.py',beam_b1,ip_b1,ef_opt,model_dir_b1,td_opt_b1, od_opt_b1,mx_opt])
    with open(b1_simu_dir + "b1_analysis.out", "wb") as b1out:  call([python_path, programs_dir + 'ActionPhaseJump.py',beam_b1,ip_b1,af_opt_b1,bpm_for_avermax_b1,model_dir_b1,left_arc_start_b1,right_arc_start_b1,left_arc_end_b1,right_arc_end_b1,left_kick_bpm_b1,right_kick_bpm_b1,od_opt_b1], stdout=b1out)   
    IPact_sw_b1 = b1_simu_dir + name_tbt_b1 + '/'+ 'IPaction_bpmsw.sdds'
    IPact_b1 = b1_simu_dir + name_tbt_b1 + '/' + 'IPaction.sdds'
    IPph_sw_b1= b1_simu_dir + name_tbt_b1 + '/' + 'IPphase_bpmsw.sdds'
    IPph_b1= b1_simu_dir + name_tbt_b1  +'/' + 'IPphase.sdds'
    call(['cp',IPact_sw_b1,IPact_b1])
    call(['cp',IPph_sw_b1,IPph_b1])
    print("APJ analysis done in",  b1_simu_dir + name_tbt_b1 + '/')


    ##############Beam2 #################################
    td_opt_b2= '-td=' + full_dir_b2 +'/'
    od_opt_b2= '-od='+ b2_simu_dir
    af_opt_b2 = '-af='+ full_dir_b2 +'/' + 'avermax.sdds.new'
    #print 'tbt_simulator5',beam_b2,ip_b2,ef_opt,model_dir_b2,td_opt_b2, od_opt_b2
    print("")
    print("")
    print("Simulating Beam2 TBT with file errors", ii, '...')
    call([python_path, programs_dir + 'tbt_simulator.py',beam_b2,ip_b2,ef_opt,model_dir_b2,td_opt_b2, od_opt_b2, mx_opt])
    with open(b2_simu_dir + "b2_analysis.out", "wb") as b2out:  call([python_path, programs_dir + 'ActionPhaseJump.py',beam_b2,ip_b2,af_opt_b2,bpm_for_avermax_b2,model_dir_b2,left_arc_start_b2,right_arc_start_b2,left_arc_end_b2,right_arc_end_b2,left_kick_bpm_b2,right_kick_bpm_b2,od_opt_b2],stdout=b2out)
    IPact_sw_b2 =  b2_simu_dir + name_tbt_b2 + '/' + 'IPaction_bpmsw.sdds'
    IPact_b2 =  b2_simu_dir + name_tbt_b2 + '/' + 'IPaction.sdds'
    IPph_sw_b2=  b2_simu_dir + name_tbt_b2 + '/'+ 'IPphase_bpmsw.sdds'
    IPph_b2=  b2_simu_dir + name_tbt_b2 + '/' + 'IPphase.sdds'
    call(['cp',IPact_sw_b2,IPact_b2])
    call(['cp',IPph_sw_b2,IPph_b2])
    print("APJ analysis done in",  b2_simu_dir + name_tbt_b2 + '/')

    ##################Plot Action and Phases both beams######################
    i1_opt_p= i1_opt_p + ','+ b1_simu_dir + name_tbt_b1 + '/'
    i2_opt_p= i2_opt_p + ','+ b2_simu_dir + name_tbt_b2 + '/'
    t_opt_p = t_opt_p + ',' + titles_s[err_files_s.index(ii)]


call([python_path, programs_dir + 'plot_ActionPhase.py',i1_opt_p,i2_opt_p,t_opt_p])


version_file = open(main_dir + 'version_plot_expVssim.dat','a')
print("Version", version, "timestamp", datetime.datetime.now(), file=version_file)
version_file.close()
