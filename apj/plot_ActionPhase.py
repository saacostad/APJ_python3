import os
import sys
from subprocess import call
from subprocess import Popen
#from math import pi, sqrt,cos,sin,atan2,fabs, asin, modf
#from subprocess import check_output
#from string import join
#import numpy as np
#import time
import argparse
import datetime

version=1
#version 1 to rename input directories. February 29, 2024
#Version 0 was copied  from /Users/jfcp/bin/plot_ActionPhase2 on February 13/2024.
print("plot_ActionPhase Version", version)

parser = argparse.ArgumentParser()
#parser.add_argument(
#        "-ip", "--ip",
#        help="ip can be 1, 2, 5 or 8",
#        choices=['1', '2', '5', '8'],
#        required=True,
#        dest="ip"
#    )
parser.add_argument(
        "-i1", "--b1_input_dirs",
        help="APJ-beam-1 tbt directories (separated by commas) that want to be plotted. Experimental data must go first",
        required=True,
        dest="b1_input_dirs"
    )
parser.add_argument(
        "-i2", "--b2_input_dirs",
        help="APJ-beam-2 tbt directories (separated by commas) that want to be plotted. Experimental data must go first",
        required=True,
        dest="b2_input_dirs"
    )
parser.add_argument(
        "-t", "--titles",
        help="Labels for every plot separated by commas",
        dest="titles"
    )
parser.add_argument(
        "-od", "--out_dir",
        help="Optional directory to store output files. By default, the program selects the directory just below b1_input_dir",
        dest="out_dir",
    )

args = parser.parse_args()
in_dirs_b1 = args.b1_input_dirs
in_dirs_b2 = args.b2_input_dirs

in_dir_b1_s = in_dirs_b1.split(',')
in_dir_b2_s = in_dirs_b2.split(',')


if (len(in_dir_b1_s) != len(in_dir_b2_s)):
    print("number of b1_input_dirs and b2_input_dirs must be equal")
    exit()


    
if (args.titles != None):
    titles=args.titles
    titles_s = titles.split(',')
    if (len(in_dir_b1_s) != len(titles_s)):
         print("number of input_dirs and titles must be equal")
         exit()
else:
    titles_s = [ ''  for i in in_dir_b1_s]

    


haf_b1 = []
hpf_b1 = []
vaf_b1 = []
vpf_b1 = []
ipaf_b1 = []
ippf_b1 = []

haf_b2 = []
hpf_b2 = []
vaf_b2 = []
vpf_b2 = []
ipaf_b2 = []
ippf_b2 = []



for i in in_dir_b1_s:
    haf = "'"+ i + 'HmaxAPJaction.sdds' + "'"
    hpf = "'"+ i + 'HmaxAPJphase.sdds' + "'"
    vaf = "'"+ i + 'VmaxAPJaction.sdds' + "'"
    vpf = "'"+ i + 'VmaxAPJphase.sdds' + "'"
    ipaf = "'"+ i + 'IPaction.sdds' + "'"
    ippf = "'"+ i + 'IPphase.sdds' + "'"
    haf_b1.append(haf)
    hpf_b1.append(hpf)
    vaf_b1.append(vaf)
    vpf_b1.append(vpf)
    ipaf_b1.append(ipaf)
    ippf_b1.append(ippf)

for i in in_dir_b2_s:
    haf = "'"+ i + 'HmaxAPJaction.sdds' + "'"
    hpf = "'"+ i + 'HmaxAPJphase.sdds' + "'"
    vaf = "'"+ i + 'VmaxAPJaction.sdds' + "'"
    vpf = "'"+ i + 'VmaxAPJphase.sdds' + "'"
    ipaf = "'"+ i + 'IPaction.sdds' + "'"
    ippf = "'"+ i + 'IPphase.sdds' + "'"
    haf_b2.append(haf)
    hpf_b2.append(hpf)
    vaf_b2.append(vaf)
    vpf_b2.append(vpf)
    ipaf_b2.append(ipaf)
    ippf_b2.append(ippf)

ip_action_b1 = ipaf_b1[0].replace("'","")
ip_action_b2 = ipaf_b2[0].replace("'","")
ip_phase_b1 = ippf_b1[0].replace("'","")
ip_phase_b2 = ippf_b2[0].replace("'","")





colors = ['red','blue','black','brown','orange', 'yellow', 'green']
start_dir = os.getcwd()
os.chdir(in_dir_b1_s[0])
full_dir_b1=os.getcwd()
os.chdir(start_dir)
full_dir_b1_s = full_dir_b1.split('/')
full_dir_b1_s.pop(0)
full_dir_b1_s.pop()
full_dir_b1_s.pop()
main_dir = '/'
for i in  full_dir_b1_s:
    main_dir = main_dir + i + '/'


action_gpl_file = open(main_dir + 'gb1b2APJaction.gpl','w')
phase_gpl_file = open(main_dir + 'gb1b2APJphase.gpl','w') 
print('reset', file=action_gpl_file)
print('set multiplot layout 2, 2', file=action_gpl_file)


print('set title "APJ actionx beam1" font ",14"', file=action_gpl_file)
print('set xlabel \'s(km)\'', file=action_gpl_file)
print('set ylabel \'APJ actionx(nm)\'', file=action_gpl_file)
for i in in_dir_b1_s:
    if(in_dir_b1_s.index(i) == 0):
        if (os.path.exists(ip_action_b1)):
            print('p',haf_b1[0], "u ($3/1000):($1==0?($4*1e9):NaN) w p pt 7 ps 0.6 notitle,",ipaf_b1[0],"u ($3/1000):($1==0?($4*1e9):NaN) w lp lt 1  pt 7 ps 1 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=action_gpl_file)
        else:
            print('p',haf_b1[0], "u ($3/1000):($1==0?($4*1e9):NaN) w p pt 7 ps 0.6 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=action_gpl_file)
    else:
        if (i == in_dir_b1_s[-1] ):
            print(","  ,haf_b1[in_dir_b1_s.index(i)], "u ($3/1000):($1==0?($4*1e9):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," notitle,",ipaf_b1[in_dir_b1_s.index(i)],"u ($3/1000):($1==0?($4*1e9):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b1_s.index(i)] +"'", file=action_gpl_file)
        else:
            print(","  ,haf_b1[in_dir_b1_s.index(i)], "u ($3/1000):($1==0?($4*1e9):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," notitle,",ipaf_b1[in_dir_b1_s.index(i)],"u ($3/1000):($1==0?($4*1e9):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b1_s.index(i)] +"'", end=' ', file=action_gpl_file)
print('\n', file=action_gpl_file)
print('unset key', file=action_gpl_file)
print('\n', file=action_gpl_file)



print('set title "APJ actiony beam1" font ",14"', file=action_gpl_file)
print('set xlabel \'s(km)\'', file=action_gpl_file)
print('set ylabel \'APJ actiony(nm)\'', file=action_gpl_file)
for i in in_dir_b1_s:
    if(in_dir_b1_s.index(i) == 0):
        if (os.path.exists(ip_action_b1)):
            print('p',vaf_b1[0], "u ($3/1000):($1==1?($4*1e9):NaN) w p pt 7 ps 0.6 notitle,",ipaf_b1[0],"u ($3/1000):($1==1?($4*1e9):NaN) w lp lt 1  pt 7 ps 1 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=action_gpl_file)
        else:
            print('p',vaf_b1[0], "u ($3/1000):($1==1?($4*1e9):NaN) w p pt 7 ps 0.6 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=action_gpl_file)
    else:
        if (i == in_dir_b1_s[-1] ):
            print(","  ,vaf_b1[in_dir_b1_s.index(i)], "u ($3/1000):($1==1?($4*1e9):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," notitle,",ipaf_b1[in_dir_b1_s.index(i)],"u ($3/1000):($1==1?($4*1e9):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b1_s.index(i)] +"'", file=action_gpl_file)
        else:
            print(","  ,vaf_b1[in_dir_b1_s.index(i)], "u ($3/1000):($1==1?($4*1e9):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," notitle,",ipaf_b1[in_dir_b1_s.index(i)],"u ($3/1000):($1==1?($4*1e9):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b1_s.index(i)] +"'", end=' ', file=action_gpl_file)
print('\n', file=action_gpl_file)

print('set title "APJ actionx beam2" font ",14"', file=action_gpl_file)
print('set xlabel \'s(km)\'', file=action_gpl_file)
print('set ylabel \'APJ actionx(nm)\'', file=action_gpl_file)
for i in in_dir_b2_s:
    if(in_dir_b2_s.index(i) == 0):
        if (os.path.exists(ip_action_b2)):
            print('p',haf_b2[0], "u ($3/1000):($1==0?($4*1e9):NaN) w p pt 7 ps 0.6 notitle,",ipaf_b2[0],"u ($3/1000):($1==0?($4*1e9):NaN) w lp lt 1  pt 7 ps 1 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=action_gpl_file)
        else:
            print('p',haf_b2[0], "u ($3/1000):($1==0?($4*1e9):NaN) w p pt 7 ps 0.6 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=action_gpl_file)
    else:
        if (i == in_dir_b2_s[-1] ):
            print(","  ,haf_b2[in_dir_b2_s.index(i)], "u ($3/1000):($1==0?($4*1e9):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," notitle,",ipaf_b2[in_dir_b2_s.index(i)],"u ($3/1000):($1==0?($4*1e9):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b2_s.index(i)] +"'", file=action_gpl_file)
        else:
            print(","  ,haf_b2[in_dir_b2_s.index(i)], "u ($3/1000):($1==0?($4*1e9):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," notitle,",ipaf_b2[in_dir_b2_s.index(i)],"u ($3/1000):($1==0?($4*1e9):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b2_s.index(i)] +"'", end=' ', file=action_gpl_file)
print('\n', file=action_gpl_file)




print('set title "APJ actiony beam2" font ",14"', file=action_gpl_file)
print('set xlabel \'s(km)\'', file=action_gpl_file)
print('set ylabel \'APJ actiony(nm)\'', file=action_gpl_file)
for i in in_dir_b2_s:
    if(in_dir_b2_s.index(i) == 0):
        if (os.path.exists(ip_action_b2)):
            print('p',vaf_b2[0], "u ($3/1000):($1==1?($4*1e9):NaN) w p pt 7 ps 0.6 notitle,",ipaf_b2[0],"u ($3/1000):($1==1?($4*1e9):NaN) w lp lt 1  pt 7 ps 1 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=action_gpl_file)
        else:
            print('p',vaf_b2[0], "u ($3/1000):($1==1?($4*1e9):NaN) w p pt 7 ps 0.6 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=action_gpl_file)
    else:
        if (i == in_dir_b2_s[-1] ):
            print(","  ,vaf_b2[in_dir_b2_s.index(i)], "u ($3/1000):($1==1?($4*1e9):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," notitle,",ipaf_b2[in_dir_b2_s.index(i)],"u ($3/1000):($1==1?($4*1e9):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b2_s.index(i)] +"'", file=action_gpl_file)
        else:
            print(","  ,vaf_b2[in_dir_b2_s.index(i)], "u ($3/1000):($1==1?($4*1e9):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," notitle,",ipaf_b2[in_dir_b2_s.index(i)],"u ($3/1000):($1==1?($4*1e9):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b2_s.index(i)] +"'", end=' ', file=action_gpl_file)
print('\n', file=action_gpl_file)








print('reset', file=phase_gpl_file)
print('set multiplot layout 2, 2', file=phase_gpl_file)

print('set title "APJ phasex beam1" font ",14"', file=phase_gpl_file)
print('set xlabel \'s(km)\'', file=phase_gpl_file)
print('set ylabel \'APJ phasex(rad)\'', file=phase_gpl_file)
for i in in_dir_b1_s:
    if(in_dir_b1_s.index(i) == 0):
        if (os.path.exists(ip_phase_b1)):
            print('p',hpf_b1[0], "u ($3/1000):($1==0?($4):NaN) w p pt 7 ps 0.6 notitle,",ippf_b1[0],"u ($3/1000):($1==0?($4):NaN) w lp lt 1  pt 7 ps 1 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=phase_gpl_file)
        else:
            print('p',hpf_b1[0], "u ($3/1000):($1==0?($4):NaN) w p pt 7 ps 0.6 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=phase_gpl_file)
    else:
        if (i == in_dir_b1_s[-1] ):
            print(","  ,hpf_b1[in_dir_b1_s.index(i)], "u ($3/1000):($1==0?($4):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," notitle,",ippf_b1[in_dir_b1_s.index(i)],"u ($3/1000):($1==0?($4):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b1_s.index(i)] +"'", file=phase_gpl_file)
        else:
            print(","  ,hpf_b1[in_dir_b1_s.index(i)], "u ($3/1000):($1==0?($4):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," notitle,",ippf_b1[in_dir_b1_s.index(i)],"u ($3/1000):($1==0?($4):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b1_s.index(i)] +"'", end=' ', file=phase_gpl_file)
print('\n', file=phase_gpl_file)
print('unset key', file=phase_gpl_file)
print('\n', file=phase_gpl_file)


print('set title "APJ phasey beam1" font ",14"', file=phase_gpl_file)
print('set xlabel \'s(km)\'', file=phase_gpl_file)
print('set ylabel \'APJ phasey(rad)\'', file=phase_gpl_file)
for i in in_dir_b1_s:
    if(in_dir_b1_s.index(i) == 0):
        if (os.path.exists(ip_phase_b1)):
            print('p',vpf_b1[0], "u ($3/1000):($1==1?($4):NaN) w p pt 7 ps 0.6 notitle,",ippf_b1[0],"u ($3/1000):($1==1?($4):NaN) w lp lt 1  pt 7 ps 1 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=phase_gpl_file)
        else:
            print('p',vpf_b1[0], "u ($3/1000):($1==1?($4):NaN) w p pt 7 ps 0.6 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=phase_gpl_file)
    else:
        if (i == in_dir_b1_s[-1] ):
            print(","  ,vpf_b1[in_dir_b1_s.index(i)], "u ($3/1000):($1==1?($4):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," notitle,",ippf_b1[in_dir_b1_s.index(i)],"u ($3/1000):($1==1?($4):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b1_s.index(i)] +"'", file=phase_gpl_file)
        else:
            print(","  ,vpf_b1[in_dir_b1_s.index(i)], "u ($3/1000):($1==1?($4):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," notitle,",ippf_b1[in_dir_b1_s.index(i)],"u ($3/1000):($1==1?($4):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b1_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b1_s.index(i)] +"'", end=' ', file=phase_gpl_file)
print('\n', file=phase_gpl_file)


print('set title "APJ phasex beam2" font ",14"', file=phase_gpl_file)
print('set xlabel \'s(km)\'', file=phase_gpl_file)
print('set ylabel \'APJ phasex(rad)\'', file=phase_gpl_file)
for i in in_dir_b2_s:
    if(in_dir_b2_s.index(i) == 0):
        if (os.path.exists(ip_phase_b2)):
            print('p',hpf_b2[0], "u ($3/1000):($1==0?($4):NaN) w p pt 7 ps 0.6 notitle,",ippf_b2[0],"u ($3/1000):($1==0?($4):NaN) w lp lt 1  pt 7 ps 1 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=phase_gpl_file)
        else:
            print('p',hpf_b2[0], "u ($3/1000):($1==0?($4*1e9):NaN) w p pt 7 ps 0.6 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=phase_gpl_file)
    else:
        if (i == in_dir_b2_s[-1] ):
            print(","  ,hpf_b2[in_dir_b2_s.index(i)], "u ($3/1000):($1==0?($4):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," notitle,",ippf_b2[in_dir_b2_s.index(i)],"u ($3/1000):($1==0?($4):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b2_s.index(i)] +"'", file=phase_gpl_file)
        else:
            print(","  ,hpf_b2[in_dir_b2_s.index(i)], "u ($3/1000):($1==0?($4):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," notitle,",ippf_b2[in_dir_b2_s.index(i)],"u ($3/1000):($1==0?($4):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b2_s.index(i)] +"'", end=' ', file=phase_gpl_file)
print('\n', file=phase_gpl_file)


print('set title "APJ phasey beam2" font ",14"', file=phase_gpl_file)
print('set xlabel \'s(km)\'', file=phase_gpl_file)
print('set ylabel \'APJ phasey(rad)\'', file=phase_gpl_file)
for i in in_dir_b2_s:
    if(in_dir_b2_s.index(i) == 0):
        
        print('p',vpf_b2[0], "u ($3/1000):($1==1?($4):NaN) w p pt 7 ps 0.6 notitle,",ippf_b2[0],"u ($3/1000):($1==1?($4):NaN) w lp lt 1  pt 7 ps 1 lc rgb", "'"+colors[0]+"'" ,"title", "'" + titles_s[0] + "'", end=' ', file=phase_gpl_file)
    else:
        if (i == in_dir_b2_s[-1] ):
            print(","  ,vpf_b2[in_dir_b2_s.index(i)], "u ($3/1000):($1==1?($4):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," notitle,",ippf_b2[in_dir_b2_s.index(i)],"u ($3/1000):($1==1?($4):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b2_s.index(i)] +"'", file=phase_gpl_file)
        else:
            print(","  ,vpf_b2[in_dir_b2_s.index(i)], "u ($3/1000):($1==1?($4):NaN) w l lw 2 lt 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," notitle,",ippf_b2[in_dir_b2_s.index(i)],"u ($3/1000):($1==1?($4):NaN) w lp lt 1 pt 7 ps 1 lc rgb ", "'"+colors[in_dir_b2_s.index(i)]+"'" ," title ", "'" + titles_s[in_dir_b2_s.index(i)] +"'", end=' ', file=phase_gpl_file)

action_gpl_file.close() 
phase_gpl_file.close()

            
gaction = main_dir + 'gb1b2APJaction.gpl'
gphase = main_dir + 'gb1b2APJphase.gpl'
call(['gnuplot', '-p', gaction])            
call(['gnuplot', '-p', gphase])

version_file = open(main_dir + 'version_plot_ActionPhase.dat','a')
print("Version", version, "timestamp", datetime.datetime.now(), file=version_file)
version_file.close()




