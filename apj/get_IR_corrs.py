import os
import sys
from subprocess import call
from math import pi, sqrt,cos,sin,atan2,fabs, asin, modf
from subprocess import check_output
# from string import join
import numpy as np
import time
import argparse
import datetime

version=2
#version 2 
#version 1 to correct bug in sbs IR_correction files. February 17/2024.
#Version 0 was copied  from /Users/jfcp/bin/get_IR_corrs on February 13/2024.
print("get_IR_corrs Version", version)

def dphi_rad(phim,phit):
    phim_dec =  float(str(phim-int(phim))[1:])*2*pi
    phit_dec =  float(str(phit-int(phit))[1:])*2*pi
    dphi = phim_dec - phit_dec
    return dphi

    


def common_corr(b1x_beam1,b1y_beam1,b1x_beam2,b1y_beam2, dkq1):
    aa=[]
    bb=[]
    for line in b1x_beam1:
        lines=line.split()
        dphim = dphi_rad(float(lines[12]),float(lines[7]))
        dphi4 = dphi_rad(float(lines[9]),float(lines[7]))
        dphi6 = dphi_rad(float(lines[11]),float(lines[7]))
        dphim=0
        dphi4=0
        dphi6=0        
        bb.append(-float(lines[0])*float(lines[1]) - float(lines[13])*dkq1)
        aa.append([float(lines[2]),float(lines[3])]) 
    for line in b1y_beam1:
        lines=line.split()
        dphim = dphi_rad(float(lines[12]),float(lines[7]))
        dphi4 = dphi_rad(float(lines[9]),float(lines[7]))
        dphi6 = dphi_rad(float(lines[11]),float(lines[7]))
        #dphim=0
        dphim=0
        dphi4=0
        dphi6=0
        bb.append(-float(lines[0])*float(lines[1]) - float(lines[13])*dkq1)
        aa.append([float(lines[2]),float(lines[3])]) 
    for line in b1x_beam2:
        lines=line.split()
        dphim = dphi_rad(float(lines[12]),float(lines[7]))
        dphi4 = dphi_rad(float(lines[9]),float(lines[7]))
        dphi6 = dphi_rad(float(lines[11]),float(lines[7]))
        #dphim=0
        dphim=0
        dphi4=0
        dphi6=0
        bb.append(-float(lines[0])*float(lines[1]) - float(lines[13])*dkq1)
        aa.append([float(lines[2]),float(lines[3])]) 
    for line in b1y_beam2:
        lines=line.split()
        dphim = dphi_rad(float(lines[12]),float(lines[7]))
        dphi4 = dphi_rad(float(lines[9]),float(lines[7]))
        dphi6 = dphi_rad(float(lines[11]),float(lines[7]))
        #dphim=0
        dphim=0
        dphi4=0
        dphi6=0
        bb.append(-float(lines[0])*float(lines[1]) - float(lines[13])*dkq1)
        aa.append([float(lines[2]),float(lines[3])]) 

    aarray = np.array(aa)
    barray = np.array(bb)
    B1pQal, B1pQbl =  np.linalg.lstsq(aarray,barray)[0]
    return B1pQal, B1pQbl

def get_equivalent_grads (dkq1, dkq2, dkq3,b1xf,b1yf):
    for line in b1xf:
        lines=line.split()
        intx_q1 = float(lines[13])
        intx_q2 = float(lines[2])
        intx_q3 =  float(lines[3])
        betx_eq = float(lines[1])
        
    b1x = -(intx_q1*dkq1 + intx_q2*dkq2 + intx_q3*dkq3)/betx_eq
    for line in b1yf:
        lines=line.split()
        inty_q1 = float(lines[13])
        inty_q2 = float(lines[2])
        inty_q3 =  float(lines[3])
        bety_eq = float(lines[1])
    b1y = -(inty_q1*dkq1 + inty_q2*dkq2 + inty_q3*dkq3)/bety_eq
    return b1x, b1y

def get_equivalent_grads_IR(dkq1l, dkq2l, dkq3l, dkq1r ,dkq2r, dkq3r,b1xf,b1yf,b1xlf,b1ylf,b1xrf,b1yrf):
    line=b1xf.readline()
    lines=line.split()
    betx_eq = float(lines[1])
    line=b1yf.readline()
    lines=line.split()
    bety_eq = float(lines[1])

    line=b1xlf.readline()
    lines=line.split()
    intx_q1l = float(lines[13])
    intx_q2l = float(lines[2])
    intx_q3l =  float(lines[3])    

    line=b1ylf.readline()
    lines=line.split()
    inty_q1l = float(lines[13])
    inty_q2l = float(lines[2])
    inty_q3l =  float(lines[3])   

    line=b1xrf.readline()
    lines=line.split()
    intx_q1r = float(lines[13])
    intx_q2r = float(lines[2])
    intx_q3r =  float(lines[3])    

    line=b1yrf.readline()
    lines=line.split()
    inty_q1r = float(lines[13])
    inty_q2r = float(lines[2])
    inty_q3r =  float(lines[3])
    
    b1x = -(intx_q1l*dkq1l + intx_q2l*dkq2l + intx_q3l*dkq3l + intx_q1r*dkq1r + intx_q2r*dkq2r + intx_q3r*dkq3r)/betx_eq
    b1y = -(inty_q1l*dkq1l + inty_q2l*dkq2l + inty_q3l*dkq3l + inty_q1r*dkq1r + inty_q2r*dkq2r + inty_q3r*dkq3r)/bety_eq

    return b1x, b1y


parser = argparse.ArgumentParser()
parser.add_argument(
        "-ip", "--ip",
        help="ip can be 1, 2, 5 or 8",
        choices=['1', '2', '5', '8'],
        required=True,
        dest="ip"
    )
parser.add_argument(
        "-i1", "--b1input_dir",
        help="Directory where APJ analysis for beam 1 is already done",
        required=True,
        dest="b1input_dir"
    )
parser.add_argument(
        "-i2", "--b2input_dir",
        help="Directory where APJ analysis for beam 2 is already done",
        required=True,
        dest="b2input_dir"
    )
parser.add_argument(
        "-od", "--out_dir",
        help="Optional directory to store output files. By default, the program selects the directory just below the directory specify by -i1",
        dest="out_dir",
    )
parser.add_argument(
        "-l1", "--dkq1l",
        help="Correction assumed in Q1L quad in 1/m^(-2)",
        dest="dkq1l",
        default=0
    )
parser.add_argument(
        "-r1", "--dkq1r",
        help="Correction assumed in Q1R quad in 1/m^(-2)",
        dest="dkq1r",
        default=0
    )
parser.add_argument(
        "-l2", "--dkq2l",
        help="Correction assumed in Q2L quad in 1/m^(-2). This correction will overwrite any calculated value",
        dest="dkq2l",
    )
parser.add_argument(
        "-r2", "--dkq2r",
        help="Correction assumed in Q2R quad in 1/m^(-2). This correction will overwrite any calculated value",
        dest="dkq2r",
    )
parser.add_argument(
        "-l3", "--dkq3l",
        help="Correction assumed in Q3L quad in 1/m^(-2). This correction will overwrite any calculated value",
        dest="dkq3l",
    )
parser.add_argument(
        "-r3", "--dkq3r",
        help="Correction assumed in Q3R quad in 1/m^(-2). This correction will overwrite any calculated value",
        dest="dkq3r",
    )




        
args = parser.parse_args()
ip = args.ip
in_dir_b1 = args.b1input_dir
in_dir_b2 = args.b2input_dir

#start_dir = os.getcwd()
#os.chdir(in_dir_b1)
#full_dir_b1=os.getcwd()
full_dir_b1=os.path.abspath(in_dir_b1)
#os.chdir(start_dir)

full_dir_b1_s = full_dir_b1.split('/')
full_dir_b1_s.pop(0)
full_dir_b1_s.pop()
apj_dir = '/'
for i in  full_dir_b1_s:
    apj_dir = apj_dir + i + '/'

if (args.out_dir != None):
    apj_dir = args.out_dir

    
dkq1l= float(args.dkq1l)
dkq1r= float(args.dkq1r)

###################Four or six correctors per IR#####################################
if (os.path.exists(in_dir_b1 + 'b1xl.dat') and os.path.exists(in_dir_b1 + 'b1xr.dat') and os.path.exists(in_dir_b1 + 'b1yl.dat') and os.path.exists(in_dir_b1 + 'b1yr.dat') and          os.path.exists(in_dir_b2 + 'b1xl.dat') and os.path.exists(in_dir_b2 + 'b1xr.dat') and os.path.exists(in_dir_b2 + 'b1yl.dat') and os.path.exists(in_dir_b2 + 'b1yr.dat')):

    b1xl_beam1 = open(in_dir_b1 + 'b1xl.dat','r')
    b1xr_beam1 = open(in_dir_b1 + 'b1xr.dat','r')
    b1yl_beam1 = open(in_dir_b1 + 'b1yl.dat','r')
    b1yr_beam1 = open(in_dir_b1 + 'b1yr.dat','r')
    b1xl_beam2 = open(in_dir_b2 + 'b1xl.dat','r')
    b1xr_beam2 = open(in_dir_b2 + 'b1xr.dat','r')
    b1yl_beam2 = open(in_dir_b2 + 'b1yl.dat','r')
    b1yr_beam2 = open(in_dir_b2 + 'b1yr.dat','r')

    B1pQal, B1pQbl= common_corr(b1xl_beam1,b1yl_beam1,b1xl_beam2,b1yl_beam2,dkq1l)
    B1pQar, B1pQbr=common_corr(b1xr_beam1,b1yr_beam1,b1xr_beam2,b1yr_beam2,dkq1r)

    if (args.dkq2l != None or args.dkq3l != None):
        if (args.dkq2l != None and args.dkq3l != None):
             B1pQal= float(args.dkq2l)
             B1pQbl= float(args.dkq3l)
        else:
            print("dkq2l and dkq3l must be given if any of the two is input.")
            exit()
            
    if (args.dkq2r != None or args.dkq3r != None):
        if (args.dkq2r != None and args.dkq3r != None):
             B1pQar= float(args.dkq2r)
             B1pQbr= float(args.dkq3r)
        else:
            print("dkq2r and dkq3r must be given if any of the two is input.")
            exit()            
    

    '''
    dkq1l= 0.0000225
    dkq1r = -0.000021
    B1pQal = 0.0000016
    B1pQar = 0.0000135
    B1pQbl = 0.0000225
    B1pQbr = -0.000021
    '''

    
    b1xl_beam1.close()
    b1xr_beam1.close()
    b1yl_beam1.close()
    b1yr_beam1.close()
    b1xl_beam2.close()
    b1xr_beam2.close()
    b1yl_beam2.close()
    b1yr_beam2.close()



    if (args.dkq1l != 0 or args.dkq1r != 0 ):
        f4q_for_apjsim=open(apj_dir +'IR'+ip+'_errors_6quads.madx','w')
        knob_file_4q= apj_dir + '6qLocalCorrectionIP'+ip+'_APJ.dat'
        f4q_for_sbssim = open(apj_dir + '6qcorrections_IP'+ip+'.madx','w')
        print("six correctors per IR")
    else:
        f4q_for_apjsim=open(apj_dir +'IR'+ip+'_errors_4quads.madx','w')
        knob_file_4q= apj_dir + '4qLocalCorrectionIP'+ip+'_APJ.dat'
        f4q_for_sbssim = open(apj_dir + '4qcorrections_IP'+ip+'.madx','w')
        print("four correctors per IR")


    print("corrq1l" + ip + " =",dkq1l,";")
    print("corrq2l" + ip + " =",B1pQal,";")
    print("corrq3l" + ip + " =",B1pQbl,";")
    print("corrq1r" + ip + " =",dkq1r,";")
    print("corrq2r" + ip + " ="  ,B1pQar,";")
    print("corrq3r" + ip + " =",B1pQbr,";")
    print(" ")       

    print("errq1l" + ip + " =",-dkq1l,";", file=f4q_for_apjsim)
    print("errq2l" + ip + " =",-B1pQal,";", file=f4q_for_apjsim)
    print("errq3l" + ip + " =",-B1pQbl,";", file=f4q_for_apjsim)
    print("errq1r" + ip + " =",-dkq1r,";", file=f4q_for_apjsim)
    print("errq2r" + ip + " ="  ,-B1pQar,";", file=f4q_for_apjsim)
    print("errq3r" + ip + " =",-B1pQbr,";", file=f4q_for_apjsim)
    f4q_for_apjsim.close()


    f4q_knob=open(knob_file_4q,'w')
    print('@TITLE' + ' 4qLocalCorrectionIP'+ip+'_APJ', file=f4q_knob)
    print('MQXA3.L'+ip+'/K1', B1pQbl, file=f4q_knob)
    print('MQXA3.R'+ip+'/K1', B1pQbr, file=f4q_knob)
    print('MQXB2.L'+ip+'/K1', B1pQal, file=f4q_knob)
    print('MQXB2.R'+ip+'/K1', B1pQar, file=f4q_knob)
    print('MQXA1.L'+ip+'/K1', dkq1l, file=f4q_knob)
    print('MQXA1.R'+ip+'/K1', dkq1r, file=f4q_knob)
    
    f4q_knob.close()

    if(-B1pQal <= 0): print('MQXB.B2L'+ip+'->K1 =  MQXB.B2L'+ip+'->K1',-B1pQal,' ;', file=f4q_for_sbssim)
    else:  print('MQXB.B2L'+ip+'->K1 =  MQXB.B2L'+ip+'->K1 +',-B1pQal,' ;', file=f4q_for_sbssim)
    if (-B1pQal <= 0):  print('MQXB.A2L'+ip+'->K1 =  MQXB.A2L'+ip+'->K1',-B1pQal,' ;', file=f4q_for_sbssim)
    else:  print('MQXB.A2L'+ip+'->K1 =  MQXB.A2L'+ip+'->K1 +',-B1pQal,' ;', file=f4q_for_sbssim)
    if (-B1pQbl <= 0): print('MQXA.3L'+ip+'->K1 =  MQXA.3L'+ip+'->K1',-B1pQbl,' ;', file=f4q_for_sbssim)
    else: print('MQXA.3L'+ip+'->K1 =  MQXA.3L'+ip+'->K1 +',-B1pQbl,' ;', file=f4q_for_sbssim)    
    if (-dkq1l <= 0): print('MQXA.1L'+ip+'->K1 =  MQXA.1L'+ip+'->K1',-dkq1l,' ;', file=f4q_for_sbssim)
    else:  print('MQXA.1L'+ip+'->K1 =  MQXA.1L'+ip+'->K1 +',-dkq1l,' ;', file=f4q_for_sbssim)

    if(-B1pQar <= 0):  print('MQXB.B2R'+ip+'->K1 =  MQXB.B2R'+ip+'->K1',-B1pQar,' ;', file=f4q_for_sbssim)   
    else: print('MQXB.B2R'+ip+'->K1 =  MQXB.B2R'+ip+'->K1 +',-B1pQar,' ;', file=f4q_for_sbssim)
    if(-B1pQar <= 0): print('MQXB.A2R'+ip+'->K1 =  MQXB.A2R'+ip+'->K1',-B1pQar,' ;', file=f4q_for_sbssim)
    else: print('MQXB.A2R'+ip+'->K1 =  MQXB.A2R'+ip+'->K1 +',-B1pQar,' ;', file=f4q_for_sbssim)
    if (-B1pQbr <= 0): print('MQXA.3R'+ip+'->K1 =  MQXA.3R'+ip+'->K1',-B1pQbr,' ;', file=f4q_for_sbssim)
    else: print('MQXA.3R'+ip+'->K1 =  MQXA.3R'+ip+'->K1 +',-B1pQbr,' ;', file=f4q_for_sbssim)
    if (-dkq1r <= 0): print('MQXA.1R'+ip+'->K1 =  MQXA.1R'+ip+'->K1',-dkq1r,' ;', file=f4q_for_sbssim)
    else: print('MQXA.1R'+ip+'->K1 =  MQXA.1R'+ip+'->K1 +',-dkq1r,' ;', file=f4q_for_sbssim)

    
 


    f4q_for_sbssim.close()
else:
    call(['rm', '-f', apj_dir +'IR'+ip+'_errors_4quads.madx'])
    call(['rm', '-f', apj_dir +'IR'+ip+'_errors_6quads.madx'])
    call(['rm', '-f',apj_dir + '4qLocalCorrectionIP'+ip+'_APJ.dat'])
    call(['rm', '-f',apj_dir + '6qLocalCorrectionIP'+ip+'_APJ.dat'])
    call(['rm', '-f',apj_dir + '4qcorrections_IP'+ip+'.madx'])
    call(['rm', '-f',apj_dir + '6qcorrections_IP'+ip+'.madx'])
    

###########################Two correctors per IR##########################################

if (os.path.exists(in_dir_b1 + 'b1x.dat') and  os.path.exists(in_dir_b1 + 'b1y.dat')  and  os.path.exists(in_dir_b2 + 'b1x.dat')  and os.path.exists(in_dir_b2 + 'b1y.dat')):
    
    b1x_beam1 = open(in_dir_b1 + 'b1x.dat','r')
    b1y_beam1 = open(in_dir_b1 + 'b1y.dat','r')
    b1x_beam2 = open(in_dir_b2 + 'b1x.dat','r')
    b1y_beam2 = open(in_dir_b2 + 'b1y.dat','r')

    B1pQal2, B1pQar2= common_corr(b1x_beam1,b1y_beam1,b1x_beam2,b1y_beam2,0)

    b1x_beam1.close()
    b1y_beam1.close()
    b1x_beam1.close()
    b1x_beam2.close()
    b1y_beam2.close()

    print("two correctors per IR ")
    print("corrq2l" + ip + " =",B1pQal2,";")
    print("corrq2r" + ip + " ="  ,B1pQar2,";")
    print(" ")   

    f2q_for_apjsim=open(apj_dir +'IR'+ip+'_errors_2quads.madx','w')
    print("errq2l" + ip + " =",-B1pQal2,";", file=f2q_for_apjsim)
    print("errq2r" + ip + " ="  ,-B1pQar2,";", file=f2q_for_apjsim)
    f2q_for_apjsim.close()

    knob_file_2q= apj_dir + '2qLocalCorrectionIP'+ip+'_APJ.dat'
    f2q_knob=open(knob_file_2q,'w')
    print('@TITLE' + ' 2qLocalCorrectionIP'+ip+'_APJ', file=f2q_knob)
    print('MQXB2.L'+ip+'/K1', B1pQal2, file=f2q_knob)
    print('MQXB2.R'+ip+'/K1', B1pQar2, file=f2q_knob)
    f2q_knob.close()


    f2q_for_sbssim = open(apj_dir + '2qcorrections_IP'+ip+'.madx','w')
    if(-B1pQal2 <= 0):
        print('MQXB.B2L'+ip+'->K1 =  MQXB.B2L'+ip+'->K1',-B1pQal2,' ;', file=f2q_for_sbssim)
        print('MQXB.A2L'+ip+'->K1 =  MQXB.A2L'+ip+'->K1',-B1pQal2,' ;', file=f2q_for_sbssim)
    else:
        print('MQXB.B2L'+ip+'->K1 =  MQXB.B2L'+ip+'->K1 +',-B1pQal2,' ;', file=f2q_for_sbssim)
        print('MQXB.A2L'+ip+'->K1 =  MQXB.A2L'+ip+'->K1 +',-B1pQal2,' ;', file=f2q_for_sbssim)
    if(-B1pQar2 <= 0):
        print('MQXB.B2R'+ip+'->K1 =  MQXB.B2R'+ip+'->K1',-B1pQar2,' ;', file=f2q_for_sbssim)
        print('MQXB.A2R'+ip+'->K1 =  MQXB.A2R'+ip+'->K1',-B1pQar2,' ;', file=f2q_for_sbssim)
    else:
        print('MQXB.B2R'+ip+'->K1 =  MQXB.B2R'+ip+'->K1 +',-B1pQar2,' ;', file=f2q_for_sbssim)
        print('MQXB.A2R'+ip+'->K1 =  MQXB.A2R'+ip+'->K1 +',-B1pQar2,' ;', file=f2q_for_sbssim)        
        
    f2q_for_sbssim.close()
else:
    print("At least one of the files b1*dat was not found")
    call(['rm', '-f', apj_dir +'IR'+ip+'_errors_2quads.madx'])
    call(['rm', '-f',apj_dir + '2qLocalCorrectionIP'+ip+'_APJ.dat'])
    call(['rm', '-f',apj_dir + '2qcorrections_IP'+ip+'.madx'])
    exit()



##########################Calculation of equivalent gradients from corrections######################
eq_gradf = open(apj_dir + 'eq_grads_from_corrs_IR'+ip+'.dat' ,'w') 
if (os.path.exists(in_dir_b1 + 'b1xl.dat') and os.path.exists(in_dir_b1 + 'b1xr.dat') and os.path.exists(in_dir_b1 + 'b1yl.dat') and os.path.exists(in_dir_b1 + 'b1yr.dat') and          os.path.exists(in_dir_b2 + 'b1xl.dat') and os.path.exists(in_dir_b2 + 'b1xr.dat') and os.path.exists(in_dir_b2 + 'b1yl.dat') and os.path.exists(in_dir_b2 + 'b1yr.dat')):
    
    b1xl_beam1 = open(in_dir_b1 + 'b1xl.dat','r')
    b1xr_beam1 = open(in_dir_b1 + 'b1xr.dat','r')
    b1yl_beam1 = open(in_dir_b1 + 'b1yl.dat','r')
    b1yr_beam1 = open(in_dir_b1 + 'b1yr.dat','r')
    b1xl_beam2 = open(in_dir_b2 + 'b1xl.dat','r')
    b1xr_beam2 = open(in_dir_b2 + 'b1xr.dat','r')
    b1yl_beam2 = open(in_dir_b2 + 'b1yl.dat','r')
    b1yr_beam2 = open(in_dir_b2 + 'b1yr.dat','r')


    #print dkq1l, B1pQal, B1pQbl,dkq1r, B1pQar, B1pQbr
    b1xl_b1, b1yl_b1 = get_equivalent_grads (dkq1l, B1pQal, B1pQbl, b1xl_beam1,b1yl_beam1)
    b1xr_b1, b1yr_b1 = get_equivalent_grads (dkq1r, B1pQar, B1pQbr, b1xr_beam1,b1yr_beam1)
    #b1xl_b2, b1yl_b2 = get_equivalent_grads (dkq1l, B1pQal, B1pQbl,b1xl_beam2,b1yl_beam2)
    #b1xr_b2, b1yr_b2 = get_equivalent_grads (dkq1r, B1pQar, B1pQbr,b1xr_beam2,b1yr_beam2)

    #print "b1xl_b2, b1yl_b2", b1xl_b2, b1yl_b2
    #print "b1xr_b2, b1yr_b2", b1xr_b2, b1yr_b2
    b1xl_beam1.close()
    b1xr_beam1.close()
    b1yl_beam1.close()
    b1yr_beam1.close()
    b1xl_beam2.close()
    b1xr_beam2.close()
    b1yl_beam2.close()
    b1yr_beam2.close()



    b1xl_beam1 = open(in_dir_b1 + 'b1xl.dat','r')
    b1xr_beam1 = open(in_dir_b1 + 'b1xr.dat','r')
    b1yl_beam1 = open(in_dir_b1 + 'b1yl.dat','r')
    b1yr_beam1 = open(in_dir_b1 + 'b1yr.dat','r')
    b1xl_beam2 = open(in_dir_b2 + 'b1xl.dat','r')
    b1xr_beam2 = open(in_dir_b2 + 'b1xr.dat','r')
    b1yl_beam2 = open(in_dir_b2 + 'b1yl.dat','r')
    b1yr_beam2 = open(in_dir_b2 + 'b1yr.dat','r')
    b1x_beam1 = open(in_dir_b1 + 'b1x.dat','r')
    b1y_beam1 = open(in_dir_b1 + 'b1y.dat','r')
    b1x_beam2 = open(in_dir_b2 + 'b1x.dat','r')
    b1y_beam2 = open(in_dir_b2 + 'b1y.dat','r')

    b1x_2q_b1 ,b1y_2q_b1 = get_equivalent_grads_IR(0, B1pQal2, 0, 0, B1pQar2, 0, b1x_beam1,b1y_beam1,b1xl_beam1,b1yl_beam1,b1xr_beam1,b1yr_beam1)
    #b1x_b2 ,b1y_b2 = get_equivalent_grads_IR(0, B1pQal2, 0, 0, B1pQar2, 0, b1x_beam2,b1y_beam2,b1xl_beam2,b1yl_beam2,b1xr_beam2,b1yr_beam2)

    b1xl_beam1.close()
    b1xr_beam1.close()
    b1yl_beam1.close()
    b1yr_beam1.close()
    b1xl_beam2.close()
    b1xr_beam2.close()
    b1yl_beam2.close()
    b1yr_beam2.close()
    b1x_beam1.close()
    b1y_beam1.close()
    b1x_beam1.close()
    b1x_beam2.close()
    b1y_beam2.close()



    b1xl_beam1 = open(in_dir_b1 + 'b1xl.dat','r')
    b1xr_beam1 = open(in_dir_b1 + 'b1xr.dat','r')
    b1yl_beam1 = open(in_dir_b1 + 'b1yl.dat','r')
    b1yr_beam1 = open(in_dir_b1 + 'b1yr.dat','r')
    b1xl_beam2 = open(in_dir_b2 + 'b1xl.dat','r')
    b1xr_beam2 = open(in_dir_b2 + 'b1xr.dat','r')
    b1yl_beam2 = open(in_dir_b2 + 'b1yl.dat','r')
    b1yr_beam2 = open(in_dir_b2 + 'b1yr.dat','r')
    b1x_beam1 = open(in_dir_b1 + 'b1x.dat','r')
    b1y_beam1 = open(in_dir_b1 + 'b1y.dat','r')
    b1x_beam2 = open(in_dir_b2 + 'b1x.dat','r')
    b1y_beam2 = open(in_dir_b2 + 'b1y.dat','r')
    b1x_6q_b1 ,b1y_6q_b1 = get_equivalent_grads_IR(dkq1l, B1pQal, B1pQbl,dkq1r, B1pQar, B1pQbr, b1x_beam1,b1y_beam1,b1xl_beam1,b1yl_beam1,b1xr_beam1,b1yr_beam1)
    #b1x_6q_b2 ,b1y_6q_b2 = get_equivalent_grads_IR(dkq1l, B1pQal, B1pQbl,dkq1r, B1pQar, B1pQbr, b1x_beam2,b1y_beam2,b1xl_beam2,b1yl_beam2,b1xr_beam2,b1yr_beam2)


    #print "b1x_b2, b1y_b2", b1x_b2, b1y_b2
    print("gradient components of equivalent kick estimated from corrections")
    print("b1xl_b1, b1yl_b1", b1xl_b1, b1yl_b1)
    print("b1xr_b1, b1yr_b1", b1xr_b1, b1yr_b1)
    print("b1x_2q_b1, b1y_2q_b1", b1x_2q_b1, b1y_2q_b1)
    print("b1x_4(6)q_b1, b1y_4(6)q_b1", b1x_6q_b1, b1y_6q_b1)
    print(" ")
    
    print("b1xl_b1, b1yl_b1", b1xl_b1, b1yl_b1, file=eq_gradf)
    print("b1xr_b1, b1yr_b1", b1xr_b1, b1yr_b1, file=eq_gradf)
    print("b1x_2q_b1, b1y_2q_b1", b1x_2q_b1, b1y_2q_b1, file=eq_gradf)
    print("b1x_4(6)q_b1, b1y_4(6)q_b1", b1x_6q_b1, b1y_6q_b1, file=eq_gradf)


    
    b1xl_beam1.close()
    b1xr_beam1.close()
    b1yl_beam1.close()
    b1yr_beam1.close()
    b1xl_beam2.close()
    b1xr_beam2.close()
    b1yl_beam2.close()
    b1yr_beam2.close()
    b1x_beam1.close()
    b1y_beam1.close()
    b1x_beam1.close()
    b1x_beam2.close()
    b1y_beam2.close()
eq_gradf.close()    


version_file = open(apj_dir + 'version_get_IR_corrs.dat','a')
print("Version", version, "timestamp", datetime.datetime.now(), file=version_file)
version_file.close()


