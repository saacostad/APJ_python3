

import sys
import os
import math
from scipy.integrate import quad
from subprocess import call
import argparse
import datetime

version=4
#Version 4 to be able to run the software in any system diferent from cern afs. Needs aditional files. March 15 2024.
#Version 3 to ln to acc-models-lhc and to run madx nominal_file.madx from the directory  where ln to acc-models-lhc was created. March 4 2024
#Version 2 to forward madx output to device or file different to the screen. March 1 2024
#Version 1 to redifine outpu_dir, default for madx program and location of apj_files. Feb 28 2024
#Version 0 was copied  from /Users/jfcp/bin/gen_nom_files3 on February 13/2024.
print("get_nom_files Version",version)

def beta_integral(L, alpha0, beta0,  KK, plane, beam):
    """calc beta av in foc magnet
    b ... beta at waist
    L ...  L*
    KK ... Quadrupole Gradient
    K ... Sqrt of Quadrupole Gradient
    KL ... Sqrt of Quadrupole Gradient times Quadrupole length
    """
    K = math.sqrt(abs(KK))
    KL = K*L
    b = beta0 / (1 + alpha0**2)
    if (((KK > 0) and (plane == '-x') and (beam == '1')) or ((KK < 0) and (plane == '-y') and (beam == '1')) or ((KK > 0) and (plane == '-y') and (beam == '2')) or ((KK < 0) and (plane == '-x') and (beam == '2'))):
        sin2KL = ((math.sin(2 * KL)) / (2 * KL))
        beta_av = 0.5 * beta0 * (1 + sin2KL) + alpha0 * ((math.sin(KL) ** 2) / (KL * K)) + (1 - sin2KL) / (2 * b * K ** 2)
    else:
        sinh2KL = ((math.sinh(2 * KL)) / (2 * KL))
        beta_av = 0.5 * beta0 * (1 + sinh2KL) + alpha0 * ((math.sinh(KL) ** 2) / (KL * K)) + (sinh2KL - 1) / (2 * b * K ** 2)
    beta_int = beta_av*L
    return beta_int

def beta_f(x, alpha0, beta0,  KK, plane, beam):
    """ 
    KK ... Quadrupole Gradient
    K ... Sqrt of Quadrupole Gradient
    """
    K = math.sqrt(abs(KK))
    gamma =  (1 + alpha0**2)/beta0
    if (((KK > 0) and (plane == '-x') and (beam == '1')) or ((KK < 0) and (plane == '-y') and (beam == '1')) or ((KK > 0) and (plane == '-y') and (beam == '2')) or ((KK < 0) and (plane == '-x') and (beam == '2'))):
        SS = (math.sin(K*x))/K
        CC = math.cos(K*x)
    else:
        SS = (math.sinh(K*x))/K
        CC = math.cosh(K*x)        
    return beta0*CC**2 +2*CC*SS*alpha0 + gamma*SS**2


def beta_sk(x,alphax, betax,alphay, betay):
    return math.sqrt((-2*alphax*x + betax)*(-2*alphay*x + betay))
    

    
def beta_integral2(L, alpha0, beta0,  KK, plane, beam):
    def beta_aux(x):
        return beta_f(x, alpha0, beta0,  KK, plane, beam)
    res, err = quad(beta_aux, 0, L)
    return res
    

    


def beta_int_skew(L, alphax, betax,alphay, betay):
    def beta_sk_aux(x):
        return beta_sk(x,alphax, betax,alphay, betay)
    res, err = quad(beta_sk_aux, 0, L)
    return res        

def betaxbetay_int(L, alphax, betax,alphay, betay, KK, beam):
    def betxbety_aux(x):
        return math.sqrt((beta_f(x, alphax, betax,  KK, '-x', beam))*(beta_f(x, alphay, betay,  KK, '-y', beam)))
    res, err = quad(betxbety_aux, 0, L)
    return res           
    




def rename_mag(magnet_name):
    mag = magnet_name.split('.')
    if ((mag[0] =='MQXA' and mag[1][0] =='3') or (mag[0] =='MQXFA' and mag[1][1] =='3')): newname = 'Q3'+ mag[1][-2]+mag[1][-1]
    if ((mag[0] =='MQXA' and mag[1][0] =='1') or (mag[0] =='MQXFA' and mag[1][1] =='1')): newname = 'Q1'+ mag[1][-2]+mag[1][-1]
    if (mag[0] =='MQXB' or mag[0] =='MQXFB' ):newname = 'Q2'+ mag[1][-2]+mag[1][-1]
    if (mag[0] =='MQSX'): newname = 'MQSX'+'.'+mag[1]
    if (mag[0] == 'MQY' and mag[1][0] =='4'):newname = 'Q4' + mag[1][-2]+mag[1][-1]
    if (mag[0] == 'MQY' and mag[1][0] =='5'): newname = 'Q5' + mag[1][-2]+mag[1][-1]
    if (mag[0] == 'MQML' and mag[1][0] =='5'): newname = 'Q5' + mag[1][-2]+mag[1][-1]
    if (mag[0] == 'MQML' and mag[1][0] =='6'): newname = 'Q6' + mag[1][-2]+mag[1][-1]   
    return newname





#from subprocess import run



#if (len(sys.argv) == 3):
#    twissfilename=sys.argv[1]
#    latticefilename=sys.argv[2]
#else:
#	if (len(sys.argv) == 2):
#		twissfilename=sys.argv[1]
#		latticefilename="lattice.asc"
#	else:
#		twissfilename="my_model"
#		latticefilename="lattice.asc"

#print "Leyendo archivo",twissfilename
parser = argparse.ArgumentParser()
parser.add_argument(
        "-b", "--beam",
        help="beam can be 1 or 2",
        choices=['1', '2'],
        required=True,
        dest="beam"
    )
parser.add_argument(
        "-id", "--input_dir",
        help="Directory where the original job.twiss.madx is",
        required=True,
        dest="input_dir"
    )
parser.add_argument(
        "-od", "--out_dir",
        help="Directory where all output files will be written",
        required=True,
        dest="out_dir"
    )
parser.add_argument(
        "-mx", "--mad_program",
        help="Location of madx executable",
        dest="mad_program",
        #default='/Users/jfcp/bin/madx507'
        default='/afs/cern.ch/user/m/mad/madx/releases/5.07.00/madx-linux64-gnu'
    )
parser.add_argument(
        "-a", "--accel",
        help="accelerator name",
        choices=['lhc','hl_lhc'],
        dest="accel",
        default='lhc'
    )
parser.add_argument(
   "-nc", "--no_cern_afs",
    help="To be able to run job.twiss.madx from any system different to cern afs. Macros and modifiers are called from apj/lhc/madx_macros/ and apj/lhc/madx_modifiers/",
    action="store_true",
    dest="no_cern_afs"
)

args = parser.parse_args()

#if (args.accel == 'lhc'): apj_files_dir='/Users/jfcp/apj/lhc/'
#if (args.accel == 'hl_lhc'):apj_files_dir='/Users/jfcp/apj/hl_lhc/' 

if (args.accel == 'lhc'): apj_files_dir = os.path.dirname(sys.argv[0]) + '/lhc/'
if (args.accel == 'hl_lhc'):apj_files_dir= os.path.dirname(sys.argv[0]) + '/hl_lhc/'    


original_model_dir=args.input_dir
beam = args.beam
out_dir=args.out_dir
#print out_dir

call(["mkdir","-p", out_dir])
job_orig_file=original_model_dir + "/" + "job.twiss.madx"
acc_models_lhc=out_dir + 'acc-models-lhc'
#call(["ln", "-s", "acc-models-lhc", acc_models_lhc])
call(["ln", "-s", original_model_dir + "acc-models-lhc", acc_models_lhc])
#call(["ln", "-s", "/Users/jfcp/Beta-Beat.src/acc-models-lhc", acc-models-lhc])
twiss_orig=original_model_dir + "/" + "twiss.dat"
#call(["cp",twiss_orig,out_dir+'/'+ 'twiss_orig.dat'])
call(["cp",job_orig_file,out_dir])
job_file =  out_dir + 'job.twiss.madx'
template_file=apj_files_dir + 'b'+ beam +'/' + 'nom.madx'
nominal_file= out_dir  + "nominal_lattice.madx"


twiss_file = '"'+out_dir  + "twiss.dat" + '"'
twiss_command='exec, do_twiss_monitors(LHCB'+ beam + ',' +  twiss_file + ', 0.0);'
jobf = open(job_file,'r')
templ=open(template_file,'r')
nomf=open(nominal_file,'w')
for line in jobf:
    if (args.no_cern_afs):
            if ('madx_macros' in line):
                lines = line.split('"')
                macro_path = lines[1]
                new_path = apj_files_dir +  'madx_macros/' + os.path.basename(macro_path)
                print('call , file = "' + new_path + '";', file=nomf)
            elif ('madx_modifiers' in line):
                lines = line.split('"')
                macro_path = lines[1]
                new_path = apj_files_dir +  'madx_modifiers/' + os.path.basename(macro_path)
                print('call , file = "' + new_path + '";', file=nomf)    
            elif ('twiss_ac' in line or 'twiss_adt' in line or 'twiss_elements' in line or 'twiss_monitors' in line):
                if ('twiss_monitors' in line):
                    print(twiss_command, file=nomf)
                else:
                    print('!', end=' ', file=nomf)
                    print(line.strip(), file=nomf)       
            else: print(line.strip(), file=nomf)
    else:
            if ('twiss_ac' in line or 'twiss_adt' in line or 'twiss_elements' in line or 'twiss_monitors' in line):
                if ('twiss_monitors' in line):
                    print(twiss_command, file=nomf)
                else:
                    print('!', end=' ', file=nomf)
                    print(line.strip(), file=nomf)       
            else: print(line.strip(), file=nomf)        

jobf.close()
nomf.close()
call(['cp',nominal_file,out_dir + 'job.twiss_out.madx'])
nomf=open(nominal_file,'a')
for line in templ:
    if ('twiss.optics' in line):
        new_file = out_dir + 'twiss.optics'
        newline=line.replace("twiss.optics",new_file)
    elif ('twiss_c.optics' in line):
        new_file = out_dir + 'twiss_c.optics'
        newline=line.replace("twiss_c.optics",new_file)
    elif ('twiss_shifted.dat' in line):
        new_file = out_dir  + 'twiss_shifted.dat'
        newline=line.replace("twiss_shifted.dat",new_file)
    elif ('Quad_KyL.txt' in line):
        new_file = out_dir  + 'Quad_KyL.txt'
        newline=line.replace("Quad_KyL.txt",new_file)
    elif ('my_model' in line):
        new_file = out_dir  + 'my_model'
        newline=line.replace("my_model",new_file)
    else:
        newline = line
    print(newline.strip(), file=nomf)
templ.close()    
nomf.close()

madx_exec=args.mad_program

#call([madx_exec, nominal_file])
print(nominal_file, " was created, running madx...")
current_dir = os.getcwd()
os.chdir(out_dir)
with open(out_dir + "madx.out", "wb") as madxout:  call([madx_exec, nominal_file], stdout=madxout)
os.chdir(current_dir)


twissfilename= out_dir + "my_model"
latticefilename=out_dir + "lattice.asc"
twissfile = open(twissfilename,'r')
latticefile = open(latticefilename,'w')
for line in twissfile:
    sline = line.split(None)
    if not(('@' in sline) or ('$' in sline) or ('#' in sline)):#header
        #print sline
    #else:
        if ('*' in sline):#column names
            colname = sline.index('NAME')-1
            cols = sline.index('S')-1
            colbetx = sline.index('BETX')-1
            colbety = sline.index('BETY')-1
            colmux = sline.index('MUX')-1
            colmuy = sline.index('MUY')-1
            colalfx = sline.index('ALFX')-1
            colalfy = sline.index('ALFY')-1             
        else:#print output
            print(sline[colname], sline[cols],sline[colbetx],'0',sline[colbety],sline[colalfx],sline[colalfy],'0','0','0',sline[colmux],sline[colmuy], file=latticefile)
twissfile.close()
latticefile.close()
print("Archivo",latticefilename," was created")



twiss_optics=out_dir  + "twiss.optics" #twiss.optics
QLyKf_file=out_dir  + "Quad_KyL.txt"
integral_out=out_dir  + "integrals.dat"

QLyKf=open(QLyKf_file,'r')
twiss=open(twiss_optics,'r')
fout = open(integral_out,'w')

name=[]
length=[]
strength=[]

for i in QLyKf:
    il = i.split(None)
    if (len(il) > 0):
        if ('MQ' in il[0]):
            name.append(il[0])
            length.append(float(il[3].strip(',')))
            strength.append(float(il[4]))





namet=[]
ss=[]
betx=[]
mux=[]
bety=[]
muy=[]
alfx=[]
alfy=[]
flag=0

for j in twiss:
    jl=j.split(None)
    if (flag):
        namet.append(jl[0].strip('"'))
        ss.append(float(jl[1]))
        betx.append(float(jl[2]))
        alfx.append(float(jl[8]))
        bety.append(float(jl[4]))
        alfy.append(float(jl[9]))
        mux.append(float(jl[3]))
        muy.append(float(jl[5]))
    if ('$' in jl[0]): flag = 1

#print namet

for k in name:

    if ('MQSX' in k ):
        beta_intx = beta_int_skew(length[name.index(k)], alfx[namet.index(k)], betx[namet.index(k)], alfy[namet.index(k)], bety[namet.index(k)])
        beta_inty = beta_intx
        betxbety = beta_intx
    else:
        #print k, strength[name.index(k)]
        if (strength[name.index(k)] == 0):
            beta_intx = 0
            beta_inty = 0
            betxbety = 0
        else:            
            beta_intx = beta_integral2(length[name.index(k)], alfx[namet.index(k)], betx[namet.index(k)], strength[name.index(k)], '-x', beam)
            beta_inty = beta_integral2(length[name.index(k)], alfy[namet.index(k)], bety[namet.index(k)], strength[name.index(k)], '-y', beam)
            betxbety =  betaxbetay_int(length[name.index(k)], alfx[namet.index(k)], betx[namet.index(k)], alfy[namet.index(k)], bety[namet.index(k)],strength[name.index(k)], beam)
    k2 = '"'+k+'"'
    print(k2,  beta_intx, beta_inty, rename_mag(k), mux[namet.index(k)], muy[namet.index(k)], betxbety, file=fout)



QLyKf.close()
twiss.close()
fout.close()
print(integral_out, ' was created')

twiss_optics_c=out_dir  + "twiss_c.optics"
twiss=open(twiss_optics,'r')
twiss_c=open(twiss_optics_c,'r')

optics_file= out_dir  + 'optics.out'
opt=open(optics_file,'w')
for line in twiss:
    lines=line.split()
    line_c = twiss_c.readline()
    line_cs=line_c.split()
    if ('MB' in line):
        s_out = float(lines[1])
        s_c = float(line_cs[1])
        s_in = s_c - (s_out -s_c)
        print(lines[0], s_in, 0, file=opt)
        print(lines[0], s_in, 1, file=opt)
        print(lines[0], s_out, 1, file=opt)
        print(lines[0], s_out, 0, file=opt)

    if ('MQ' in line):
        s_out = float(lines[1])
        s_c = float(line_cs[1])
        s_in = s_c - (s_out -s_c)
        print(lines[0], s_in, 0, file=opt)
        print(lines[0], s_in, 2, file=opt)
        print(lines[0], s_out, 2, file=opt)
        print(lines[0], s_out, 0, file=opt)

    if ('BPM' in line):
        s_out = float(lines[1])
        s_c = float(line_cs[1])
        s_in = s_c - (s_out -s_c)
        print(lines[0], s_in, 0, file=opt)
        print(lines[0], s_in, 0.25, file=opt)
        print(lines[0], s_out, 0.25, file=opt)
        print(lines[0], s_out, 0, file=opt)

twiss.close()
twiss_c.close()

print(optics_file, ' was created \n')


version_file = open(out_dir + 'version_get_nom_files.dat','a')
print("Version", version, "timestamp", datetime.datetime.now(), file=version_file)
version_file.close()
