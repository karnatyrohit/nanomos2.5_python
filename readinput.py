###########################################################################
##A python function to read the input file using parser #################
###########################################################################

import time
import sys
import numpy as np
from parser import parser

"""# DEVICE DIRECTIVE
Nsd = 2e20; Nbody = 0; Lg_top = 9; Lg_bot = 9; Lsd = 10
overlap_s = -4; overlap_d = -4
dopslope_s=1; dopslope_d=1
t_si = 3; t_top = 1.0; t_bot = 1.0; Te = 300  #tsi, tox_top,tox_bot, temp


# GRID DIRECTIVE
dx=0.3; dy=0.1; refine=1.0


# TRANSPORT DIRECTIVE
transport_model1 = 'clbte'; mu_low = 300; beta = 2; Vel_sat=1e7
ELE_TAUW = 1e-13; ELE_CQ = 1

# BIAS DIRECTIVE
Vg1 = 0.0   #vgtop
Vg2 = 0.0   #vgbot
Vs = 0.0; Vdi = 0.05; Vg_step = 0.1; Vd_step = 0.35
Ng_step = 4; Nd_step = 1; Vd_initial = 0.1

# MATERIAL DIRECTIVE
phi_top = 4.188; phi_bot = 4.188 #phi = wfunc
m_l = 0.98; m_t = 0.19  # mlong and mtran
eps_top = 3.9; eps_bot = 3.9  #kox = eps
bar_top = 3.34; bar_bot = 3.34 # bar = dec
eps_si = 11.7 #ksi

# SOLUTION DIRECTIVE
criterion_outer = 0.001; criterion_inner = 1e-4 # dvmax and dvpois

# OPTIONS DIRECTIVE
t_vall = 3 # valleys= 1 if unprimed  3 if all
max_subband = 1
#True =1, False =0
DG_flag = 1
fermiflag = 1
ox_pnt_flag = 0

# PLOTTING CAPABILITES - yes =1. no =0
plot_IV = 1; plot_Ec3d = 0; plot_Ne3d = 0; plot_Ec_sub = 0; plot_Nesub = 0
plot_Te = 0; plot_Ec_IV = 0; plot_Ne_IV = 0"""

filename = raw_input('Enter filename to be run')
# dirname = input('Enter name output directory')



class Var:
    name = ''
    type = ''
    nval = 0
    val =  0

    """def add_val(self,j):
        for i_val in range(j):
            self.val.append(0)
            print 'ran'"""


class P:
    err = 1
    ncard = ''
    nvar = 0
    var = []
    errmess = ''

    def add_var(self):
        self.var.append(Var())

    """def reset(self):
        err = 1
        ncard = ''
        nvar = 0
        var = []"""

p = P()

fin = open(filename, 'r')

while p.err != -1:
    p.ncard = ''
    p.nvar = 0
    p.var = []
    #p.reset()
    parser(fin, p)

    if p.err == 1 or p.err == 999:
        if p.ncard == 'device':
            for i in np.arange(0, p.nvar):
                if p.var[i].name == 'nsd':
                    if p.var[i].type == 'number':
                        Nsd = p.var[i].val
                    else:
                        print 'Invalid assignment for nsd'
                        exit()

                elif p.var[i].name == 'nbody':
                    if p.var[i].type == 'number':
                        Nbody = p.var[i].val
                    else:
                        print 'Invalid assignment for nbody'
                        exit()

                elif p.var[i].name == 'lgtop':
                    if p.var[i].type == 'number':
                        Lg_top = p.var[i].val
                    else:
                        print 'Invalid assignment for lgtop'
                        exit()

                elif p.var[i].name == 'lgbot':
                    if p.var[i].type == 'number':
                        Lg_bot = p.var[i].val
                    else:
                        print 'Invalid assignment for lgbot'
                        exit()

                elif p.var[i].name == 'lsd':
                    if p.var[i].type == 'number':
                        Lsd=p.var[i].val
                    else:
                        print 'Invalid assignment for lsd'
                        exit()

                elif p.var[i].name == 'overlap_s':
                    if p.var[i].type == 'number':
                        overlap_s=p.var[i].val
                    else:
                        print 'Invalid assignment for overlap_s'
                        exit()

                elif p.var[i].name == 'overlap_d':
                    if p.var[i].type == 'number':
                        overlap_d=p.var[i].val
                    else:
                        print 'Invalid assignment for overlap_d'
                        exit()

                elif p.var[i].name == 'dopslope_s':
                    if p.var[i].type == 'number':
                        dopslope_s=p.var[i].val
                    else:
                        print 'Invalid assignment for dopslope_s'
                        exit()

                elif p.var[i].name == 'dopslope_d':
                    if p.var[i].type == 'number':
                        dopslope_d=p.var[i].val
                    else:
                        print 'Invalid assignment for dopslope_d'
                        exit()

                elif p.var[i].name == 'tsi':
                    if p.var[i].type == 'number':
                        t_si=p.var[i].val
                    else:
                        print 'Invalid assignment for tsi'
                        exit()

                elif p.var[i].name == 'tox_top':
                    if p.var[i].type == 'number':
                        t_top=p.var[i].val
                    else:
                        print 'Invalid assignment for tox_top'
                        exit()

                elif p.var[i].name == 'tox_bot':
                    if p.var[i].type == 'number':
                        t_bot=p.var[i].val
                    else:
                        print 'Invalid assignment for tox_bot'
                        exit()

                elif p.var[i].name == 'temp':
                    if p.var[i].type == 'number':
                        Te=p.var[i].val
                    else:
                        print 'Invalid assignment for temp'
                        exit()

        elif p.ncard == 'grid':
            for i in np.arange(0,p.nvar):
                if p.var[i].name == 'dx':
                    if p.var[i].type == 'number':
                        dx=p.var[i].val
                    else:
                        print 'Invalid assignment for dx'
                        exit()

                elif p.var[i].name == 'dy':
                    if p.var[i].type == 'number':
                        dy=p.var[i].val
                    else:
                        print 'Invalid assignment for dy'
                        exit()

                elif p.var[i].name == 'refine':
                    if p.var[i].type == 'number':
                        refine=p.var[i].val
                    else:
                        print 'Invalid assignment for refine'
                        exit()

        elif p.ncard == 'transport':
            for i in np.arange(0,p.nvar):
                if p.var[i].name == 'model':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'dd':
                            transport_model1 = 2
                        elif p.var[i].val =='clbte':
                            transport_model1 = 3
                        elif p.var[i].val == 'qbte':
                            transport_model1 = 4
                        elif p.var[i].val == 'qdte':
                            transport_model1 = 5
                        elif p.var[i].val =='et':
                            transport_model1 = 6
                        else:
                            print '****** ERROR !!! MODEL CAN ONLY BE DD/CLBTE/QBTE/ET/QDTE *******'
                            exit()

                    else:
                        print 'Invalid assignment for model'
                        exit()

                elif p.var[i].name == 'mu_low':
                    if p.var[i].type == 'number':
                        mu_low = p.var[i].val
                    else:
                        print 'Invalid assignment for mu_low'
                        exit()

                elif p.var[i].name == 'beta':
                    if p.var[i].type == 'number':
                        beta=p.var[i].val
                    else:
                        print 'Invalid assignment for beta'
                        exit()

                elif p.var[i].name == 'vsat':
                    if p.var[i].type == 'number':
                        Vel_sat=p.var[i].val
                    else:
                        print 'Invalid assignment for vsat'
                        exit()

                elif p.var[i].name == 'ELE_TAUW':
                    if p.var[i].type == 'number':
                        ELE_TAUW=p.var[i].val
                    else:
                        print 'Invalid assignment for ELE_TAUW'
                        exit()

                elif p.var[i].name == 'ELE_CQ':
                    if p.var[i].type == 'number':
                        ELE_CQ=p.var[i].val
                    else:
                        print 'Invalid assignment for ELE_CQ'
                        exit()

        elif p.ncard =='bias':
            for i in np.arange(0,p.nvar):
                if p.var[i].name == 'vgtop':
                    if p.var[i].type == 'number':
                        Vg1=p.var[i].val
                    else:
                        print 'Invalid assignment for vgtop'
                        exit()

                elif p.var[i].name == 'vgbot':
                    if p.var[i].type == 'number':
                        Vg2=p.var[i].val
                    else:
                        print 'Invalid assignment for vgbot'
                        exit()

                elif p.var[i].name == 'vs':
                    if p.var[i].type == 'number':
                        Vs=p.var[i].val
                    else:
                        print 'Invalid assignment for vs'
                        exit()

                elif p.var[i].name == 'vd':
                    if p.var[i].type == 'number':
                        Vdi=p.var[i].val
                    else:
                        print 'Invalid assignment for vd'
                        exit()

                elif p.var[i].name == 'vgstep':
                    if p.var[i].type == 'number':
                        Vg_step=p.var[i].val
                    else:
                        print 'Invalid assignment for vgstep'
                        exit()

                elif p.var[i].name == 'vdstep':
                    if p.var[i].type == 'number':
                        Vd_step=p.var[i].val
                    else:
                        print 'Invalid assignment for vdstep'
                        exit()

                elif p.var[i].name == 'ngstep':
                    if p.var[i].type == 'number':
                        Ng_step=p.var[i].val
                    else:
                        print 'Invalid assignment for ngstep'
                        exit()

                elif p.var[i].name == 'ndstep':
                    if p.var[i].type == 'number':
                        Nd_step=p.var[i].val
                    else:
                        print 'Invalid assignment for ndstep'
                        exit()

                elif p.var[i].name == 'vd_initial':
                    if p.var[i].type == 'number':
                        Vd_initial=p.var[i].val
                    else:
                        print 'Invalid assignment for vd_initial'
                        exit()




        elif p.ncard == 'material':
            for i in np.arange(0,p.nvar):
                if p.var[i].name == 'wfunc_top':
                    if p.var[i].type == 'number':
                        phi_top =p.var[i].val
                    else:
                        print 'Invalid assignment for wfunc_top'
                        exit()

                elif p.var[i].name == 'wfunc_bot':
                    if p.var[i].type == 'number':
                        phi_bot=p.var[i].val
                    else:
                        print 'Invalid assignment for wfunc_bot'
                        exit()

                elif p.var[i].name == 'mlong':
                    if p.var[i].type == 'number':
                        m_l=p.var[i].val
                    else:
                        print 'Invalid assignment for mlong'
                        exit()

                elif p.var[i].name == 'mtran':
                    if p.var[i].type == 'number':
                        m_t=p.var[i].val
                    else:
                        print 'Invalid assignment for mtran'
                        exit()

                elif p.var[i].name == 'kox_top':
                    if p.var[i].type == 'number':
                        eps_top=p.var[i].val
                    else:
                        print 'Invalid assignment for kox_top'
                        exit()

                elif p.var[i].name == 'kox_bot':
                    if p.var[i].type == 'number':
                        eps_bot=p.var[i].val
                    else:
                        print 'Invalid assignment for kox_bot'
                        exit()

                elif p.var[i].name == 'dec_top':
                    if p.var[i].type == 'number':
                        bar_top=p.var[i].val
                    else:
                        print 'Invalid assignment for dec_top'
                        exit()

                elif p.var[i].name == 'dec_bot':
                    if p.var[i].type == 'number':
                        bar_bot=p.var[i].val
                    else:
                        print 'Invalid assignment for dec_bot'
                        exit()

                elif p.var[i].name == 'ksi':
                    if p.var[i].type == 'number':
                        eps_si=p.var[i].val
                    else:
                        print 'Invalid assignment for ksi'
                        exit()




        elif p.ncard == 'solve':
            for i in np.arange(0,p.nvar):
                if p.var[i].name == 'dvmax':
                    if p.var[i].type == 'number':
                        criterion_outer=p.var[i].val
                    else:
                        print 'Invalid assignment for dvmax'
                        exit()


                if p.var[i].name == 'dvpois':
                    if p.var[i].type == 'number':
                        criterion_inner=p.var[i].val
                    else:
                        print 'Invalid assignment for dvmax'
                        exit()



#*****************************************************************************
        elif p.ncard == 'plots':
            for i in np.arange(0,p.nvar):
                #Iv plot
                if p.var[i].name == 'I_V':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'y':
                            plot_IV=1
                        else:
                            plot_IV=0

                    else:
                        print 'Invalid assignment for I_V'
                        exit()


                #Ec plot
                if p.var[i].name == 'Ec3d':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'y':
                            plot_Ec3d=1
                        else:
                            plot_Ec3d=0

                    else:
                        print 'Invalid assignment for Ec'
                        exit()


                    #Ne plot
                if p.var[i].name == 'Ne3d':
                    if p.var[i].type == 'string':
                        if p.var[i].val =='y':
                            plot_Ne3d=1
                        else:
                            plot_Ne3d=0

                    else:
                        print 'Invalid assignment for Ne'
                        exit()


                    #Ec_sub plot
                if p.var[i].name == 'Ec_sub':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'y':
                            plot_Ecsub=1
                        else:
                            plot_Ecsub=0

                    else:
                        print 'Invalid assignment for Ecsub'
                        exit()


                    #Ne_sub plot
                if p.var[i].name == 'Ne_sub':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'y':
                            plot_Nesub=1
                        else:
                            plot_Nesub=0

                    else:
                        print 'Invalid assignment for Ne_sub'
                        exit()


                    #Te plot
                if p.var[i].name == 'Te':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'y':
                            plot_Te=1
                        else:
                            plot_Te=0

                    else:
                        print 'Invalid assignment for Te'
                        exit()


                    #Ec_IV plot
                if p.var[i].name == 'Ec_IV':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'y':
                            plot_Ec_IV=1
                        else:
                            plot_Ec_IV=0

                    else:
                        print 'Invalid assignment for Ec_IV'
                        exit()


                    #Ne_IV plot
                if p.var[i].name == 'Ne_IV':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'y':
                            plot_Ne_IV=1
                        else:
                            plot_Ne_IV=0

                    else:
                        print 'Invalid assignment for Ne_IV'
                        exit()



#*****************************************************************************
        elif p.ncard == 'options':
            for i in np.arange(0,p.nvar):
                if p.var[i].name == 'valleys':
                    if p.var[i].type == 'string':
                        if p.var[i].val == 'unprimed':
                            t_vall=1
                        else:
                            t_vall=3

                    else:
                        print 'Invalid assignment for valleys'
                        exit()

                elif p.var[i].name == 'num_subbands':
                    if p.var[i].type == 'number':
                        max_subband=p.var[i].val
                    else:
                        print 'Invalid assignment for num_subbands'
                        exit()

                elif p.var[i].name == 'dg':
                    if p.var[i].type == 'string':
                        if p.var[i].val.lower() == 'true':
                            DG_flag = 1
                        elif p.var[i].val.lower() == 'false':
                            DG_flag = 0
                        else:
                            print 'Invalid assignment for dg'
                            exit()
                    else:
                        print 'Invalid assignment for dg'
                        exit()

                elif p.var[i].name == 'fermi':
                    if p.var[i].type == 'string':
                        if p.var[i].val.lower() == 'true':
                            fermiflag = 1
                        elif p.var[i].val.lower() == 'false':
                            fermiflag = 0
                        else:
                            print 'Invalid assignment for fermi'
                            exit()
                    else:
                        print 'Invalid assignment for fermi'
                        exit()

                elif p.var[i].name == 'ox_penetrate':
                    if p.var[i].type == 'string':
                        if p.var[i].val.lower() == 'true':
                            ox_pnt_flag = 1
                        elif p.var[i].val.lower() == 'false':
                            ox_pnt_flag = 0
                        else:
                            print 'Invalid assignment for ox_penetrate'
                            exit()
                    else:
                        print 'Invalid assignment for ox_penetrate'
                        exit()

fin.close()
p.err = 1

if transport_model1 != 3:
        if (Vg1>1 or Vg2>1 or Vs>1 or Vdi>1):
            print 'Too High of Voltage'
            exit()

#if (Nd_step>0&Ng_step>0)
#    dump='Only one voltage sweep can be defined'
#    return


if Ng_step > 20 or Nd_step > 20:
    print 'More than 20 bias points specified'
    exit()



#physics constants
eps_o = 8.85e-12
psi_si = 4.05
q = 1.6e-19
k_B = 1.38e-23
h_bar = 1.05e-34
m_e = 0.91e-30
Nc = 2.8e25
Ncc = 2*m_e*m_t*k_B*Te/(np.pi*(h_bar**2))

###############modifications##########################################
Nsd *= 1e6
Nbody *= 1e6
Lg_top *= 1e-9
Lg_bot *= 1e-9
Lsd *= 1e-9
overlap_s *= 1e-9
overlap_d *= 1e-9
dopslope_s *= 1e-9
dopslope_d *= 1e-9
t_si *= 1e-9
t_bot *= 1e-9
t_top *= 1e-9
dx *= 1e-9
dy *= 1e-9
mu_low *= 1e-4
Vel_sat *= 1e-2


Lsd = round(Lsd/dx)*dx
Lg_top = round(Lg_top/dx)*dx
Lg_bot = round(Lg_bot/dx)*dx
t_top = round(t_top/dy)*dy
t_bot = round(t_bot/dy)*dy
t_si = round(t_si/dy)*dy


class Nsd1:
    value = Nsd


class Nbody1:
    value = Nbody


"""if transport_model1 == 'dd':
    transport_model1 = 2
elif transport_model1 == 'clbte':
    transport_model1 = 3
elif transport_model1 == 'qbte':
    transport_model1 = 4
elif transport_model1 == 'ddte':
    transport_model1 = 5
elif transport_model1 == 'et':
    transport_model1 = 6
else:
    print '****** ERROR !!! MODEL CAN ONLY BE DD/CLBTE/QBTE/ET/QDTE *******'
    exit()"""
######################################################################


class transportmodel:
    value = transport_model1

#####################################################################
# CHECKING ALL INPUT VARIABLES
#####################################################################
if DG_flag != 1 and DG_flag != 0:
    sys.exit('******ERROR, DG_flag can only be 0 or 1!!!******')

if fermiflag != 1 and fermiflag != 0:
    sys.exit('******ERROR, fermi_flag can only be 0 or 1!!!******')

if ox_pnt_flag != 1 and ox_pnt_flag != 0:
    sys.exit('******ERROR, ox_pnt_flag can only be 0 or 1!!!******')

if Lsd > 15e-9:
    Lsd = 15e-9
    ver = '******NOTE! LSD is reduced to 15nm !!!******'
    print ver
    time.sleep(2)

if Lg_top > 50e-9:
    Lg_top = 50e-9
    ver = '******NOTE! LGTOP is reduced to 50nm !!!******'
    print ver
    time.sleep(2)

if Lg_bot > 1.2*Lg_top:
    Lg_bot = 1.2*Lg_top
    ver = '******NOTE! LGBOT is reduced to 1.2*LGTOP !!!******'
    print ver
    time.sleep(2)

if t_si > 5e-9:
    ver = '******Please specify a TSI less than 5nm !!!******'
    print ver
    time.sleep(2)

if t_vall != 1 and t_vall != 3:
    sys.exit('******ERROR, VALLEYS can only be UNPRIMED or ALL!!!******')

if max_subband > 3:
    max_subband = 3
    ver = '******Note! NUM-SUBBAND is limited to 3 for each VALLEY!!!******'
    print ver
    time.sleep(2)


# if t_si<=2e-9
#  t_vall=1
#  max_subband=1
#  ver='******Note! For TSI<2nm, only one subband is assumed!!!******'
#  print ver
#  time.sleep(2)
#

if criterion_outer <= 1e-4 and (transport_model1 == 2 or transport_model1 == 3):
    criterion_outer = 1e-4
    ver = '******Note! DVMAX for models CLBTE and QBTE is limited to 0.1meV!!!******'
    print ver
    time.sleep(2)


if t_si < 2e-9 and Te > 250:
    dummy_flag = 0
else:
    dummy_flag = 0.5

if transport_model1 == 6:
    fermiflag = 0

class fermiflag1:
    value = fermiflag

class Vdc:
    value = Vdi
