###########################################################################
##A python function to read the input file using parser #################
###########################################################################

import time
import sys
import numpy as np

# DEVICE DIRECTIVE
Nsd = 2e20; Nbody = 0; Lg_top = 9; Lg_bot = 9; Lsd = 10
overlap_s = -4; overlap_d = -4
dopslope_s=1; dopslope_d=1
t_si = 3; t_top = 1.0; t_bot = 1.0; Te = 300  #tsi, tox_top,tox_bot, temp


#

# GRID DIRECTIVE
dx=0.3; dy=0.1; refine=1


# TRANSPORT DIRECTIVE
transport_model = 'clbte'; mu_low = 300; beta = 2; Vel_sat=1e7
ELE_TAUW = 1e-13; ELE_CQ = 1

# BIAS DIRECTIVE
Vg1 = 0.0   #vgtop
Vg2 = 0.0   #vgbot
Vs = 0.0; Vd = 0.05; Vg_step = 0.1; Vd_step = 0.35
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
fermi_flag = 1
ox_pnt_flag = 0

# PLOTTING CAPABILITES - yes =1. no =0
plot_IV = 1; plot_Ec3d = 0; plot_Ne3d = 0; plot_Ec_sub = 0; plot_Nesub = 0
plot_Te = 0; plot_Ec_IV = 0; plot_Ne_IV = 0

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
dopslope_d *= 1e-9
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


if transport_model == 'dd':
    transport_model = 2
elif transport_model == 'clbte':
    transport_model = 3
elif transport_model == 'qbte':
    transport_model = 4
elif transport_model == 'ddte':
    transport_model = 5
elif transport_model == 'et':
    transport_model = 6
else:
    print '****** ERROR !!! MODEL CAN ONLY BE DD/CLBTE/QBTE/ET/QDTE *******'
######################################################################

#####################################################################
# CHECKING ALL INPUT VARIABLES
#####################################################################
if DG_flag != 1 and DG_flag != 0:
    sys.exit('******ERROR, DG_flag can only be 0 or 1!!!******')

if fermi_flag != 1 and fermi_flag != 0:
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

if criterion_outer <= 1e-4 and (transport_model == 2 or transport_model == 3):
    criterion_outer = 1e-4
    ver = '******Note! DVMAX for models CLBTE and QBTE is limited to 0.1meV!!!******'
    print ver
    time.sleep(2)


if t_si < 2e-9 and Te > 250:
    dummy_flag = 0
else:
    dummy_flag = 1/2


if transport_model == 6:
    fermi_flag = 0
