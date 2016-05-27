# A M file to specify the function that generates the doping profile
# With different overlaps and different doping gradient on the drain and source side
# By sebastien Goasguen October 2001/ Purdue

from readinput import *
import math

def doping(Nx, Ny, Ntotal, junction_l, junction_r, Nd):
    Lsda = round(Lsd/dx)
    Lg_topa = round(Lg_top/dx)
    Lg_bota = round(Lg_bot/dx)
    t_topa = round(t_top/dy)
    t_bota = round(t_bot/dy)
    t_sia = round(t_si/dy)

    if (dopslope_s!=0 or dopslope_d!=0):
        decay_lens = (2/math.sqrt(2.3)*dopslope_s)
        decay_lend = (2/math.sqrt(2.3)*dopslope_d)

    ########################ASSUMING GAUSSIAN DISTRIBUTION########################
    for iii_row in range(((Nx*(t_topa+1))/Nx),((Ntotal-Nx*t_bota)/Nx-2)):
        for iii_col in range(1,junction_l):
            i_node = iii_row*Nx+iii_col
            Nd[i_node] = N_sd-N_body


      for iii_col=junction_l+1:junction_r-1
        i_node=iii_row*Nx+iii_col
        if (dopslope_s~=0 & dopslope_d~=0)
            Nd(i_node)=N_sd*exp(-((iii_col-junction_l)*dx/decay_lens)^2)+...
                   N_sd*exp(-((iii_col-junction_r)*dx/decay_l)^2)-...
                   N_body
        elseif (dopslope_s~=0 & dopslope_d==0)
             Nd(i_node)=N_sd*exp(-((iii_col-junction_l)*dx/decay_lens)^2)-...
                   N_body;
        elseif (dopslope_s==0 & dopslope_d~=0)
             Nd(i_node)=N_sd*exp(-((iii_col-junction_r)*dx/decay_l)^2)-...
                   N_body;




      for iii_col=junction_r:Nx
        i_node=iii_row*Nx+iii_col
        Nd(i_node)=N_sd-N_body



     if ox_pnt_flag==1

      Nd((Nx*t_topa+1):(Nx*t_topa+Nx))=...
          Nd((Nx*t_topa+1+Nx):(Nx*t_topa+2*Nx));
      Nd((Ntotal-Nx*(t_bota+1)+1):(Ntotal-Nx*(t_bota+1)+Nx))=...
          Nd((Ntotal-Nx*(t_bota+1)+1-Nx):(Ntotal-Nx*(t_bota+1)));



