###########################################################################
####################A function to evaluate F_prime once####################
#########################September 2001 - Purdue###########################
###########################################################################

from readinput import *
import numpy as np

def fprime():
    global  delta_x, Lsd, Lg_top, Lg_bot, t_top, t_bot, t_si, dx, dy
    global  eps_top, eps_bot, eps_si
    global  Nx, Ny, Ntotal
    global  F_prime, transport_model

    Lsda=round(Lsd/dx)
    Lg_topa=round(Lg_top/dx)
    Lg_bota=round(Lg_bot/dx)
    t_topa=round(t_top/dy)
    t_bota=round(t_bot/dy)
    t_sia=round(t_si/dy)
    ###########################################################################
    ####################Top gate insulator region##############################
    ###########################################################################
    for i_node=1:(Nx*(t_topa+1))

      if(i_node>=1 & i_node<=Lsda+1-1)
        F_prime(i_node,i_node)=1
        F_prime(i_node,i_node+Nx)=-1;
      elseif(i_node>=Lsda+1 & i_node<=((Lsda+Lg_topa)+1))
        F_prime(i_node,i_node)=1
      elseif(i_node>=((Lsda+Lg_topa)+1)+1 & i_node<=Nx)
        F_prime(i_node,i_node)=1
        F_prime(i_node,i_node+Nx)=-1
      elseif(i_node>=(Nx*t_topa+1) & i_node<=(Nx*(t_topa+1)))
        F_prime(i_node,i_node-Nx-1)=-eps_top/eps_si*dy/dx/8
        F_prime(i_node,i_node-Nx)=-eps_top/eps_si*(dx/dy-dy/dx/4)
        F_prime(i_node,i_node-Nx+1)=-eps_top/eps_si*dy/dx/8
        F_prime(i_node,i_node-1)=-(eps_top/eps_si+1)*dy/dx*3/8
        F_prime(i_node,i_node)=(eps_top/eps_si+1)*(dy/dx*3/4+dx/dy)
        F_prime(i_node,i_node+1)=-(eps_top/eps_si+1)*dy/dx*3/8
        F_prime(i_node,i_node+Nx-1)=-dy/dx/8
        F_prime(i_node,i_node+Nx)=-(dx/dy-dy/dx/4)
        F_prime(i_node,i_node+Nx+1)=-dy/dx/8
      else
        F_prime(i_node,i_node-Nx)=-eps_top/eps_si*dx/dy
        F_prime(i_node,i_node-1)=-eps_top/eps_si*dy/dx
        F_prime(i_node,i_node)=2*(dy/dx+dx/dy)*eps_top/eps_si
        F_prime(i_node,i_node+1)=-eps_top/eps_si*dy/dx
        F_prime(i_node,i_node+Nx)=-eps_top/eps_si*dx/dy



    %Bottom gate insulator region
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i_node=(Ntotal-Nx*(t_bota+1)+1):Ntotal
      if(i_node>=(Ntotal-Nx*(t_bota+1)+1) & i_node<=(Ntotal-Nx*t_bota))
        F_prime(i_node,i_node-Nx-1)=-dy/dx/8
        F_prime(i_node,i_node-Nx)=-(dx/dy-dy/dx/4)
        F_prime(i_node,i_node-Nx+1)=-dy/dx/8
        F_prime(i_node,i_node-1)=-(eps_bot/eps_si+1)*dy/dx*3/8
        F_prime(i_node,i_node)=(eps_bot/eps_si+1)*(dy/dx*3/4+dx/dy)
        F_prime(i_node,i_node+1)=-(eps_bot/eps_si+1)*dy/dx*3/8
        F_prime(i_node,i_node+Nx-1)=-eps_bot/eps_si*dy/dx/8
        F_prime(i_node,i_node+Nx)=-eps_bot/eps_si*(dx/dy-dy/dx/4)
        F_prime(i_node,i_node+Nx+1)=-eps_bot/eps_si*dy/dx/8
      elseif(i_node>=(Ntotal-Nx+1) & i_node<=(Ntotal-Nx+1+Lsda)-1)
        F_prime(i_node,i_node)=1
        F_prime(i_node,i_node-Nx)=-1
      elseif(i_node>=(Ntotal-Nx+1+Lsda) & i_node<=(Ntotal-Nx+1+Lsda+Lg_bota))
        F_prime(i_node,i_node)=1
      elseif(i_node>=(Ntotal-Nx+1+Lsda+Lg_bota)+1 & i_node<=Ntotal)
        F_prime(i_node,i_node)=1
        F_prime(i_node,i_node-Nx)=-1
      else
        F_prime(i_node,i_node-Nx)=-eps_bot/eps_si*dx/dy
        F_prime(i_node,i_node-1)=-eps_bot/eps_si*dy/dx
        F_prime(i_node,i_node)=2*(dx/dy+dy/dx)*eps_bot/eps_si
        F_prime(i_node,i_node+1)=-eps_bot/eps_si*dy/dx
        F_prime(i_node,i_node+Nx)=-eps_bot/eps_si*dx/dy




    %Specify the F_prime matrix in
    %the silicon film region
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i_node=(Nx*(t_topa+1))+1:(Ntotal-Nx*(t_bota+1)+1)-1
      F_prime(i_node,i_node-Nx)=-dx/dy
      F_prime(i_node,i_node-1)=-dy/dx
      F_prime(i_node,i_node)=2*(dx/dy+dy/dx)
      F_prime(i_node,i_node+1)=-dy/dx
      F_prime(i_node,i_node+Nx)=-dx/dy


    %Modify the F_prime matrix at
    %the right and left boundaries
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    i_node_l=1
    i_node_r=Nx
    for iii=1:Ny
      if(iii==1)
        F_prime(i_node_l,:)=0
        F_prime(i_node_l,i_node_l)=2
        F_prime(i_node_l,i_node_l+1)=-1
        F_prime(i_node_l,i_node_l+Nx)=-1
        F_prime(i_node_r,:)=0
        F_prime(i_node_r,i_node_r)=2
        F_prime(i_node_r,i_node_r-1)=-1
        F_prime(i_node_r,i_node_r+Nx)=-1
      elseif(iii>1 & iii<round((Nx*(t_topa+1))/Nx))
        F_prime(i_node_l,:)=0
        F_prime(i_node_l,i_node_l)=1
        F_prime(i_node_l,i_node_l+1)=-1
        F_prime(i_node_r,:)=0
        F_prime(i_node_r,i_node_r)=1
        F_prime(i_node_r,i_node_r-1)=-1
      elseif(iii>=round((Nx*(t_topa+1))/Nx) & iii<=round((Ntotal-Nx*t_bota)/Nx))

            F_prime(i_node_l,:)=0
            F_prime(i_node_l,i_node_l)=1
            F_prime(i_node_l,i_node_l+1)=-1
            F_prime(i_node_r,:)=0
            F_prime(i_node_r,i_node_r)=1
            F_prime(i_node_r,i_node_r-1)=-1

      elseif(iii>round((Ntotal-Nx*t_bota)/Nx) & iii<Ny)
        F_prime(i_node_l,:)=0
        F_prime(i_node_l,i_node_l)=1
        F_prime(i_node_l,i_node_l+1)=-1
        F_prime(i_node_r,:)=0
        F_prime(i_node_r,i_node_r)=1
        F_prime(i_node_r,i_node_r-1)=-1
      elseif(iii==Ny & ((Ntotal-Nx+1)<(Ntotal-Nx+1+Lsda)) &...
              ((Ntotal-Nx+1+Lsda+Lg_bota)<Ntotal))
        F_prime(i_node_l,:)=0
        F_prime(i_node_l,i_node_l)=2
        F_prime(i_node_l,i_node_l+1)=-1
        F_prime(i_node_l,i_node_l-Nx)=-1
        F_prime(i_node_r,:)=0
        F_prime(i_node_r,i_node_r)=2
        F_prime(i_node_r,i_node_r-1)=-1
        F_prime(i_node_r,i_node_r-Nx)=-1

      i_node_l=1+iii*Nx
      i_node_r=(1+iii)*Nx


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	END OF SPECIFIING F_prime
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%