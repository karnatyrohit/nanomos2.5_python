###########################################################################
####################A function to evaluate F_prime once####################
#########################September 2001 - Purdue###########################
###########################################################################

from readinput import *


def fprime(Nx, Ny, Ntotal, F_prime):

    Lsda = round(Lsd/dx)
    Lg_topa = round(Lg_top/dx)
    Lg_bota = round(Lg_bot/dx)
    t_topa = round(t_top/dy)
    t_bota = round(t_bot/dy)
    t_sia = round(t_si/dy)

    ###########################################################################
    ####################Top gate insulator region##############################
    ###########################################################################
    for i_node in range(0,(Nx*(t_topa+1))):
        if i_node >= 0 and i_node < Lsda:
            F_prime[i_node,i_node] = 1
            F_prime[i_node,i_node+Nx] = -1
        elif (i_node >= Lsda and i_node < ((Lsda+Lg_topa)+1)):
            F_prime[i_node,i_node] = 1
        elif (i_node >= ((Lsda+Lg_topa)+1) and i_node < Nx):
            F_prime[i_node,i_node] = 1
            F_prime[i_node,i_node+Nx] = -1
        elif(i_node >= (Nx*t_topa) and i_node < (Nx*(t_topa+1))):
            F_prime[i_node,i_node-Nx-1] = -eps_top/eps_si*dy/dx/8
            F_prime[i_node,i_node-Nx] = -eps_top/eps_si*(dx/dy-dy/dx/4)
            F_prime[i_node,i_node-Nx+1] = -eps_top/eps_si*dy/dx/8
            F_prime[i_node,i_node-1] = -(eps_top/eps_si+1)*dy/dx*3/8
            F_prime[i_node,i_node] = (eps_top/eps_si+1)*(dy/dx*3/4+dx/dy)
            F_prime[i_node,i_node+1] = -(eps_top/eps_si+1)*dy/dx*3/8
            F_prime[i_node,i_node+Nx-1] = -dy/dx/8
            F_prime[i_node,i_node+Nx] = -(dx/dy-dy/dx/4)
            F_prime[i_node,i_node+Nx+1] = -dy/dx/8
        else:
            F_prime[i_node,i_node-Nx] = -eps_top/eps_si*dx/dy
            F_prime[i_node,i_node-1] = -eps_top/eps_si*dy/dx
            F_prime[i_node,i_node] = 2*(dy/dx+dx/dy)*eps_top/eps_si
            F_prime[i_node,i_node+1] = -eps_top/eps_si*dy/dx
            F_prime[i_node,i_node+Nx] = -eps_top/eps_si*dx/dy

    # Bottom gate insulator region
    ###########################
    for i_node in range((Ntotal-Nx*(t_bota+1)),Ntotal):
        if(i_node >= (Ntotal-Nx*(t_bota+1)) and i_node < (Ntotal-Nx*t_bota)):
            F_prime[i_node,i_node-Nx-1] = -dy/dx/8
            F_prime[i_node,i_node-Nx] = -(dx/dy-dy/dx/4)
            F_prime[i_node,i_node-Nx+1] = -dy/dx/8
            F_prime[i_node,i_node-1] = -(eps_bot/eps_si+1)*dy/dx*3/8
            F_prime[i_node,i_node] = (eps_bot/eps_si+1)*(dy/dx*3/4+dx/dy)
            F_prime[i_node,i_node+1] = -(eps_bot/eps_si+1)*dy/dx*3/8
            F_prime[i_node,i_node+Nx-1] = -eps_bot/eps_si*dy/dx/8
            F_prime[i_node,i_node+Nx] = -eps_bot/eps_si*(dx/dy-dy/dx/4)
            F_prime[i_node,i_node+Nx+1] = -eps_bot/eps_si*dy/dx/8
        elif(i_node >= (Ntotal-Nx) and i_node < (Ntotal-Nx+Lsda)):
            F_prime[i_node,i_node] = 1
            F_prime[i_node,i_node-Nx] = -1
        elif(i_node >=  (Ntotal-Nx+Lsda) and i_node < (Ntotal-Nx+1+Lsda+Lg_bota)):
            F_prime[i_node,i_node] = 1
        elif(i_node >= (Ntotal-Nx+1+Lsda+Lg_bota) and i_node < Ntotal):
            F_prime[i_node,i_node] =1
            F_prime[i_node,i_node-Nx] = -1
        else:
            F_prime[i_node,i_node-Nx] = -eps_bot/eps_si*dx/dy
            F_prime[i_node,i_node-1] = -eps_bot/eps_si*dy/dx
            F_prime[i_node,i_node] = 2*(dx/dy+dy/dx)*eps_bot/eps_si
            F_prime[i_node,i_node+1] = -eps_bot/eps_si*dy/dx
            F_prime[i_node,i_node+Nx] = -eps_bot/eps_si*dx/dy

    # Specify the F_prime matrix in
    # the silicon film region
    ###########################
    for i_node in range((Nx*(t_topa+1)),(Ntotal-Nx*(t_bota+1))):
        F_prime[i_node,i_node-Nx] = -dx/dy
        F_prime[i_node,i_node-1] = -dy/dx
        F_prime[i_node,i_node] = 2*(dx/dy+dy/dx)
        F_prime[i_node,i_node+1] = -dy/dx
        F_prime[i_node,i_node+Nx] = -dx/dy

    # Modify the F_prime matrix at
    # the right and left boundaries
    ###########################
    i_node_l = 0
    i_node_r = Nx - 1
    for iii in range(0, Ny):
        if iii == 0:
            F_prime[i_node_l, :] = 0
            F_prime[i_node_l, i_node_l] = 2
            F_prime[i_node_l, i_node_l+1] = -1
            F_prime[i_node_l, i_node_l+Nx] = -1
            F_prime[i_node_r, :] = 0
            F_prime[i_node_r, i_node_r] = 2
            F_prime[i_node_r, i_node_r-1] = -1
            F_prime[i_node_r, i_node_r+Nx] = -1
        elif(iii > 0 and iii < round((Nx*(t_topa))/Nx)):
            F_prime[i_node_l, :] = 0
            F_prime[i_node_l, i_node_l] = 1
            F_prime[i_node_l, i_node_l+1] = -1
            F_prime[i_node_r, :] = 0
            F_prime[i_node_r, i_node_r] = 1
            F_prime[i_node_r, i_node_r-1] = -1
        elif(iii >= round((Nx*(t_topa))/Nx) and iii < round((Ntotal-Nx*t_bota)/Nx)):
            F_prime[i_node_l, :] = 0
            F_prime[i_node_l, i_node_l] = 1
            F_prime[i_node_l, i_node_l+1] = -1
            F_prime[i_node_r, :] = 0
            F_prime[i_node_r, i_node_r] = 1
            F_prime[i_node_r, i_node_r-1] = -1
        elif(iii >= round((Ntotal-Nx*t_bota)/Nx) and iii<Ny):
            F_prime[i_node_l, :] = 0
            F_prime[i_node_l, i_node_l] = 1
            F_prime[i_node_l, i_node_l+1] = -1
            F_prime[i_node_r, :] = 0
            F_prime[i_node_r, i_node_r] = 1
            F_prime[i_node_r, i_node_r-1] = -1
        elif iii == Ny-1 and ((Ntotal-Nx+1) < (Ntotal-Nx+1+Lsda)) and ((Ntotal-Nx+1+Lsda+Lg_bota) < Ntotal):
            F_prime[i_node_l, :] = 0
            F_prime[i_node_l, i_node_l] = 2
            F_prime[i_node_l, i_node_l+1] = -1
            F_prime[i_node_l, i_node_l-Nx] = -1
            F_prime[i_node_r, :] = 0
            F_prime[i_node_r, i_node_r] = 2
            F_prime[i_node_r, i_node_r-1] = -1
            F_prime[i_node_r, i_node_r-Nx] = -1

        i_node_l = 1+iii*Nx
        i_node_r = (1+iii)*Nx

    return F_prime
    #####################################################
    #	END OF SPECIFYING F_prime
    #####################################################