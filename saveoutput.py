###########################################################################
##A python file function to save and plot the output of nanomos ###########
###########################################################################

from readinput import *
import datetime
from matplotlib.pyplot import *
import numpy as np
from globvars import globvars


def saveoutput(Ec,Ne,Ie,Ne_sub,E_sub,Te_sub,converge,Vd_temp):
    transport_model = transportmodel.value
    N_sd = Nsd1.value
    N_body = Nbody1.value
    fermi_flag = fermiflag1.value
    Vd = Vdc.value
    # Rohit- N_dos, Trans and E

    ###########################################################
    #SPECIFY XI, YI, Vg_bias, and Vd_bias for plotting use only
    ###########################################################

    Nx = round((2*Lsd+Lg_top)/dx)+1
    XI = np.linspace(-(Lsd+Lg_top/2)*1e9, (Lsd+Lg_top/2)*1e9, Nx)  # in nanometers
    Ny = round((t_top+t_bot+t_si)/dy)+1
    YI = np.linspace(-t_top*1e9, (t_si+t_bot)*1e9, Ny)
    Vg_bias = np.linspace(Vg1, Vg1+Vg_step*Ng_step, Ng_step+1)
    Vd_bias = np.linspace(Vd_temp, Vd_temp+Vd_step*Nd_step, Nd_step+1)

    ###########################################################
    #       POST PROCESSING
    ###########################################################
    MEc = np.reshape(Ec[Nd_step, :, Ng_step], (Ny, Nx)).transpose()
    trMEc = MEc.transpose()
    MNe = np.reshape(Ne[Nd_step, :, Ng_step], (Ny, Nx)).transpose()
    trMNe = MNe.transpose()

    #mkdir(dirname)
    #cd(dirname)

    #raw_data is always stored!!!!!
    #--------------------------------------------------
    #save output_rawdata
    #Save convergence data in a file
    #save DOS.dat N_dos -ascii
    Trans  = globvars.Trans
    #save trans.dat Trans -ascii
    E = globvars.E
    #save E.dat E -ascii
    N_dos = globvars.N_dos

    fid = open("convergence.dat", "w")

    fid.write('****************************************************************\n')
    fid.write('*******You ran a Nanomos simulation on the %s *********\n' % datetime.datetime.now())
    fid.write('********************* Using Nanomos 2.0 ************************\n')
    fid.write('********** The input file you used  was the following **********\n')
    fid.write('****************************************************************\n')

    #cd ..
    # fid1=fopen(filename,'r')
    #     while 1
    #         tline = fgetl(fid1)
    #         if ~ischar(tline), break, end
    #         count=fprintf(fid,'%s\n',tline);
    #     end
    # status=fclose(fid1)
    # cd(dirname)

    fid.write('****************************************************************\n')
    fid.write('********** Convergence data for simulated bias points **********\n')
    fid.write('****************************************************************\n')

    for i in np.arange(0,Ng_step+1):
        for j in np.arange(0,Nd_step+1):
            fid.write('\n')
            fid.write('**** Gate voltage = %f ****Drain voltage = %f' % (Vg_bias[i], Vd_bias[j]))
            fid.write('\n')
            for item in converge[i,j]:
                fid.write('%e\n' % item)

    fid.write('****************************************************************')
    fid.write('********** The end of the convergence data file*****************')
    fid.write('****************************************************************')

    fid.close()

    ############################################# IV plots ##############################################
    if plot_IV == 1:
    # ID-VD CHARACTERISTICS (A/m), SAVED TO "id_vd.dat"
    # -------------------------------------------------
        if (Ng_step>=1 and Nd_step>=1):
            Ng = np.size(Ie,0)
            figure(1)
            plot(Vd_bias,Ie[0,:],'o-')
            grid(True)
            hold(True)
            for i in np.arange(1,Ng):
                plot(Vd_bias, Ie[i, :], 'o-')

            xlabel('V_{DS} [V]')
            ylabel('I_{DS} [\muA/\mum]')
            title('I_{DS} vs. V_{DS}')
            savefig('ID_VD.png')

            temmm = Vd_bias
            Ie2 = Ie.transpose()
            fid = open('ID_VD.dat','w')
            ind = 0
            for item1 in temmm:
                fid.write("%e " % item1)
                for item2 in Ie2[ind,:]:
                    fid.write("%e " % item2)
                fid.write('\n')
                ind += 1
            fid.close()

    # ID-VG CHARACTERISTICS (A/m), SAVED TO "id_vg.dat"
    # -------------------------------------------------
    if (Ng_step>=1 and Nd_step==0):
        figure(1)
        semilogy(Vg_bias, Ie[:,0],'o-')
        grid(True)
        xlabel('V_{GS} [V]')
        ylabel('I_{DS} [\muA/\mum]')
        title('I_{DS} vs. V_{GS}')
        savefig('ID_VG.png')

        temmm=[Vg_bias, Ie[:,0]]
        np.savetxt( 'ID_VG.dat', temmm, fmt='%e', delimiter=';')

    # ID-VD CHARACTERISTICS (A/m), SAVED TO "id_vd.dat"
    # -------------------------------------------------
    if (Nd_step>=1 and Ng_step==0):
        figure(1)
        plot(Vd_bias, Ie[0,:], 'o-')
        grid(True)
        xlabel('V_{DS} [V]')
        ylabel('I_{DS} [\muA/\mum]')
        title('I_{DS} vs. V_{DS}')
        savefig('ID_VD.png')

        temmm = [Vd_bias, Ie[0,:]]
        np.savetxt('ID_VD.dat', temmm, fmt='%e', delimiter=';')
    #if plot_Iv==1 end

    #***************************************************************************************
    # Ec(X,Y)
    # -------------------------------------------
    """if plot_Ec3d==10:
        figure(6)
        surf(XI,YI,trMEc)
        #shading interp commented out to reduce size of the .ps file
        title('3D Conduction band edge potential profile')
        hx=xlabel('X [nm]')
        hy=ylabel('Y [nm]')
        hz=zlabel('Ec [eV]')
        savefig('Ec_X_Y.png')

        XII = [0, XI]
        tem1 = [YI, trMEc]
        tem2 = [XII, tem1]
        np.savetxt('Ec_X_Y.dat', tem2, fmt='%e', delimiter=';')"""

   #*******************************************************************************************
    if (plot_Ecsub==1 and max_subband>=1):
        figure(8)
        for iii in np.arange(0,max_subband):
            plot(XI, E_sub[0, Ng_step, Nd_step, :, iii],'r-')
            hold(True)
            grid(True)
            if (t_vall==3):
                plot(XI, E_sub[1, Ng_step, Nd_step, :,iii],'k-')
                plot(XI, E_sub[2, Ng_step, Nd_step, :,iii],'-')

        title('The Subbands energy profile along the channel')
        xlabel('X [nm]')
        ylabel('E_{SUB} [eV]')
        savefig('Ec_sub_X.png')

    ############################################################################################
    if (plot_Nesub==1 and max_subband>=1):
        figure(9)
        for iii in np.arange(0,max_subband):
            semilogy(XI, Ne_sub[0, Ng_step, Nd_step, :, iii],'r-')
            hold(True)
            grid(True)
            if (t_vall==3):
                semilogy(XI, Ne_sub[1, Ng_step, Nd_step, :, iii],'k-')
                semilogy(XI, Ne_sub[2, Ng_step, Nd_step, :, iii],'-')
        title('2D electron density of the subbands along the channel ')
        xlabel('X [nm]')
        ylabel('N2D [cm^{-2}]')

        savefig('Ne_sub_X.png')

    if (Ng_step==0 and Nd_step==0):
        #	PLOT the graph of Transmission coefficient versus Energy
        #------------------------------------------------------------------
        if transport_model==4:
            figure(11)
            plot(Trans, E)
            xlabel('Transmission Coefficient')
            ylabel('Energy (eV)')
            savefig('Trans.png')
        #	PLOT the Density of states versus energy along the device
        #------------------------------------------------------------------
        if transport_model==4:
            figure(12)
            pcolor(XI, E, np.sqrt(N_dos[:,0:len(E)].transpose()))
            #shading interp
            # Square root of the DOS to get better visualization
            xlabel('Distance along the device (nm)')
            ylabel('Energy (eV)')
            savefig('DOS.png')
            #print -depsc2 DOS.ps


    return
