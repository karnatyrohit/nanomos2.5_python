###########################################################################
##A python file function to save and plot the output of nanomos ###########
###########################################################################

from readinput import *
import datetime
from matplotlib.pyplot import *
import numpy as np
from globvars import globvars
from mpl_toolkits.mplot3d import Axes3D as ax
from matplotlib import cm



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
    if plot_Ec3d == 1:
        figure(6)
        [X, Y] = np.meshgrid(XI, YI)
        Z = trMEc
        ax = gca(projection = '3d')
        surf = ax.plot_surface(X, Y, Z, rstride=3, cstride=3, cmap=cm.coolwarm, linewidth=0.5, antialiased = True)

        #surf(XI,YI,trMEc)
        #shading interp commented out to reduce size of the .ps file
        title('3D Conduction band edge potential profile')
        ax.set_xlabel('X [nm]')
        ax.set_ylabel('Y [nm]')
        ax.set_zlabel('Ec [eV]')
        ax.view_init(elev=60, azim=50)
        ax.dist=8
        savefig('Ec_X_Y.png')

        XII = (0, XI)
        tem1 = (YI, trMEc)
        tem2 = (XII, tem1)
        #np.savetxt('Ec_X_Y.dat', tem2, fmt='%e', delimiter=';')
        #f1 = open('Ec_X_Y','w')
        #writer = csv.writer(f1, delimiter = ',')
        #writer.writerows(tem2)



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
    ################################################################################################
    if plot_Ne_IV==1:
    # SUBBAND CHARGE DENSITY (/cm^2) vs X for Diff. Vg
    # --------------------------------------------------------------
        if (Ng_step>=1 or (Ng_step==0 and Nd_step==0)):

            figure(2)
            for iii in np.arange(0,Ng_step+1):
                Ne1_temp = np.squeeze(Ne_sub[:, iii, Nd_step, :, :])
                #if len(np.shape(Ne1_temp)) == 1:
                Ne1 = Ne1_temp
                #else:
                #    Ne1 = np.sum(Ne1_temp,0)
                Ne2 = np.zeros((len(Ne1), Ng_step+1))
                Ne2[:,iii]=(Ne1)*1e-4
                plot(XI, Ne2[:, iii], 'r-')
                hold(True)
                grid(True)

            if Ng_step>=1:
                title('2D electron density along the channel at different Vg')
            elif Ng_step==0:
                title('2D electron density along the channel')
            xlabel('X [nm]')
            ylabel('N2D [cm^{-2}]')
            savefig('N2D_X1.png')

    #figure(11)
    #Ne1=sum(squeeze(Ne_sub(:,:,:,1,Nd_step+1)),3);
    #Ne2(:,1)=sum(Ne1,2)*1e-4;
    #semilogy(XI',Ne2(:,1),'r-');
    #hold on
    #grid on
    #if Ng_step>0
    #  for iii=2:Ng_step+1
    #    Ne1=sum(squeeze(Ne_sub(:,:,:,iii,Nd_step+1)),3);
    #    Ne2(:,iii)=sum(Ne1,2)*1e-4;
    #    semilogy(XI',Ne2(:,iii),'r-');
    #  end
    #end
    #title('2D electron density along the channel @Diff. VG');
    #xlabel('X [nm]');
    #ylabel('N2D [cm^{-2}]');
    #print -depsc ./output/LOG_N2D_X.ps

            temmm=[XI,Ne2]
            #save N2D_X1.dat temmm -ascii;


    # SUBBAND CHARGE DENSITY (/cm^2) vs X for Diff. Vd
    # --------------------------------------------------------------
        if Nd_step>=1:
            figure(3)
            for iii in np.arange(0,Nd_step+1):
                Ne1_temp = np.squeeze(Ne_sub[:, Ng_step, iii, :, :])
                #if len(np.size(Ne1_temp)) == 1:
                Ne1 = Ne1_temp
                #else:
                #Ne1 = np.sum(Ne1_temp,0)
                Ne2 = np.zeros((len(Ne1), Nd_step+1))
                Ne2[:, iii] = (Ne1)*1e-4
                plot(XI, Ne2[:, iii], 'r-')
                hold(True)
                grid(True)
            title('2D electron density along the channel at different Vd')
            xlabel('X [nm]')
            ylabel('N2D [cm^{-2}]')
            savefig('N2D_X2.png')

        temmm=[XI,Ne2]
        #save N2D_X2.dat temmm -ascii

    #******************************************************************************************
    if plot_Ec_IV==1:
    # The First SUBBAND ENERGY PROFILE vs X for Diff. Vg
    # ------------------------------------------------------

        if (Ng_step>=1 or (Ng_step==0 and Nd_step==0)):
            figure(4)
            for iii in np.arange(0,Ng_step+1):
                plot(XI, E_sub[0, iii, Nd_step, :, 0],'r-')
                hold(True)
                grid(True)

            if Ng_step>=1:
                title('The First Subband energy profile along the channel at different Vg')
            elif Ng_step==0:
                title('The First Subband energy profile along the channel')
            xlabel('X [nm]')
            ylabel('E_{SUB} [eV]')
            savefig('Ec_X1.png')

            tem = E_sub[0, iii, Nd_step, :, 0]
            sq_tem = np.squeeze(tem)
            temmm = [XI, sq_tem]
            #save Ec_X1.dat temmm -ascii

    # The First SUBBAND ENERGY PROFILE vs X for Diff. Vd
    # ------------------------------------------------------

        if Nd_step>=1:
            figure(5)
            for iii in np.arange(0,Nd_step+1):
                plot(XI,E_sub[0,Ng_step,iii, :, 0],'r-')
                hold(True)
                grid(True)

            title('The First Subband energy profile along the channel at different Vd')
            xlabel('X [nm]')
            ylabel('E_{SUB} [eV]')
            savefig('Ec_X2.png')

            tem = E_sub[0, Ng_step, :, :, 0]
            sq_tem = np.squeeze(tem)
            temmm = [XI,sq_tem]
            #save Ec_X2.dat temmm -ascii

    # 3D CHARGE DENSITY N(X,Y)
    # ------------------------------------------------------------
    if plot_Ne3d==1:
        figure(7)
        #surf(XI,YI,trMNe)
        #shading interp commented out to reduce size of the .ps file
        [X, Y] = np.meshgrid(XI, YI)
        Z = trMNe
        ax = gca(projection = '3d')
        surf = ax.plot_surface(X, Y, Z, rstride=3, cstride=3, cmap=cm.coolwarm, linewidth=0.5, antialiased = True)

        title('3D Electron density profile')
        ax.set_xlabel('X [nm]')
        ax.set_ylabel('Y [nm]')
        ax.set_zlabel('Ne [m^{-3}]')
        ax.view_init(elev=60, azim=50)
        ax.dist=8
        savefig('Ne_X_Y.png')

        XII = [0, XI]
        tem1 = [YI, trMNe]
        tem2 = [XII, tem1]
        #save Ne_X_Y.dat tem2 -ascii

    return
