###########################################################################
##A python file function to save and plot the output of nanomos ###########
###########################################################################

from readinput import *
import datetime
from matplotlib.pyplot import *
import numpy as np


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
    #save trans.dat Trans -ascii
    #save E.dat E -ascii

    fid = open("convergence.dat","w")

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
            grid()
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

    # ID-VD CHARACTERISTICS (A/m), SAVED TO "id_vd.dat"
    # -------------------------------------------------
        if (Nd_step>=1 and Ng_step==0):
            figure(1)
            plot(Vd_bias,Ie[0,:],'o-')
            grid()
            xlabel('V_{DS} [V]')
            ylabel('I_{DS} [\muA/\mum]')
            title('I_{DS} vs. V_{DS}')
            savefig('ID_VD.png')

            temmm = [Vd_bias, Ie[1,:]]
            np.savetxt('ID_VD.dat', temmm, fmt='%.6f', delimiter=';')
    #if plot_Iv==1 end

    return