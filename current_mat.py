#****************************************************************************************
#***************CALCULATE THE CURRENT MATRIX Ii(muj) GIVEN TRANSMISSIONS*****************
#***********************and Fermi energies at each scatters******************************
#*****************************(Zhibin Ren 10-12-00)**************************************
#****************************************************************************************

from readinput import *
from globvars import  globvars
from fermi import fermi
from scipy import sparse
from scipy.sparse.linalg import spsolve

def current_mat(mu_scatter_old, T_E, E):
    #**************************FUNDAMENTAL physical constants********************************
    fermi_flag = fermiflag1.value
    Nx = globvars.Nx
    Temp = Te
    nu_scatter = globvars.nu_scatter
    mx = globvars.mx
    my = globvars.my
    mz = globvars.mz

    #**********************************INITIALIZATION****************************************
    dmu = 1e-6
    I_criterion = 1e-2 #1E-8A/um

    Iin = np.zeros((nu_scatter, 1))
    mu_scatter_new = np.zeros((nu_scatter+2, 1))

    delta_mu = np.zeros((nu_scatter+2, 1))
    delta_mu = np.squeeze(delta_mu)
    I_tem = np.zeros(2)
    IMU = np.zeros((nu_scatter, nu_scatter))
    IMU_dummy = np.zeros((nu_scatter, nu_scatter))
    T_dummy = np.sum(T_E, 0)
    #print sparse.csr_matrix(T_dummy[:,0])
    #print mu_scatter_old

    #***************************************************************************************
    #*************************************EVALUATE Iin**************************************
    #***************************************************************************************
    for i_node in np.arange(0,nu_scatter):
        I_dummy1 = np.dot(fermi(((mu_scatter_old[i_node]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0),T_dummy[:,i_node])
        #print fermi(((mu_scatter_old[i_node]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0)
        #print I_dummy1
        #print T_dummy[:,i_node]
        I_dummy2 = 0
        #print np.shape(fermi(((mu_scatter_old[i_node]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0))
        #print np.shape(T_dummy[:,i_node])
        #Idum = fermi(((mu_scatter_old[i_node]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0) * np.reshape(T_dummy[:,i_node],(799,1))
        #print Idum
        #exit()
        for j_node in np.arange(0,nu_scatter+2):
            I_dummy2 = I_dummy2 + np.dot(fermi(((mu_scatter_old[j_node]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0),T_E[j_node, :, i_node])
            #print I_dummy2
        Iin[i_node] = I_dummy1-I_dummy2
        #print I_dummy1
        #print 'idum2'
        #print I_dummy2
        #print Iin[i_node]
        #exit()
    Isc = np.max(abs(Iin))
    #print Iin

    #***************************************************************************************
    #*************************************EVALUATE IMU**************************************
    #***************************************************************************************

    if Isc>=I_criterion:
        for i_node in np.arange(0, nu_scatter):
            IMU_dummy[i_node,i_node]=np.dot(((fermi(((mu_scatter_old[i_node]+dmu-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0)-fermi(((mu_scatter_old[i_node]-dmu-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0))/(2.0*dmu)),T_dummy[:,i_node])
            for j_node in np.arange(0, nu_scatter):
                IMU[i_node,j_node]=np.dot(((fermi(((mu_scatter_old[j_node]+dmu-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0)-fermi(((mu_scatter_old[j_node]-dmu-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0))/(2.0*dmu)),T_E[j_node, :,i_node])
        IMU=IMU_dummy-IMU

    #***************************************************************************************
    #**********************************END OF EVALUATE IMU**********************************
    #***************************************************************************************

    mu_scatter_new = mu_scatter_old

    #*********************************Newton searching loop*********************************

    iiii = 0
    #print Isc
    #print I_criterion
    while(Isc>=I_criterion):
        print 'Entering Jacobian loop in Current_mat'
        delta_mu[0:nu_scatter] = -spsolve(sparse.csr_matrix(IMU),sparse.csr_matrix(Iin))

        # The following IF statement performs a check to insure that the correction does not
        # bring the fermi level in the device above the source fermi level or below the drain fermi level.
        # If it does, then we averaged the fermi level of the scatterer to make it fit within that physical range of fermi levels.
        # If we don't perform that check, in some cases we can get a scatterer fermi level well below the drain fermi level
        # which forces charges to flow from the drain to fill the available state.......as a result -> NO convergence.

        mu_scatter_new=mu_scatter_new+delta_mu

        for i_node in np.arange(0, nu_scatter):
            if mu_scatter_new[i_node]>mu_scatter_new[nu_scatter]:
                mu_scatter_new[i_node] = mu_scatter_new[nu_scatter]
            elif mu_scatter_new[i_node]<mu_scatter_new[nu_scatter+1]:
                mu_scatter_new[i_node]=mu_scatter_new[nu_scatter+1]
            else:
                mu_scatter_new[i_node]=mu_scatter_new[i_node]

        #if (max(mu_scatter_new(1:nu_scatter))>mu_scatter_new(nu_scatter+1) ...
        #| min(mu_scatter_new(1:nu_scatter))<mu_scatter_new(nu_scatter+2))

        #   mu_scatter_new(1:nu_scatter)=(mu_scatter_old(1:nu_scatter)+mu_scatter_new(1:nu_scatter))/2.0;
        #   fprintf(1,'#s#s#e\n','AVERAGED ','MAX CHANGE ',max(abs(mu_scatter_new(1:nu_scatter)-mu_scatter_old(1:nu_scatter))));

        #else
        #    fprintf(1,'#s#s#e\n','NOT AVERAGED ','MAX CHANGE ',max(abs(mu_scatter_new(1:nu_scatter)-mu_scatter_old(1:nu_scatter))));

        #end

        #for i nu_scatter

        #  for i_node=1:nu_scatter
        #    if(abs(delta_mu(i_node))<=1)
        #       delta_mu(i_node)=delta_mu(i_node);
        #    elseif(1<abs(delta_mu(i_node)) & abs(delta_mu(i_node)) <3.7)
        #       delta_mu(i_node)=sign(delta_mu(i_node))*power(abs(delta_mu(i_node)),1/5);
        #    elseif(abs(delta_mu(i_node))>=3.7)
        #       delta_mu(i_node)=sign(delta_mu(i_node))*log(abs(delta_mu(i_node)));
        #    end
        #  end

        #  mu_scatter_new=mu_scatter_new+delta_mu;

        #****************************************************************************************
        #************************************ EVALUATE Iin***************************************
        #****************************************************************************************
        for i_node in np.arange(0,nu_scatter):
            I_dummy1 = np.dot(fermi(((mu_scatter_new[i_node]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0),T_dummy[:,i_node])
            I_dummy2=0
            for j_node in np.arange(0,nu_scatter+2):
                I_dummy2=I_dummy2+np.dot(fermi(((mu_scatter_new[j_node]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0),T_E[j_node, :,i_node])
            Iin[i_node] = I_dummy1-I_dummy2
        Isc = max(abs(Iin))
        print 'Isc =', Isc

        #****************************************************************************************
        #***********************************END OF EVALUATE Iin**********************************
        #****************************************************************************************

        if Isc>=I_criterion:

        #****************************************************************************************
        #**************************************EVALUATE IMU**************************************
        #****************************************************************************************
            for i_node in np.arange(0,nu_scatter):
                IMU_dummy[i_node,i_node]=np.dot(((fermi(((mu_scatter_new[i_node]+dmu-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0)-fermi(((mu_scatter_new[i_node]-dmu-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0))/(2.0*dmu)),T_dummy[:,i_node])
                for j_node in np.arange(0,nu_scatter):
                    IMU[i_node,j_node]=np.dot(((fermi(((mu_scatter_new[j_node]+dmu-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0)-fermi(((mu_scatter_new[j_node]-dmu-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0))/(2.0*dmu)),T_E[j_node, :,i_node])
            IMU=IMU_dummy-IMU

            #****************************************************************************************
            #***********************************END OF EVALUATE IMU**********************************
            #****************************************************************************************
        iiii = iiii+1
        print 'iiii = ', iiii

        # Copy old vals
        mu_scatter_old=mu_scatter_new

    for i_node in np.arange(0, 2):
        I_dummy1 = np.dot(fermi(((mu_scatter_new[i_node+nu_scatter]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0),T_dummy[:,i_node+nu_scatter])
        #print I_dummy1
        #exit()
        I_dummy2=0
        for j_node in np.arange(0, nu_scatter+2):
            I_dummy2 = I_dummy2+np.dot(fermi(((mu_scatter_new[j_node]-E)/(k_B*Temp/q)),fermi_flag,-1.0/2.0),T_E[j_node, :, i_node+nu_scatter])
        I_tem[i_node]=I_dummy1-I_dummy2

    Is=I_tem[0]
    #print Is
    Id=I_tem[1]
    #print Id



    return [Is, Id, Iin, mu_scatter_new]
#******************************************************************************************
#*************************** THE END OF FUNCTION CURRENT_MAT*******************************
#******************************************************************************************
