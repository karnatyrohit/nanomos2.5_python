import pylab

list_of_files = [('Ec_X2.dat', 'python'), ('Ec_X2m.dat', 'Matlab')]


datalist = [ ( pylab.loadtxt(filename), label ) for filename, label in list_of_files]

i=1
for data, label in datalist:
    if i==1:
        for ind in range(0,13):
            pylab.plot( data[:,0], data[:,1:84],'b')
        i=i+1
    else:
        for ind in range(0,13):
            pylab.plot( data[:,0], data[:,1:84],'r--')


pylab.legend()
pylab.title('The First Subband energy profile along the channel at different Vd')
pylab.xlabel('X [nm]')
pylab.ylabel('$E_{SUB}$ [eV]')
pylab.plt.savefig('Ecx2j.png')
pylab.plt.show()