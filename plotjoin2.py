import pylab

list_of_files = [('N2D_X2.dat', 'python'), ('N2D_X2m.dat', 'Matlab')]


datalist = [ ( pylab.loadtxt(filename), label ) for filename, label in list_of_files]

i=1
for data, label in datalist:
    if i==1:
        for ind in range(1,13):
            pylab.plot( data[:,0], data[:,1:84],'b')
        i=i+1
    else:
        for ind in range(1,13):
            pylab.plot( data[:,0], data[:,1:84],'r--')


pylab.legend()
pylab.title('2D electron density along the channel at different Vd')
pylab.xlabel('X [nm]')
pylab.ylabel('N2D [$cm^{-2}$]')
pylab.plt.savefig('N2Dx2j.png')
pylab.plt.show()
