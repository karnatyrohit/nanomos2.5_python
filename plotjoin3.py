import pylab

list_of_files = [('N2D_X1.dat', 'python'), ('N2D_X1m.dat', 'Matlab')]


datalist = [ ( pylab.loadtxt(filename), label ) for filename, label in list_of_files]

i=1
for data, label in datalist:
    if i==1:
        pylab.plot( data[:,0], data[:,1],'b',label = label)
        i=i+1
    else:
        pylab.plot( data[:,0], data[:,1],'ro', label = label)


pylab.legend()
pylab.title('2D electron density along the channel')
pylab.xlabel('X [nm]')
pylab.ylabel('N2D $[cm^{-2}]$')
pylab.plt.savefig('N2Dx1j.png')
pylab.plt.show()
