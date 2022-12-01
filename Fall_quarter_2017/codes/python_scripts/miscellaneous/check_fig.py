import matplotlib as mpl
mpl.use('Agg')
import pylab as py

py.figure()
# ax = fig.add_subplot(111)
py.plot(range(10,20),range(10))
py.savefig('temp.pdf')

