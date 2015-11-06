# note: python3...
#
import numpy
import pylab as plt
#
import scipy.fftpack
# note might also use numpy.fft ???
#
def test1():
	N = 600
	# sample spacing
	T = 1.0 / 800.0
	x = numpy.linspace(0.0, N*T, N)
	y = numpy.sin(50.0 * 2.0*numpy.pi*x) + 0.5*numpy.sin(80.0 * 2.0*numpy.pi*x)
	yf = scipy.fftpack.fft(y)
	xf = numpy.linspace(0.0, 1.0/(2.0*T), N/2)

	fig, ax = plt.subplots()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.plot(xf, 2.0/N * numpy.abs(yf[0:N/2]), '.-')
#
#
def fft_demo1():
	# some basic FFT experiment and play-stuff, to remember how to do FFTs and to figure out how to do cool stuff with them.
	#
	pass

if __name__=='__main__':
	pass
else:
	plt.ion()
