# note: python3...
#
import numpy
import pylab as plt
import random
import scipy
import scipy.optimize as spo
import math
import functools
#
import scipy.fftpack
# note might also use numpy.fft ???
#
pi2 = math.pi*2.0		# can't name it 2pi
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
def fft_demo1(n_modes=2, n_cycles=5, points_per_cycle=500, amplitudes=None, phis=None, omega_factor=10., fignum=0):
	# some basic FFT experiment and play-stuff, to remember how to do FFTs and to figure out how to do cool stuff with them.
	# #so it might be possible to use the complex part of the amplitudes to phase-shift, but you can't just take the real/Im part to do that.
	#
	# make a couple synthetic sequences; take FFT, reconstruct, subtract modes, etc.
	#
	fs_label = 18.
	#
	amplitudes = (amplitudes or [1.0 for j in range(n_modes)])
	phis       = (phis or [0. for j in range(n_modes)])
	omega_factor = (omega_factor or float(n_cycles)*pi2)
	R = random.Random()
	#
	Ts = numpy.arange(0., n_cycles*pi2, pi2/points_per_cycle)
	#
	datas = []
	for k, (a,phi,omega) in enumerate(zip(amplitudes, phis, omega_factor*numpy.random.random(n_modes))):
		datas += [{'A':a, 'phi':phi, 'omega':omega, 'data':numpy.array(numpy.cos(Ts*omega + phi)), 'data2':numpy.array(numpy.cos(Ts*omega))}]
	#
	D_cum = functools.reduce(lambda a,b: a+b, [D['data'] for D in datas])
	D_cum2 = functools.reduce(lambda a,b: a+b, [D['data2'] for D in datas])
	#
	# now, get FFT of D_cum:
	#y_fft = numpy.fft.fft(D_cum)[:.5*len(D_cum)]
	#x_fft = numpy.linspace(0., 1.0/(2.*n_cycles*pi2), n_cycles*points_per_cycle/2.)
	y_fft = numpy.fft.fft(D_cum)
	x_fft = numpy.linspace(0., 1.0/(n_cycles*pi2), n_cycles*points_per_cycle)
	#
	y_prime = numpy.fft.ifft(numpy.real(y_fft))	# this yields an excellend reproduction.
	#
	# reocnstruct from only dominant phase:
	fft_basis_prime = [y if y==max(y_fft) else 0. for y in y_fft]
	y_single_phase = numpy.fft.ifft(fft_basis_prime)
	#
	plt.figure(fignum)
	plt.clf()
	#
	for D in datas:
		plt.plot(Ts/pi2, D['data'], '-', label='$A:%.2f, phi:%.2f, \\omega: %.2f$' % (D['A'], D['phi'], D['omega']), lw=1.5, zorder=4, alpha=.5 )
	plt.plot(Ts/pi2, D_cum, '.-', lw=2.5, zorder=5, alpha=.9, label='cum.')
	plt.legend(loc=0, numpoints=1)
	plt.xlabel('$x$ (or $t$) in factors of $pi$ radians', size=fs_label)
	plt.ylabel('$y=f(x)$', size=fs_label)
	#
	plt.figure(fignum+1)
	plt.clf()
	plt.plot(range(len(y_fft)),  numpy.real(y_fft), '.-', lw=2., label='fft_real')
	plt.plot(range(len(y_fft)),  numpy.imag(y_fft), '.-', lw=2., label='fft_imag')
	plt.legend(loc=0, numpoints=1)
	plt.xlabel('frequency $\\omega$ or wave number $k \\sim 1/x$', size=fs_label)
	plt.xlabel('fft amplitude $A$', size=fs_label)	
	#
	# reconstruct from dominant phase:
	plt.figure(fignum+2)
	plt.clf()
	#
	max_0 = max([abs(x) for x in D_cum])
	max_sp = max([abs(x) for x in y_single_phase])
	max_ratio = max_0/max_sp
	plt.plot(Ts/pi2, D_cum, '.-', lw=2.5, zorder=5, alpha=.9, label='cum.')
	plt.plot(Ts/pi2, [max_ratio * y for y in y_single_phase], '.-', lw=2, label='single_phase')
	plt.legend(loc=0, numpoints=1)
	'''
	#plt.plot(Ts/pi2, numpy.array(y_prime)-D_cum, '.-')
	plt.plot(Ts/pi2, numpy.array(y_prime), 'b.-')
	plt.plot(Ts/pi2, numpy.real(numpy.array(y_prime)), 'r-')
	plt.plot(Ts/pi2, numpy.imag(numpy.array(y_prime)), 'c-')
	plt.plot(Ts/pi2, D_cum, 'g.-', lw=1.5, zorder=5, alpha=.9, label='cum.')
	plt.plot(Ts/pi2, D_cum2, 'm.-', lw=1.5, zorder=5, alpha=.9, label='cum.')
	'''

def ffdemo2(phases=[0., .25, .5, 1.], fignum=0, n_cycles=4., n_per_cycle=500):
	# demo to examine phase and complex component of fft amplitudes.
	#
	#
	my_axes=[]
	dy=(1.0/len(phases))-.02
	#
	figses = []
	for fn in range(fignum, fignum+2):		
		f=plt.figure(fn)
		plt.clf()
		for k in range(len(phases)):
			f.add_axes([.03, .02 + k*dy, .45, dy])
			if k==0:
				f.add_axes([.5, .02 + k*dy, .45, dy])
			else:
				f.add_axes([.5, .02 + k*dy, .45, dy], sharex=f.axes[1+2*(k-1)], sharey=f.axes[1+2*(k-1)])
		plt.draw()
		#
		figses+=[f]
		
	#
	X = [x*pi2/float(n_per_cycle) for x in numpy.arange(n_cycles*n_per_cycle)]
	X_fft = numpy.linspace(0., 1.0/(n_cycles*pi2), n_cycles*n_per_cycle)
	As = []
	#
	for k, phi in enumerate(phases):
		Y = [math.sin(x + phi) for x in X]
		A_fft = numpy.fft.fft(Y)
		As += [A_fft]
		#	
		#
		figses[0].axes[2*k].plot(X,Y, '-', lw=2.0)
		figses[0].axes[2*k+1].plot(X_fft[0:len(X_fft)/2], numpy.real(A_fft)[0:len(X_fft)/2], 'b-', lw=2.0)
		figses[0].axes[2*k+1].plot(X_fft[0:len(X_fft)/2], numpy.imag(A_fft)[0:len(X_fft)/2], 'g-', lw=2.0)
		#
	#for ax in figses[0].axes: ax.draw()
	plt.show()
	#
	# now, reconstruct from the fourier series. how do we eliminate the phase?
	#
	plt.figure(1)
	for k, A in enumerate(As):
		A_prime = [math.sqrt(a*numpy.conj(a)) for a in A]
		Y_prime = numpy.fft.ifft(A_prime)
		#
		#N_len = len(list(X_fft))
		figses[1].axes[2*k].plot(X,Y_prime, '-', lw=2.0)
		figses[1].axes[2*k+1].plot(X_fft, A_prime, 'b-', lw=2.0)
	plt.show()
	

if __name__=='__main__':
	pass
else:
	plt.ion()
