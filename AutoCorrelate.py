import scipy
import itertools
from matplotlib import pylab as plt

class AutoCorrelate:
	class Trace:
		def __init__(self, trace, lag):
			self.trace   = trace
			self.lag     = lag
			self.mew     = scipy.mean(trace)
			self.var     = scipy.var (trace)
			self.xy_x_y  = AutoCorrelate._xy_x_y_(trace)
			self.weights = scipy.arange(len(trace), 1, -1)

		def __len__(self):
			return (len(self.trace))

	class SubTrace(Trace):
		def __init__(self, trace, lag):
			assert(lag >= trace.lag    )
			assert(lag %  trace.lag ==0)
			self.step   = trace.lag
			self.trace  = trace.trace
			self.lag    = lag
			self.mew    = trace.mew
			self.var    = trace.var
			self.xy_x_y = trace.xy_x_y [::int(lag/trace.lag)]
			self.weights= trace.weights[::int(lag/trace.lag)]



	@staticmethod
	def global_mean_global_variance(traces): 
		##Calculates mean and variance from individual mean/variance using the definitions of mean and variance
		##http://stats.stackexchange.com/questions/43159/how-to-calculate-pooled-variance-of-two-groups-given-known-group-variances-mean
		longest = max( (len(trace) for trace in traces) )
		length  = sum([len(t)                          for t in traces])
		mew     = sum([len(t) *  t.mew                 for t in traces]) / length
		sigma   = sum([len(t) * (t.mew ** 2 + t.var)   for t in traces]) / length - mew ** 2

		decoded = scipy.array([AutoCorrelate._foil_(trace.xy_x_y, trace.weights, mew) for trace in traces])
		decoded = scipy.array([scipy.pad(ac, (0,longest-len(ac)),'constant')  for ac in decoded])
		return(scipy.sum(decoded, axis=0) / sigma / length)

	@staticmethod
	def local_mean_local_variance(traces):
		longest = max((len(trace.weights) for trace in traces))
		decoded = scipy.array([AutoCorrelate._foil_(trace.xy_x_y, trace.weights, trace.mew) / (trace.var * len(trace)) for trace in traces])
		decoded = scipy.array([scipy.pad(ac           , (0, longest - len(ac           )), 'constant')  for ac in decoded])
		weights = scipy.array([scipy.pad(trace.weights, (0, longest - len(trace.weights)), 'constant') for trace in traces])
		return(scipy.average(decoded, weights=weights , axis=0))


	@staticmethod
	def _lag_(array, lag):
		"""
		Calculates the triplet representation of a trace for a single lag time.
		[sum(X*Y),sum(X),sum(Y)]
		"""
		n = len(array)
		return (sum(array[:n-lag]*array[lag:n]),
				sum(array[:n-lag]),
				sum(array[lag:n]))   

	@staticmethod
	def _xy_x_y_(array):
		"""
		Calculates the triplet representation of a trace for all lag times.
		[[sum(X*Y),sum(X),sum(Y)],... for each lag]
		"""
		n = len(array)
		return(scipy.array(list(map(lambda i:AutoCorrelate._lag_(array, i) ,range(n-1)))))

	@staticmethod
	def _foil_(XY_X_Y, weights, mew):
		"""
		Given a triplet representation of a trace for all lag times Calculates
		the Autocorrelation.
		"""
		XY=          XY_X_Y[:,0]
		UX=mew     * XY_X_Y[:,1]
		UY=mew     * XY_X_Y[:,2]
		UU=weights * mew ** 2
		return(XY-UX-UY+UU)

def __sample__(data_size, samples = 1):
	##ONLEY WORKS ON SCIPY ARRAYS
	samples = (scipy.random.randint(0,data_size, size=data_size) for i in range(samples))
	for sample in samples:
		yield (sample)

"""
https://stackoverflow.com/questions/47850760/using-scipy-fft-to-calculate-autocorrelation-of-a-signal-gives-different-answer?rq=1
import scipy.fftpack
def autocorrelation(x) :
    xp = scipy.fftpack.ifftshift((x - scipy.average(x))/scipy.std(x))
    n,  = xp.shape
    xp = scipy.r_[xp[:n//2], scipy.zeros_like(xp), xp[n//2:]]
    f  = scipy.fftpack.fft(xp)
    p  = scipy.absolute(f)**2
    pi = scipy.fftpack.ifft(p)
    return scipy.real(pi)[:n//2]/(scipy.arange(n//2)[::-1]+n//2)

#np.arange(n,0,-1)[:n//2]
#scipy.arange(n//2)[::-1]+n//2

# def autocorrelation(x):
#     xp = ifftshift((x - np.average(x))/np.std(x))
#     n, = xp.shape
#     xp = np.r_[xp[:n//2], np.zeros_like(xp), xp[n//2:]]
#     f = fft(xp)
#     p = np.absolute(f)**2
#     pi = ifft(p)
#     return np.real(pi)[:n//2]/(np.arange(n//2)[::-1]+n//2)


ac = [autocorrelation(trace) for trace in traces[selection]]
x  = scipy.array([scipy.pad(c,pad_width=(0,360 - len(c)),mode='constant', constant_values=scipy.nan) for c in ac])
#plt.plot(scipy.mean(scipy.apply_along_axis(autocorrelation,axis=1,arr=padded),axis=0))
plt.plot(scipy.nanmean(x, axis=0),c='red')
plt.plot(AutoCorrelate.AutoCorrelate.local_mean_local_variance(atraces))
plt.xlim([0,150])
"""