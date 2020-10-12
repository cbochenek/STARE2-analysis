import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sigpyproc.Readers import FilReader
from scipy import fftpack
import scipy.ndimage as snd

def lorentzian(x,a,x0,gamma):        
    return a/np.pi*.5*gamma/((x-x0)**2+(.5*gamma)**2)
def get_delay(arr1,arr2):            
          min_length = np.min(np.array([len(arr1),len(arr2)]))
          print "min_length = %s" %min_length
          a = arr1[:min_length]
          b = arr2[:min_length]                                   
          A = fftpack.fft(a)                   
          B = fftpack.fft(b)     
          Br = -B.conjugate()                                
          autocorr = np.abs(fftpack.ifft(Br*A))
          delay = np.argmax(autocorr)
          return autocorr,delay 
def dedisperse(data_arr,DM,tsamp,fs):
        #fs = np.linspace(1.28,1.53,2048)
        fs = fs/1000.
        dts = (0.00415/tsamp*DM*(fs**-2-(fs[-1])**-2)).astype(int)
        out_data_arr = np.zeros(data_arr.shape)
        for i in range(len(dts)):
                x = np.roll(data_arr[i],-1*dts[len(fs)-1 -i])
                out_data_arr[i] = x
        return out_data_arr
def disperse(data_arr,DM,tsamp,fs):
        #fs = np.linspace(1.28,1.53,2048)
        fs = fs/1000.
        dts = (0.00415/tsamp*DM*(fs**-2-(fs[-1])**-2)).astype(int)
        out_data_arr = np.zeros(data_arr.shape)
        for i in range(len(dts)):
                x = np.roll(data_arr[i],dts[len(fs)-1 -i])
                out_data_arr[i] = x
        return out_data_arr

window = 3

cands = np.genfromtxt("/home/user/cand_times_sync/heimdall_3.cand")
cands_snr = cands[:,0]
cands_dm = cands[:,5]
cands_ind = cands[:,1]
cands_width = cands[:,3]

cands2 = np.genfromtxt("/home/user/cand_times_sync_od/heimdall_2.cand")
cands2_snr = cands2[:,0]
cands2_dm = cands2[:,5]
cands2_ind = cands2[:,1]
cands2_width = cands2[:,3]


#ovro_int = 20355632#112793376
#gdscc_int = 103012504 #20355632
#ovro = FilReader("/home/user/coincidences/candidate_%s.fil" %ovro_int)
#gdscc = FilReader("/home/user/coincidences/candidate_%s.fil" %gdscc_int)
ovro = FilReader("candidate_ovro_200428.fil")
gdscc = FilReader("candidate_gdscc_200428.fil")
gdscc_arr = gdscc.readBlock(0,gdscc.header['nsamples'])
ovro_arr = ovro.readBlock(0,ovro.header['nsamples'])
ovro_ind = np.where(cands_ind == ovro_int)[0][0].astype(int)
gdscc_ind = np.where(cands2_ind == gdscc_int)[0][0].astype(int)
av_window = int(np.min([2**cands2_width[gdscc_ind],2**cands_width[ovro_ind]]))
#av_window = 26
dm = 333.3#cands_dm[ovro_ind]/2. + cands2_dm[gdscc_ind]/2.

gdscc_arr_conv = snd.uniform_filter1d(gdscc_arr,av_window,axis=1)
ovro_arr_conv = snd.uniform_filter1d(ovro_arr,av_window,axis=1)
fs = np.linspace(1280,1530,2048)
gdscc_dedisp = dedisperse(gdscc_arr_conv,dm,0.000065536,fs)
ovro_dedisp = dedisperse(ovro_arr_conv,dm,0.000065536,fs)

#low = np.where(gdscc_dedisp.mean(axis=0) == np.max(gdscc_dedisp.mean(axis=0)))[0]-2**cands2_width[gdscc_ind]-128
#high = np.where(gdscc_dedisp.mean(axis=0) == np.max(gdscc_dedisp.mean(axis=0)))[0]+2**cands2_width[gdscc_ind]+128
#low = int(low[0])
#high = int(high[0])
low = 500
high = 1500
print low,high
gdscc_corr = gdscc_dedisp[500:].mean(axis=0)[low:high]
ovro_corr = ovro_dedisp[500:].mean(axis=0)[low:high]
autocorr,delay = get_delay(ovro_corr,gdscc_corr)
delay = 0
print (high-low)-window,(high-low+window)
plt.plot(autocorr)
plt.show()
fit_autocorr = np.roll(autocorr,(high-low)/2)[delay+(high-low)/2-window-2:delay+(high-low)/2+window+1-2]
print len(autocorr),len(fit_autocorr)
xs = np.linspace(delay-window-2,delay+window-2,2*window+1)
plt.plot(xs,fit_autocorr)
plt.show()
popt, pcov = curve_fit(lorentzian,xs,fit_autocorr,p0=[32657,delay+(high-low)/2,5.])
print popt
print pcov
plt.figure(1)
plt.plot(np.roll(autocorr,(high-low)/2))
plt.plot(xs+(high-low)/2,lorentzian(xs,*popt))
plt.figure(2)
plt.plot(ovro_dedisp[500:].mean(axis=0))
plt.plot(gdscc_dedisp[500:].mean(axis=0))
plt.vlines([low,high],plt.ylim()[0],plt.ylim()[0])
plt.show()
