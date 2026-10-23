import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
M=np.load("md/bulk_dense/dipoles.npy"); dt=2e-15; T=2250.0
M=M-M.mean(0)
def ac(x):
    n=len(x); f=np.fft.rfft(x,n=2*n); a=np.fft.irfft(f*np.conj(f),n=2*n)[:n].real
    a/=(n-np.arange(n)); return a
C=sum(ac(M[:,k]) for k in range(3))
maxlag=min(len(C),6000); Cw=C[:maxlag]*np.hanning(2*maxlag)[maxlag:]
spec=np.fft.rfft(Cw,n=2**16).real
freq=np.fft.rfftfreq(2**16,dt); wn=freq/2.99792458e10
hbar=1.054571817e-34; kB=1.380649e-23; w=2*np.pi*freq
alpha=w*(1-np.exp(-hbar*w/(kB*T)))*spec; alpha[0]=0
m=(wn>0)&(wn<4500)
plt.figure(figsize=(8,5)); plt.plot(wn[m],alpha[m]/alpha[m].max())
for x,l in [(600,"libration"),(1600,"bend"),(3400,"stretch")]:
    plt.axvline(x,ls="--",c="grey",alpha=.5); plt.text(x,1.02,l,rotation=90,va="bottom",fontsize=8)
plt.xlabel("wavenumber (cm$^{-1}$)"); plt.ylabel("absorption (norm.)")
plt.title("dipole-ACF spectrum  2250 K / 1256 kg m$^{-3}$"); plt.tight_layout()
plt.savefig("md/spectrum_2250.png",dpi=130)
print("wrote md/spectrum_2250.png  peak=%.0f cm-1"%wn[m][np.argmax(alpha[m])])
