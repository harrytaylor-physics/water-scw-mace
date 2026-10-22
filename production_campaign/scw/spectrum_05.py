import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M=np.load("md/bulk_05fs/dipoles.npy"); dt=0.5e-15; T=2250.0; M=M-M.mean(0)
def ac(x):
    n=len(x); f=np.fft.rfft(x,n=2*n); a=np.fft.irfft(f*np.conj(f),n=2*n)[:n].real; a/=(n-np.arange(n)); return a
C=sum(ac(M[:,k]) for k in range(3))
print("C(0.5fs)/C(0) = %.4f   (~0.997+ => smooth; ~0.98 => white model-noise floor)"%(C[1]/C[0]))
maxlag=min(len(C),8000); Cw=C[:maxlag]*np.hanning(2*maxlag)[maxlag:]
I=np.fft.rfft(Cw,n=2**17).real
freq=np.fft.rfftfreq(2**17,dt); wn=freq/2.99792458e10
fl=np.median(np.abs(I[(wn>3900)&(wn<5500)]))   # floor where no water modes exist
pk=np.abs(I[(wn>5)&(wn<5500)]).max()
print("white floor = %.2e  (%.1f%% of librational peak)"%(fl,100*fl/pk))
fig,ax=plt.subplots(2,1,figsize=(9,7))
m=(wn>5)&(wn<4500)
ax[0].semilogy(wn[m],np.abs(I[m])); ax[0].axhline(fl,ls=":",c="r"); ax[0].set_ylabel("|I| log"); ax[0].set_title("full density")
mz=(wn>900)&(wn<4200); Iz=np.clip(np.abs(I[mz])-fl,0,None)
ax[1].plot(wn[mz], Iz/(Iz.max() if Iz.max()>0 else 1))
ax[1].set_xlabel("cm-1"); ax[1].set_ylabel("density-floor (norm to THIS window)"); ax[1].set_title("vibrational window, floor-subtracted")
for a_ in ax:
    for x in (600,1600,3400): a_.axvline(x,ls="--",c="grey",alpha=.5)
plt.tight_layout(); plt.savefig("md/spectrum_05fs.png",dpi=120)
for lo,hi,nm in [(1200,2000,"bend"),(2800,3800,"stretch")]:
    mm=(wn>lo)&(wn<hi); print("%s: peak %.0f cm-1, %.2fx above floor"%(nm,wn[mm][np.argmax(np.abs(I[mm]))],np.abs(I[mm]).max()/fl))
print("wrote md/spectrum_05fs.png")
