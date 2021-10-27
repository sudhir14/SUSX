##### STEPS 1. ESCAN 2 PYHTON
#### ms, tb vs mh  scatter
#### ms, tb vs Sh scatter
#### other input, output scatters.
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import binned_statistic_2d
from scipy.interpolate import make_interp_spline, BSpline
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
fig1, (ax1, ax2, ax3)         =plt.subplots(1,3,sharey=True,constrained_layout=True,figsize=(8,2))
fig2, (ax4, ax5)              =plt.subplots(1,2,sharey=True,constrained_layout=True,figsize=(7,2))
fig3, (ax6, ax7, ax8)         =plt.subplots(1,3,sharey=True,constrained_layout=True,figsize=(8,2))
fig4, (ax9, ax10, ax11)       =plt.subplots(1,3,sharey=True,constrained_layout=True,figsize=(8,2))
fig5, (ax19)                  =plt.subplots(1,1,sharey=True,constrained_layout=True,figsize=(8,2))

n_bins=200; nmode=13; nn=1e5; tev=1.0e3
c1_m = 103.5; mh_m = 114.4; n1_m =  45.6;
drho_m = 0.0008; drho_s = 0.0017; amu_m = 2.68e-9; amu_s = 0.43e-9;
bsg_m = 3.55e-4; bsg_s = 0.24e-4; b2mu_m = 3e-9; b2mu_s =0.6e-9;
omg_m = 0.1186; omg_s = 0.002;
nsig=2.5;
mu, m1, m2, tb, mx, c1, c2, gl, haa, hbb, hcc, hgg, hll, hmm, hss, htt, hww, hza, hzz, m3,mh, ms, n1, n2, n3, n4, wh,drho,amu,b2mu, omg,psi,nsi,psd,nsd = np.loadtxt('RandomData_micro.txt', unpack=True)
sh= haa* hbb* hcc* hgg* hll* hmm* hss*  hww* hza* hzz *hmm *hcc *hss
sh=np.sqrt(pow(2.0*nn*np.pi,nmode-1)*np.exp(nmode-1)*sh)
#sh=np.log(sh)
#sh = np.exp(-sh)
sh=sh/sh.max()
#sc1 = (np.abs(n1) >  n1_m) &(np.abs(c1) > c1_m) & (np.abs(mh) > mh_m)
#sc2 = (np.abs(n1) >  n1_m) &(np.abs(c1) > c1_m) & (np.abs(mh) > mh_m) & (np.abs(amu - amu_m)/amu_s < nsig)& (np.abs(drho - drho_m)/drho_s < nsig)
#sc3 = (np.abs(n1) >  n1_m) &(np.abs(c1) > c1_m) & (np.abs(mh) > mh_m) & (np.abs(amu - amu_m)/amu_s < nsig)& (np.abs(drho - drho_m)/drho_s < nsig)&(np.abs(b2mu - b2mu_m)/b2mu_s < nsig)
#sc4 = (np.abs(n1) >  n1_m) &(np.abs(c1) > c1_m) & (np.abs(mh) > mh_m) & (np.abs(amu - amu_m)/amu_s < nsig)& (np.abs(drho - drho_m)/drho_s < nsig)&(np.abs(b2mu - b2mu_m)/b2mu_s < nsig)& (np.abs(omg) < omg_m)
sc5 = (n1 > n1_m) &(c1 > c1_m) &(mh > mh_m) &(omg < omg_m)
sc = sc5

btby, btbx, btbn= stats.binned_statistic(tb[sc],sh[sc],statistic='mean', bins=n_bins)
bmsy, bmsx, bmsn= stats.binned_statistic(np.log10(ms[sc]),sh[sc],statistic='mean', bins=n_bins)
bmhy, bmhx, bmhn= stats.binned_statistic(mh[sc],sh[sc],statistic='mean', bins=n_bins)
bhsy, bhsx, bhsn= stats.binned_statistic(ms[sc],mh[sc],statistic='mean', bins=n_bins)
bthy, bthx, bthn= stats.binned_statistic(tb[sc],mh[sc],statistic='mean', bins=n_bins)
bn1y, bn1x, bn1n= stats.binned_statistic(n1[sc]/tev,sh[sc],statistic='mean', bins=n_bins)
bc1y, bc1x, bc1n= stats.binned_statistic(c1[sc]/tev,sh[sc],statistic='mean', bins=n_bins)
bc2y, bc2x, bc2n= stats.binned_statistic(c2[sc]/tev,sh[sc],statistic='mean', bins=n_bins)
bmuy, bmux, bmun= stats.binned_statistic(mu[sc]/tev,sh[sc],statistic='mean', bins=n_bins)
bm1y, bm1x, bm1n= stats.binned_statistic(m1[sc]/tev,sh[sc],statistic='mean', bins=n_bins)
bm2y, bm2x, bm2n= stats.binned_statistic(m2[sc]/tev,sh[sc],statistic='mean', bins=n_bins)
omgy, omgx, omgn= stats.binned_statistic(np.log10(omg[sc]), sh[sc], statistic='mean', bins=n_bins)



ax1.plot(btbx[1:], smooth(btby,1))
ax2.plot(bmsx[1:], smooth(bmsy,10))
ax3.plot(bmhx[1:], smooth(bmhy,1))
ax4.plot(bhsx[1:], smooth(bhsy,1))
ax5.plot(bthx[1:], smooth(bthy,1))
ax6.plot(bn1x[1:], smooth(bn1y,5))
ax7.plot(bc1x[1:], smooth(bc1y,5))
ax8.plot(bc2x[1:], smooth(bc2y,5))
ax9.plot(bmux[1:], smooth(bmuy,5))
ax10.plot(bm1x[1:],smooth(bm1y,5))
ax11.plot(bm2x[1:],smooth(bm2y,5))
ax19.plot(omgx[1:], smooth(omgy,5))
#ax2.legend()
#ax19 = plt.axes(xscale='log')
ax4.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax1.set_xlim(0,70)
ax2.set_xlim(4,7.92)
ax3.set_xlim(90,140)
ax4.set_xlim(1e4,1e8)
ax5.set_xlim(0,70)
ax6.set_xlim(0.1,8)
ax7.set_xlim(0.1,8)
ax8.set_xlim(0.1,8)
ax9.set_xlim(0.1,8)
ax10.set_xlim(0.1,8)
ax11.set_xlim(0.1,8)
ax19.set_xlim(-2,4)

#ax1.set_ylim(0,0.001)

ax1.minorticks_on()
ax2.minorticks_on()
ax3.minorticks_on()
ax4.minorticks_on()
ax5.minorticks_on()
ax6.minorticks_on()
ax7.minorticks_on()
ax8.minorticks_on()
ax9.minorticks_on()
ax10.minorticks_on()
ax11.minorticks_on()
ax19.minorticks_on()

ax1.set(xlabel= '$\\tan\\beta$'                     , ylabel='$S/S_{max}$')
ax2.set(xlabel='$\\log_{10}M_S$ (GeV)'              , ylabel=''           )
ax3.set(xlabel='$m_h$ (GeV)'                        , ylabel=''           )
ax4.set(xlabel='$\\log_{10}M_S$ (GeV)'              , ylabel='$m_h$ (GeV)')
ax5.set(xlabel='$\\tan\\beta$'                      , ylabel=''           )
ax6.set(xlabel= '$M_{_{\\tilde\chi^0_{1}}}$ (TeV)'  , ylabel='$S/S_{max}$')
ax7.set(xlabel= '$M_{_{\\tilde\chi^\pm_{1}}}$ (TeV)', ylabel=''           )
ax8.set(xlabel= '$M_{_{\\tilde\chi^\pm_{2}}}$ (TeV)', ylabel=''           )
ax9.set(xlabel= '$\\mu$ (TeV)'                      , ylabel='$S/S_{max}$')
ax10.set(xlabel='$m_1$ (TeV)'                       , ylabel=''           )
ax11.set(xlabel='$m_2$ (TeV)'                       , ylabel=''           )
ax19.set(xlabel= '$\\Omega$ (TeV)'                      , ylabel='S/S_{max}$')

#ax2.grid()
fig1.savefig("stbmsmh.png")
fig2.savefig("mhtbms.png")
fig3.savefig("mn1c1c2.png")
fig4.savefig("mum1m2.png")
fig5.savefig("test.png")

np.savetxt("new.txt", sh, delimiter=";",fmt='%.8f')
