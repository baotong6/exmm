import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path='/Users/baotong/data_GalDisc/data/match_e_xmm/'
tablename='e_xmmdr14s_merge_spec_starinfo.xlsx'
df=pd.read_excel(path+tablename,sheet_name='xmm_matches')
kT=np.array(df['kT'])
hamstarindex=np.array(df['hamstarindex'])
Ttype=np.array(df['Ttype'])
star=np.where(Ttype=='star')[0]
notstar=np.where(Ttype=='not')[0]
print(len(kT[star]),len(kT[notstar]))
plt.hist(kT[star],bins=np.logspace(-1,1.9,30),histtype='step',color='gray',lw=2,label='kT (stars)')
plt.hist(kT[notstar],bins=np.logspace(-1,1.9,10),histtype='step',color='green',lw=2,label='kT (not stars)')
plt.xlabel('kT (keV) in tbabs*apec',fontsize=14)
plt.ylabel('Number of Sources',fontsize=14)
plt.tick_params(labelsize=14)
plt.legend()
plt.loglog()
plt.savefig(path+'useful_figs/1apec_kT_hist.pdf',bbox_inches='tight', pad_inches=0.1)
plt.show()
