import glob
import scipy.io as ff
import pandas as pd
import yt
import numpy as np
import os as os
from tqdm import tqdm

level_max = 23
cic_levelmax = 17


#====================

props=['M','x','y','z','vx','vy','vz','jx_g','jy_g','jz_g','dMB','dME','dM','rho','cs','dv','Esave','jx_bh','jy_bh','jz_bh','spinmag','eps_sink', 'rho_stars', 'rho_dm',  'vx_stars', 'vy_stars', 'vz_stars', 'vx_dm', 'vy_dm', 'vz_dm', 'n_stars', 'n_dm', 'rho_lowspeed_stars', 'rho_lowspeed_dm', 'fact_fast_stars', 'fact_fast_dm', 'sigma2_stars', 'sigma2_dm','rho_stars_far', 'dMtde','gamma','ru','rhou','rc','mrc','rhorc','trel','tdee','tdef','tde']
df={tmpprop:[] for tmpprop in props}
df=pd.DataFrame(data=df)

if not os.path.exists('./sinks'):
    os.system('mkdir sinks')
files=glob.glob('./sink_*dat')
if len(files) > 0:
    os.system('mv sink_*dat ./sinks')

yt.funcs.mylog.setLevel(40)
files=glob.glob('output_*/info*')
ds=yt.load(files[-1])
if ds.cosmological_simulation==1:
    Lbox = float(ds.length_unit.to('pccm'))
else:
    Lbox = float(ds.length_unit.to('pc'))

df=pd.concat([pd.DataFrame(data={'a':[]}),df],axis=1)
df=pd.concat([pd.DataFrame(data={'t':[]}),df],axis=1)
df=pd.concat([pd.DataFrame(data={'ID':[]}),df],axis=1)

files=glob.glob('./sinks/sink*dat')
files.sort()
dx=Lbox/2**level_max
dx_dm=Lbox/2**cic_levelmax
print('dx is {} pc for stars and {} pc for DM'.format(dx,dx_dm))
sink_props = pd.DataFrame({'t':[],'nsink':[]}) 

data = list(np.zeros(len(files)))

for i_file,f in tqdm(enumerate(files)):
    size = sink_props.index.size

    p=ff.FortranFile(f)
    sink_props.loc[size, 'nsink'] = p.read_ints()
    p.read_ints()
    a=list(p.read_reals('d'))
    sink_props.loc[size,'t'] = a[0]
    scale_l=p.read_reals('d')
    scale_d=p.read_reals('d')
    scale_t=p.read_reals('d')
    bhid=p.read_ints()
    d={tmpprop:p.read_reals('d') for tmpprop in props[:30]}
    d=pd.DataFrame(data=d, index=bhid)
    d = pd.concat([d, pd.DataFrame(data={tmpprop:p.read_ints() for tmpprop in props[30:32]}, index=bhid)], axis=1)
    d = pd.concat([d, pd.DataFrame(data={tmpprop:p.read_reals('d') for tmpprop in props[32:]}, index=bhid)], axis=1)
    t=list(p.read_reals('d'))

    d['M']*=scale_d*scale_l**3/2e33
    d['vx']*=scale_l/1e5/scale_t
    d['vy']*=scale_l/1e5/scale_t
    d['vz']*=scale_l/1e5/scale_t
    d['dMB']*=scale_d*scale_l**3/2e33 /scale_t * 3600*24*365
    d['dME']*=scale_d*scale_l**3/2e33 /scale_t * 3600*24*365
    d['dM']*=scale_d*scale_l**3/2e33
    d['dMtde']*=scale_d*scale_l**3/2e33
    d['rho']*=scale_d/1.67e-24
    d['cs']*=scale_l/1e5/scale_t
    d['dv']*=scale_l/1e5/scale_t
    d['Esave']*=scale_l/1e5/scale_t
    d['vx_stars']*=scale_l/1e5/scale_t
    d['vy_stars']*=scale_l/1e5/scale_t
    d['vz_stars']*=scale_l/1e5/scale_t
    d['vx_dm']*=scale_l/1e5/scale_t
    d['vy_dm']*=scale_l/1e5/scale_t
    d['vz_dm']*=scale_l/1e5/scale_t
    d['rho_stars']*=scale_d/1.67e-24
    d['rho_dm']*=scale_d/1.67e-24
    d['rho_lowspeed_stars']*=scale_d/1.67e-24
    d['rho_lowspeed_dm']*=scale_d/1.67e-24
    d['fact_fast_stars']*=scale_d/1.67e-24
    d['fact_fast_dm']*=scale_d/1.67e-24
    d['sigma2_stars']*=(scale_l/1e5/scale_t)**2
    d['sigma2_dm']*=(scale_l/1e5/scale_t)**2
    d['rho_stars_far']*=scale_d/1.67e-24
    d['ru'] *= scale_l/3.08e18
    d['rhou'] *= scale_d/1.67e-24
    d['rc'] *= scale_l /3.08e18
    d['mrc'] *= scale_l**3*scale_d/2e33
    d['rhorc'] *= scale_d/1.67e-24
    d['trel']  *= scale_t/3600/24/365
    d['tdee']  /= scale_t/3600/24/365
    d['tdef']  /=scale_t/3600/24/365
    d['tde']  *= scale_d*scale_l**3/2e33 /scale_t*3600*24*365

    if ds.cosmological_simulation==1:
        d['a'] = a[0]
        d['t'] = ds.cosmology.t_from_z(1/a[0]-1).in_units('Gyr') 
    else:
        d['t'] = t*scale_t/(1e9*365*24*3600)
        d['a'] = 1
    d['ID'] = d.index
    data[i_file] = d
    
df = pd.concat(data, ignore_index=True)

df['x']/=ds['boxlen']
df['y']/=ds['boxlen']
df['z']/=ds['boxlen']

df['vsink_rel_stars'] = np.sqrt((df['vx_stars']-df['vx'])**2+(df['vy_stars']-df['vy'])**2+(df['vz_stars']-df['vz'])**2)
df['vsink_rel_dm'] = np.sqrt((df['vx_dm']-df['vx'])**2+(df['vy_dm']-df['vy'])**2+(df['vz_dm']-df['vz'])**2)
df['rinf_stars'] = (df.M / 1e7) / (df.vsink_rel_stars / 200)**2
df['rinf_dm'] = (df.M / 1e7) / (df.vsink_rel_dm / 200)**2

CoulombLog = np.maximum(np.zeros(len(df.t)), np.log(4*dx/df.rinf_stars))
df['a_stars_slow']=4*np.pi*(6.67e-8)**2*df.M*2e33*df.rho_lowspeed_stars*1.67e-24*CoulombLog/(df.vsink_rel_stars*1e5)**2*3600*24*365*1e6/1e5
CoulombLog = np.minimum(np.zeros(len(df.t)), df.rinf_stars-4*dx) / (df.rinf_stars - 4*dx) 
df['a_stars_fast']=4*np.pi*(6.67e-8)**2*df.M*2e33*df.fact_fast_stars*1.67e-24*CoulombLog/(df.vsink_rel_stars*1e5)**2*3600*24*365*1e6/1e5

CoulombLog = np.maximum(np.zeros(len(df.t)), np.log(4*dx_dm/df.rinf_dm))
df['a_dm_slow']=4*np.pi*(6.67e-8)**2*df.M*2e33*df.rho_lowspeed_dm*1.67e-24*CoulombLog/(df.vsink_rel_dm*1e5)**2*3600*24*365*1e6/1e5
CoulombLog = np.minimum(np.zeros(len(df.t)), df.rinf_dm-4*dx) / (df.rinf_dm - 4*dx) 
df['a_dm_fast']=4*np.pi*(6.67e-8)**2*df.M*2e33*df.fact_fast_dm*1.67e-24*CoulombLog/(df.vsink_rel_dm*1e5)**2*3600*24*365*1e6/1e5

M=df.dv / df.cs
df['rinf_gas'] = (df.M / 1e7) / (df.dv**2 + df.cs**2)/200**2
CoulombLog = np.minimum(np.zeros(len(df.t)), df.rinf_gas-4*dx) / (df.rinf_gas - 4*dx) 
fudge=M
fudge[M < 0.95] = 1/M**2*(0.5*np.log((1+M)/(1-M)) - M)
fudge[(M >= 0.95) & (M <= 1.007)] = 1
fudge[M > 1.007] = 1/M**2*(0.5*np.log(M**2-1) + 3.2)
df['a_gas']=4*np.pi*(6.67e-8)**2*df.M*2e33*df.rho*1.67e-24/(df.cs*1e5)**2*fudge*(3600*24*365*1e6)/1e5*CoulombLog


df.sort_values(by=['ID','t'], inplace=True)
df.set_index('ID', drop=False, inplace=True)

for iBH in tqdm(df.ID.unique()):
    df.loc[iBH, 'dM']=df.loc[iBH].dM/np.gradient(df.loc[iBH].t)/1e9
    df.loc[iBH, 'dMtde']=df.loc[iBH].dMtde/np.gradient(df.loc[iBH].t)/1e9
    tmp = df.loc[iBH]
    tmp.drop(columns='ID', inplace=True)
    tmp.to_csv('./sinks/BH{:05}'.format(iBH)+'.csv', index=False)


tmp = df.loc[df.ID == iBH+1]
tmp.drop(columns='ID', inplace=True)
tmp.to_csv('./sinks/BH00000.csv', index=False)

sink_props.to_hdf('./sinks/sink_props.hdf','hdf')
df.set_index(['ID','t'], drop=False, inplace=True)
df.to_hdf('./sinks/all_sinks.hdf','hdf')
print('done')
