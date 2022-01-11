import os
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy import units as u
from shutil import copy
from thor.orbits import Orbits
from thor.backend import FINDORB
backend=FINDORB()

def sear(di):
    li=os.listdir(di)
    for i in li:
        #print(i)
        if 'run' in i:
            return i
    return 'no_dir'
def ch(di):
    li=os.listdir(di)
    for i in li:
        if '.' not in i:
            return f"'{i}'"
        
    return 'no_dir'

def bashgen(out_dir):
    od = os.path.join(out_dir,'masterRunTEST.sh')
    val=Time.now().utc.value.isoformat(timespec='seconds')
    #print(out_dir)
    copy(os.path.dirname(os.path.abspath(__file__))+'/eph2ades.py',out_dir+'/eph2ades.py')
    with open(od,'wt') as f:
        f.write('#!/bin/sh'+"\n")
        ed1=os.path.join(out_dir,'ephemeris/500/')
        #print(ed1)
        f.write('cd ./ephemeris/500/'+'\n')
        f.write('bash '+f"'{sear(ed1)}'"+'\n'+'cd ..'+'\n'+'cd ..'+'\n')
        psvp=f'gen_{val}_observations.psv' # will integrate this into eph2ades.py
        f.write('python eph2ades.py '
                +os.path.join('./ephemeris/500/',ch(ed1),'ephemeris.txt')+
                ' '+os.path.join('./orbit_determination',psvp)+' '
               +f'--astrometric_error={0.} --time_error={0.01}'+'\n')
        f.write('cd orbit_determination/'+'\n')
        ed=os.path.join(out_dir,'orbit_determination/')
        f.write(f'fo "{psvp}" -O "{ch(ed1)}" -tEjd2459029.5008007498 -j -D "environ.dat"'+'\n'+'cd ..'+'\n') # possibly hardcoded with -tEjd time
        f.write('cd propagation/'+'\n')
        ed=os.path.join(out_dir,'propagation/')
        f.write('bash '+f"'{sear(ed)}'")

def testFO(orbit, observatory_code, t0, dts, astrometric_error=None, backend=FINDORB(),out=None):
    #np.random.seed(3)
    if out is None:
        outd=None
       # youd = None
    else:
        #youd = str(uuid.uuid4())
        #outd=os.path.join(out,youd)
        whq1=orbit.ids[0].split(' ')
        nam='_'.join(whq1)
        outd=os.path.join(out,f"{nam}/{dts.max()}days_{10.:0.0f}mas_{Time.now().utc.value.isoformat(timespec='seconds')}")
    MAS_TO_DEG = 2.777777777777778e-07
    DEG_TO_MAS = 1/MAS_TO_DEG
    observation_times = t0 + dts
    obs = {observatory_code:observation_times}
    ephemeris,ret1= backend._generateEphemeris(orbits=orbit,observers =obs, out_dir=outd)
    #writeCall(ret1,'ephemeris',outd)
    if astrometric_error is None:
        astrometric_error = 0

    ephemeris["RA_sigma_deg"] = astrometric_error*MAS_TO_DEG
    ephemeris["Dec_sigma_deg"] = astrometric_error*MAS_TO_DEG
    ephemeris["mjd_sigma_seconds"] = [0.01 for i in range(len(ephemeris))]
    ephemeris["obs_id"] = [f"o{i:07d}" for i in range(len(ephemeris))]
    #ephemeris.rename(columns={"mjd_utc" : "mjd", "orbit_id" : "obj_id"}, inplace=True)

    # Add a simple astrometric errors to both RA and Dec
    ephemeris.loc[:, "RA_deg"] += np.random.normal(loc=0, scale=astrometric_error*MAS_TO_DEG, size=len(ephemeris)) / np.cos(np.radians(ephemeris['Dec_deg'].values))
    ephemeris.loc[:, "Dec_deg"] += np.random.normal(loc=0, scale=astrometric_error*MAS_TO_DEG, size=len(ephemeris))
    #observationsToADES(ephemeris,'./',orbit.ids[0].split()[0]+'error')
    # Run OD with the desired backend
    od_orbit_df, residuals,ret2 = backend._orbitDetermination(ephemeris, out_dir=outd)
    #writeCall(ret2,'orbit_determination',outd)
    od_orbit = Orbits.from_df(od_orbit_df)
    
    prop_orbit,ret3 = backend._propagateOrbits(orbit, observation_times[-1:], out_dir=outd)
    #writeCall(ret3,'propagation',outd)
    delta_state = od_orbit_df[["x", "y", "z", "vx", "vy", "vz"]].values - prop_orbit[["x", "y", "z", "vx", "vy", "vz"]].values
    
    delta_r = np.linalg.norm(delta_state[:,0:3])
    delta_v = np.linalg.norm(delta_state[:,3:6])
    
    # os.rmtree, rmdir add parameter check and dir deletion
    
    result = pd.DataFrame()
    
    result['orbit_id'] = od_orbit_df.orbit_id.values
    result["num_obs_fit"] = residuals["incl"].sum().astype(int)
    result['epoch [mjd]'] = prop_orbit.mjd_tdb.values
    result["delta epoch [mjd]"] = od_orbit.epochs.tdb.mjd - prop_orbit["mjd_tdb"].values
    result["delta r [km]"] = (delta_r * u.AU).to(u.km).value
    result["delta v [m/s]"] = (delta_v * u.AU / u.d).to(u.m / u.s).value 
    result["delta x [km]"] = (delta_state[:,0] * u.AU).to(u.km).value
    result["delta y [km]"] = (delta_state[:,1] * u.AU).to(u.km).value
    result["delta z [km]"] = (delta_state[:,2] * u.AU).to(u.km).value
    result["delta vx [m/s]"] = (delta_state[:,3] * u.AU / u.d).to(u.m / u.s).value
    result["delta vy [m/s]"] = (delta_state[:,4] * u.AU / u.d).to(u.m / u.s).value
    result["delta vz [m/s]"] = (delta_state[:,5] * u.AU / u.d).to(u.m / u.s).value
    result["rms delta ra"] = np.sqrt(np.mean(residuals["dRA"].values**2))
    result["rms delta dec"] = np.sqrt(np.mean(residuals["dDec"].values**2))
    result["rms delta time"] = np.sqrt(np.mean(residuals["dTime"].values**2))
   
    result["covariance"] = od_orbit_df["covariance"].values
    
    arc_length = observation_times.utc.mjd.max() - observation_times.utc.mjd.min()
    
    result.insert(1, "observatory_code",  observatory_code)
    result.insert(2, "arc_length [days]", arc_length)
    result.insert(5, "astrometric_error [mas]", astrometric_error)
    result.insert(3, "num_obs", len(ephemeris))
    if out is not None:
        bashgen(outd)
    
    return result#,ret1,ret2,ret3#,residuals