import os
import numpy as np
import pandas as pd
import astropy
from astropy.time import Time
from astropy import units as u
from astroquery.jplhorizons import Horizons
from shutil import copy
from .raiden import Orbits
from .wrapper import FINDORB
backend=FINDORB()

def sear(di):
    li=os.listdir(di)
    for i in li:
        if 'run' in i:
            return i
    return 'no_dir'
def ch(di):
    li=os.listdir(di)
    for i in li:
        if '.' not in i:
            return f"'{i}'"
        
    return 'no_dir'

def bashgen(out_dir, error, t1, obs_code):
    # t1 needs to be in UTC jd time format
    od = os.path.join(out_dir,'masterRunTEST.sh')
    val=Time.now().utc.value.isoformat(timespec='seconds')
    copy(os.path.dirname(os.path.abspath(__file__))+'/eph2ades.py',out_dir+'/eph2ades.py')
    with open(od,'wt') as f:
        f.write('#!/bin/sh'+"\n")
        ed1=os.path.join(out_dir,f'ephemeris/{obs_code}/')
        f.write(f'cd ./ephemeris/{obs_code}/'+'\n')
        f.write('bash '+f"'{sear(ed1)}'"+'\n'+'cd ..'+'\n'+'cd ..'+'\n')
        psvp=f'gen_{val}_observations.psv' # will integrate this into eph2ades.py
        f.write('python eph2ades.py '
                +os.path.join(f'./ephemeris/{obs_code}/',ch(ed1),'ephemeris.txt')+
                ' '+os.path.join('./orbit_determination',psvp)+' '
               +f'--astrometric_error={error} --time_error={0.01} --observatory_code={obs_code}'+'\n')
        f.write('cd orbit_determination/'+'\n')
        ed=os.path.join(out_dir,'orbit_determination/')
        f.write(f'fo "{psvp}" -O "{ch(ed1)}" -tEjd{t1} -j -D "environ.dat"'+'\n'+'cd ..'+'\n') # possibly hardcoded with -tEjd time
        f.write('cd propagation/'+'\n')
        ed=os.path.join(out_dir,'propagation/')
        f.write('bash '+f"'{sear(ed)}'")

def runFO(orbit, observatory_code, dts, astrometric_error=None, backend=FINDORB(),out=None):
    t0 = orbit.epochs[0]
    if isinstance(dts, astropy.time.core.Time):
        #assert t0 == dts[0], "The first time in the dts should be the observed time of your initial orbits"
        observation_times = dts
    elif isinstance(dts, (np.ndarray,list)):
        observation_times = t0 + dts
    if astrometric_error is None:
        astrometric_error = 0
    if out is None:
        outd=None
    else:
        whq1=orbit.ids[0].split(' ')
        nam='_'.join(whq1)
        outd=os.path.join(out,f"{nam}/{dts.max():0.0f}_{astrometric_error:0.0f}mas_{Time.now().utc.value.isoformat(timespec='seconds')}")
    MAS_TO_DEG = 2.777777777777778e-07
    DEG_TO_MAS = 1/MAS_TO_DEG
    obs = {observatory_code:observation_times}
    ephemeris,ret1= backend._generateEphemeris(orbits=orbit,observers =obs, out_dir=outd)

    ephemeris["RA_sigma_deg"] = astrometric_error*MAS_TO_DEG
    ephemeris["Dec_sigma_deg"] = astrometric_error*MAS_TO_DEG
    ephemeris["mjd_sigma_seconds"] = [0.01 for i in range(len(ephemeris))]
    ephemeris["obs_id"] = [f"o{i:07d}" for i in range(len(ephemeris))]
    #ephemeris.rename(columns={"mjd_utc" : "mjd", "orbit_id" : "obj_id"}, inplace=True)

    # Add a simple astrometric errors to both RA and Dec
    ephemeris.loc[:, "RA_deg"] += np.random.normal(loc=0, scale=astrometric_error*MAS_TO_DEG, size=len(ephemeris)) / np.cos(np.radians(ephemeris['Dec_deg'].values))
    ephemeris.loc[:, "Dec_deg"] += np.random.normal(loc=0, scale=astrometric_error*MAS_TO_DEG, size=len(ephemeris))

    od_orbit_df, residuals,ret2 = backend._orbitDetermination(ephemeris, out_dir=outd)
 
    od_orbit = Orbits(od_orbit_df)
    
    prop_orbit,ret3 = backend._propagateOrbits(orbit, observation_times[-1:], out_dir=outd)

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
    result["rms delta ra [arcsec]"] = np.sqrt(np.mean(residuals["dRA"].values**2))
    result["rms delta dec [arcsec]"] = np.sqrt(np.mean(residuals["dDec"].values**2))
    result["rms delta time [seconds]"] = np.sqrt(np.mean(residuals["dTime"].values**2))
   
    result["covariance"] = od_orbit_df["covariance"].values
    
    arc_length = observation_times.utc.mjd.max() - observation_times.utc.mjd.min()
    
    result.insert(1, "observatory_code",  observatory_code)
    result.insert(2, "arc_length [days]", arc_length)
    result.insert(5, "astrometric_error [mas]", astrometric_error)
    result.insert(3, "num_obs", len(ephemeris))
    if out is not None:
        bashgen(outd,astrometric_error,observation_times[-1:].utc.jd,observatory_code)
    
    return result

def testFO(orbits, observatory_code, dts, astrometric_error=None, out=None):
    '''
    Runs the end-to-end Find_Orb test for a list of orbits. Generates ephemeris, 
    conducts orbit determination, and propagates the initial orbits to the final 
    time and returns the comparison.

    Parameters
    ----------
    orbits : orbit object `~validate_findorb.raiden.Orbits`
        Orbits to be tested.
    observatory_code : str or list of strings
        MPC Observatory code(s) for the observatory to be tested. (500 for Geocenter)
    dts : array of floats, 2D array of floats, or list of astropy.time.core.Time objects with the same length as `orbits`
        List of observation times after the initial time to test the orbits over. Measured in days.
        NOTE: Anything passed to this parameter must have values in ASCENDING ORDER.
    astrometric_error : float or list of floats, optional
        Astrometric error to be added to the generated observations. 
        If None, no astrometric error is added. Units are milliarcseconds.
    out : str, optional
        Path to the output directory for saving necessary files for this test to be
        run independently. This includes configuration files, generated files, and bash scripts. 
        It has the file structure '{out}/{orbit_id}/{days propagated}days_{error}mas_{timestamp}'. 
        If None, no files are saved.
    
    Returns
    -------
    result : pandas.DataFrame
        DataFrame with the results of the test.
        Has the following columns:
        orbit_id : str
            Orbit ID of the object being tested.
        observatory_code : str
            MPC Observatory code for the observatory to be tested. (500 for Geocenter)
        arc_length [days] : float
            Length of time in days over which the orbits are tested.
        num_obs : int
            Number of observations generated in the test.
        num_obs_fit : int
            Number of observations used from the ephemeris file in the orbit determination.
        epoch [mjd] : float
            Time at the end of the arc in mjd (t0 + final dt day).
        astrometric_error [mas] : float
            Astrometric error added to the observations.
        delta epoch [mjd] : float
            Difference between the final epochs of the orbit determination and the propagated orbit.
        delta r [km] : float
            The absolute distance between the orbit determination result and the propagated orbit.
            Measured in kilometers.
        delta v [m/s] : float
            The absolute velocity difference between the orbit determination result and the propagated orbit.
            Measured in meters per second.
        delta x [km] : float
            The x position difference between the orbit determination result and the propagated orbit.
            Measured in kilometers.
        delta y [km] : float
            The y position difference between the orbit determination result and the propagated orbit.
            Measured in kilometers.
        delta z [km] : float
            The z position difference between the orbit determination result and the propagated orbit.
            Measured in kilometers.
        delta vx [m/s] : float
            The x velocity difference between the orbit determination result and the propagated orbit.
            Measured in meters per second.
        delta vy [m/s] : float
            The y velocity difference between the orbit determination result and the propagated orbit.
            Measured in meters per second.
        delta vz [m/s] : float
            The z velocity difference between the orbit determination result and the propagated orbit.
            Measured in meters per second.
        rms delta ra [arcsec] : float
            The root mean square difference of the Right Acension of the Residuals from Orbit Determination.
            Measured in arcseconds.
        rms delta dec [arcsec] : float
            The root mean square difference of the Declination of the Residuals from Orbit Determination.
            Measured in arcseconds.
        rms delta time [seconds] : float
            The root mean square difference of the Time of the Residuals from Orbit Determination.
            Measured in seconds.
        covariance : 2D array converted to 3D array
            Covariance matrix of the orbit determination result.
    '''
    if isinstance(observatory_code,(list,np.ndarray)) and isinstance(observatory_code[0],str):
        pass
    elif isinstance(observatory_code,str):
        observatory_code = [observatory_code]
    
    try:
        astrometric_error[0]
        errors = astrometric_error
    except:
        errors = [astrometric_error]
    
    #1-D cases
    #array of ints np.ndarray: [1,2,3,...]
    if isinstance(dts, (list,np.ndarray)) and isinstance(dts[0], (int,float,np.int64,np.float64)):
        dts = [dts]
    #list in time object: <Time object: [t1,t2,...]>
    elif isinstance(dts, astropy.time.core.Time) and isinstance(dts[0], astropy.time.core.Time):
        dts = [dts]
    #2-D cases
    #2d array of ints: [[1,2,3,...],[...]]
    elif isinstance(dts, (list,np.ndarray)) and isinstance(dts[0], (list,np.ndarray)):
        pass
    #2d array of time objects: list([<Time object: [t1,t2,...]>,<Time object: [t1,t2,...]>,...])
    elif isinstance(dts, (list,np.ndarray)) and isinstance(dts[0], astropy.time.core.Time) and len(dts) == orbits.num_orbits:
        
        results_i = []
        for i in range(orbits.num_orbits):
        #this will iterate over unique dts for each orbit
            for k in range(len(errors)):
                for l in range(len(observatory_code)):
                    result = runFO(orbits[i],observatory_code[l],dts[i],astrometric_error=errors[k],out=out)
                    results_i.append(result)
        results = pd.concat(
            results_i,
            ignore_index=True
        )
        return results
    else:
        raise ValueError('dts must be a list of floats, array of time objects, 2D array of floats, or 2D array of time objects corresponding to each orbit being tested.')
    #try:
        # check for 2d array
      #  dts[0][0]
    #except:
    #    dts = [dts]
    results_i = []
    for i in range(orbits.num_orbits):
        for j in range(len(dts)):
            for k in range(len(errors)):
                for l in range(len(observatory_code)):
                    result = runFO(orbits[i],observatory_code[l],dts[j],astrometric_error=errors[k],out=out)
                    results_i.append(result)
    results = pd.concat(
        results_i,
        ignore_index=True
    )
    return results

def loadOrb(data):
    '''
    Load orbit data from a .csv file or pandas DataFrame.

    Parameters
    ----------
    data: str or pandas.DataFrame
        Path to the .csv file or pandas DataFrame containing the orbit data.
        Must have the following columns:
        x: float
            x element of the state vector in Astronomical Units.
        y: float
            y element of the state vector in Astronomical Units.
        z: float
            z element of the state vector in Astronomical Units.
        vx: float
            x velocity element of the state vector in Astronomical Units per day.
        vy: float
            y velocity element of the state vector in Astronomical Units per day.
        vz: float
            z velocity element of the state vector in Astronomical Units per day.
        epoch or mjd_tdb: float, also can be a astropy.time.core.Time object converted to float
            Time of the state vector in mjd with a tdb scale.
    
    Returns
    -------
    orbits : orbit object `~validate_findorb.raiden.Orbits`
    '''
    return Orbits(data)

def getOrbHorizons(target, t0):
    '''
    Gets the orbital state vector from JPL Horizons for a given target and time.

    Parameters
    ----------
    target : str or list of str
        Name(s) of the target to get the orbital state vector for.
    t0 : astropy.time.core.Time
        Time object with scale='tdb' format='mjd' for the time of the state vector.

    Returns
    -------
    orbits: orbit object `~validate_findorb.raiden.Orbits`
    '''
    
    if type(target) == str:
        target = [target]
    targets_i = []
    for i in target:
        hobj = Horizons(id=i,epochs=t0.tdb.mjd,location='@sun').vectors(refplane="ecliptic",aberrations="geometric",)
        targets_i.append(hobj.to_pandas())
    targets = pd.concat(targets_i, ignore_index=True)
    return Orbits(targets,ids=targets['targetname'].values,epochs=t0+np.zeros(len(target)))
