import numpy as np
import sys
import os
sys.path.insert(0, '../magcolloids')
import magcolloids as mgc
ureg = mgc.ureg
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd 
import scipy.optimize as spo

def unwrap_trj(trj,bounds):
    
    trj2 = trj.copy(deep=True)
    
    def unwrap(p):
        p.iloc[:] = np.unwrap(p,axis=0)
        return p

    for c in trj.columns:
        trj2[c] = (trj2[c] - bounds[c+"_min"].values)/(bounds[c+"_max"].values - bounds[c+"_min"].values)
        
    trj2 = (trj2*2*np.pi).groupby("id",group_keys=False).apply(unwrap)/(2*np.pi)

    for c in trj.columns:
        trj2[c] = trj2[c]*(bounds[c+"_max"].values - bounds[c+"_min"].values) + bounds[c+"_min"].values

    return trj2

def calculate_velocities(trj,bounds):
    trj_unwrap = unwrap_trj(trj[["x","y","z"]],bounds)
    diffs = trj_unwrap.groupby("id",group_keys=False).diff()
    velocities = diffs/trj[["time"]].groupby("id",group_keys=False).diff().values
    velocities["time"] = trj.time

    return velocities

def load_trj_single_freq(name,directory,freq):
    # To load a video, we use the lazy read from the magcolloids package:
    trj_read = mgc.trj_lazyread(os.path.join(directory,name+".lammpstrj"), 
                     output = ["x","y","z","mux","muy","muz","fx","fy"])

    # Within the simulation, the frequency is changed with time. 
    # To obtain the frequency as a function of time we use the following definition.
    fmax = 7.5*ureg.Hz # 7.5Hz
    dt = 15*ureg.sec # 15sec
    df = 0.125*ureg.Hz

    ratio = df/dt
    total_time = fmax/ratio

    frequency = lambda time: (ratio.to((ureg.MHz/ureg.us))*(np.ceil(time.to(ureg.us)/dt.to(ureg.us))*dt)).to(ureg.Hz)

    timestep = 1e-4*ureg.s
    steps = np.array(list(trj_read.T.keys()))
    time = pd.DataFrame(data = {"time":steps*timestep.to("sec").magnitude,
                                "steps":steps,
                                "frame":range(len(steps))})
    time = time.set_index("frame")

    # We add the frequency function to the time dataset. 
    time["frequency"] = frequency(time["time"].values*ureg.sec).magnitude

    # Then we can load from the time dataset only the slice of time with the desired frequency. 
    
    trj = trj_read[time[time.frequency==freq].index]
    bounds = trj_read.get_bounds(time[time.frequency==freq].index)
    
    return trj, bounds

def subset(times, dt, rnge):
    """ Defines a repeated range of times. It returns the same range for every value of frequency"""
    fun = lambda times: (np.mod(times[:],dt))/dt
    thresh1 = (min(rnge))/dt
    thresh2 = (max(rnge))/dt
    return (fun(times)>thresh1) & (fun(times)<thresh2)
    
def load_trj_all_freq(name,directory):
    # To load a video, we use the lazy read from the magcolloids package:
    trj_read = mgc.trj_lazyread(os.path.join(directory,name+".lammpstrj"), 
                     output = ["x","y","z","mux","muy","muz","fx","fy"])

    # Within the simulation, the frequency is changed with time. 
    # To obtain the frequency as a function of time we use the following definition.
    fmax = 7.5*ureg.Hz # 7.5Hz
    dt = 15*ureg.sec # 15sec
    df = 0.125*ureg.Hz

    ratio = df/dt
    total_time = fmax/ratio

    frequency = lambda time: (ratio.to((ureg.MHz/ureg.us))*(np.ceil(time.to(ureg.us)/dt.to(ureg.us))*dt)).to(ureg.Hz)

    timestep = 1e-4*ureg.s
    steps = np.array(list(trj_read.T.keys()))
    time = pd.DataFrame(data = {"time":steps*timestep.to("sec").magnitude,
                                "steps":steps,
                                "frame":range(len(steps))})
    time = time.set_index("frame")

    # We add the frequency function to the time dataset. 
    time["frequency"] = frequency(time["time"].values*ureg.sec).magnitude

    # Then we can load from the time dataset only the slice of time with the desired frequency. 
    
    trj = trj_read[:]
    bounds = trj_read.get_bounds()
    
    return trj, bounds, time
    
def load_trj_time_range(name,directory,t_range):
    # To load a video, we use the lazy read from the magcolloids package:
    trj_read = mgc.trj_lazyread(os.path.join(directory,name+".lammpstrj"), 
                     output = ["x","y","z","mux","muy","muz","fx","fy"])

    # Within the simulation, the frequency is changed with time. 
    # To obtain the frequency as a function of time we use the following definition.

    timestep = 1e-4*ureg.s
    steps = np.array(list(trj_read.T.keys()))
    time = pd.DataFrame(data = {"time":steps*timestep.to("sec").magnitude,
                                "steps":steps,
                                "frame":range(len(steps))})
    time = time.set_index("frame")

    time_extract = time[(time.time>t_range[0]) & (time.time<t_range[1])]
    trj = trj_read[time_extract.index]
    bounds = trj_read.get_bounds(time_extract.index)
    
    return trj, bounds
    
def make_video_slowmo(name,trj,bounds):
    fig,ax = plt.subplots(1,1,figsize=(3,3),dpi=150);
    video = mgc.animate_trj(trj, region = bounds.iloc[0].values, radius = 1.4,
                            framerate = 75, start = 1, end = 2, step=1, speedup=0.1,
                            ax = ax)

    video.save(name+"_slow.gif", writer="imagemagick")
    return "<img src=\"%s\">"%(name+"_slow.gif")

def make_video_fastmo(name,trj,bounds):
    fig,ax = plt.subplots(1,1,figsize=(3,3),dpi=150);
    video = mgc.animate_trj(trj, region = bounds.iloc[0].values, radius = 1.4,
                            framerate = 75, start = 1, end = 11, step=10, speedup=1,
                            ax = ax)

    video.save(name+"_fast.gif", writer="imagemagick")
    return "<img src=\"%s\">"%(name+"_fast.gif")
    
def dipole_dipole(r,m1,m2, scaling = None):
    
    m1units = m1.units
    m2units = m2.units
    
    r = r.to(ureg.um).magnitude
    m1 = m1.magnitude
    m2 = m2.magnitude
        
    rmag = np.linalg.norm(r)
    
    if scaling is None:
        u0 = 4e-7*np.pi*ureg.N/ureg.A**2
        scaling = 3/4/np.pi*u0

    first_term = (np.cross(np.cross(r,m1),m2) + np.cross(np.cross(r,m2),m1))/rmag**5
    second_term = -(2*np.dot(m1,m2)*r)/rmag**5
    third_term = 5*np.dot(np.cross(r,m1),np.cross(r,m2))*r/rmag**7
    
    force = scaling * (first_term+second_term+third_term) * (m1units*m2units) / ureg.um**4
    
    return force.to(ureg.pN).magnitude
    
def get_circle(trj):
    
    def calc_R(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((trj.x.values-xc)**2 + (trj.y.values-yc)**2)

    def f_2(c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()    

    center_estimate = 0, 0
    center, ir = spo.leastsq(f_2, center_estimate)
    radius = np.mean(calc_R(*center))
    
    return center, radius

def dimers_array(dim):
    dim["member_a"] = np.array([list(m) for m in dim.members])[:,0]
    dim["member_b"] = np.array([list(m) for m in dim.members])[:,1]

    dim["x"] = np.array([m for m in dim.center])[:,0]
    dim["y"] = np.array([m for m in dim.center])[:,1]
    dim["z"] = np.array([m for m in dim.center])[:,2]

    dim["dx"] = np.array([m for m in dim.direction])[:,0]
    dim["dy"] = np.array([m for m in dim.direction])[:,1]
    dim["dz"] = np.array([m for m in dim.direction])[:,2]
    return dim.filter(["member_a","member_b","x","y","z","dx","dy","dz"])
