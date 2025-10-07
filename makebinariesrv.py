import numpy as np
from astropy.constants import G, M_sun, au

# Git check

#--------------------------------------------------------------------------------
def project(a, ecc, inc, phi, mean, m1, m2):
    '''
    Project the orbit onto the sky and return the radial velocity component.
    Parameter:
        a    -> Semi-major axis
        ecc  -> Ecentricity
        inc  -> Inclination
        phi  -> Orientation
        mean -> Mean anomaly 
    Returns:
        sky  -> Co-ordinates of companion on the sky 
        sep  -> Separation of primary and secondary on the sky (au)
    '''
    # Start value for eccentric anomaly Ee is the mean anomaly M
    # which is random uniform 0-2pi
    Ee = mean
    old = mean
    # Eccentric anomaly solver
    while True:
        old = Ee
        Ee = Ee - (Ee-ecc*np.sin(Ee)-mean)/(1.-ecc*np.cos(Ee))
        if (old-Ee)<1.E-6:
            break

    # Calculate true anomaly
    theta = 2*np.arctan(np.sqrt((1+ecc)/(1-ecc))*np.tan(0.5*Ee))

    # Calculate true distance
    rt = a*(1-ecc*ecc)/(1+ecc*np.cos(theta))

    # Orbital velocity
    vt = np.sqrt(G.value*(m1+m2)*M_sun.value * (2/(rt*au.value) - 1/(a*au.value)))

    # Flight path angle nu
    nu = np.arccos( np.sqrt( (1+ecc*np.cos(theta))/(2-rt/a) ) )

    # Get x and y positions and velocity components
    r = [0,0,0]
    v = [0,0,0]
    r[0] = rt*np.cos(theta)
    r[1] = rt*np.sin(theta)
    v[0] = -vt*np.cos(np.pi/2 - theta + nu)
    v[1] = -vt*np.sin(np.pi/2 - theta + nu)

    # Orientation around the z axis
    rx = r[0]
    ry = r[1]
    vx = v[0]
    vy = v[1]
    r[0] = rx*np.cos(phi) - ry*np.sin(phi)
    r[1] = rx*np.sin(phi) + ry*np.cos(phi)
    v[0] = vx*np.cos(phi) - vy*np.sin(phi)
    v[1] = vx*np.sin(phi) + vy*np.cos(phi)

    # Modify inclination
    rx = r[0]
    ry = r[1]
    vx = v[0]
    vy = v[1]
    r[0] = rx
    r[1] = ry*np.cos(inc)
    r[2] = ry*np.sin(inc)
    v[0] = vx
    v[1] = vy*np.cos(inc)
    v[2] = vy*np.sin(inc)
    
    # Positions on the sky
    # If these were visual binaries, with a known separation, we would return:
    # sep = np.sqrt(r[0]**2+r[2]**2)/dist
    
    # Return radial velocity (component in y direction)
    return(v[1])