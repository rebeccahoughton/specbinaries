import numpy as np

#--------------------------------------------------------------------------------
#Sky projection function
def project(a, ecc, inc, phi, mean, m1, m2):
    '''
    Numerically solve for the eccentric anomaly using the Newton-Raphson method.
    Variables:
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
    G = 6.67E-11
    msun = 1.989e30
    au = 1.5496e11
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
    #theta = mean

    # Calculate true distance
    rt = a*(1-ecc*ecc)/(1+ecc*np.cos(theta))

    # Orbital velocity
    vt = np.sqrt(G*(m1+m2)*msun * (2/(rt*au) - 1/(a*au)))
    # print(vt)

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
    #sep = np.sqrt(r[0]**2+r[2]**2)/dist
    
    # Return radial velocity (component in y direction)
    return(v[1])