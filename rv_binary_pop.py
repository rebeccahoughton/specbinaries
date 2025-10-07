import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import makebinariesrv as rvbin # This is another file in the same directory
from astropy.table import Table
from scipy.stats import ks_2samp
from sys import exit

# # Binary tree python files
# # Likely not needed here, as we're comparing 1D distributions.
# import binary_tree as btree
# import btree_plot as btplot
# from binary_tree import Region

# Plot formatting
mpl.rcParams['text.usetex'] = True
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.size'] = 14.4

# ------------------------------------------------------------------------------
class Params:
    '''
    Generates in the initial conditions for the fake binary population. 
    Variables:
        n           -> Number of test systems
        mean        -> Mean of the semi-major axis distribution
        std         -> Standard deviation of semi-major axis distribution
        ashape      -> The shape of the semi-major axis distribution. Either
                       "lognormal", "gaussian", or "flat". 
        dshape      -> The shape of the distance distribution. Either "flat"
                       or "obs" (for random sampling from observed distribution). 
        alim        -> Limits for the flat semi-major axis distribution.
        dlim        -> Limits on distances in the sample. 
    Returns:
        self.a, self.ecc, self.inc, self.phi, self.meani, self.dmag, self.dist
    '''

    def __init__(self, n, mean, std, ashape, max_ecc=None, dshape=None, alim=None):
        self.ashape = ashape
        self.dshape = dshape
        self.alim = alim
        self.max_ecc = max_ecc
        self.mean = mean
        self.std = std
        self.n = n


    def semimajoraxis(self):
        '''Select semi-major axis from eith a lognormal, gaussian or flat
           distribution.'''
        if self.ashape == "lognormal":
            self.var = self.std**2
            mu = np.log(self.mean**2/(np.sqrt(self.var+self.mean**2)))
            sigma = np.sqrt(np.log(1+(self.var/(self.mean**2))))
            self.a = np.exp(mu + sigma * np.random.normal(size=self.n))

        elif self.ashape == "gaussian":
            mu = self.mean
            sigma - self.std
            self.a = np.random.normal(loc=mu, scale=sigma, ize=self.n)

        elif self.ashape == "flat":
            self.a = np.random.uniform(self.alim[0], self.alim[1], size=self.n)
        

    def eccentricity(self):
        '''Select value for eccentricity.'''
        self.ecc = np.random.uniform(0, 0.9, self.n) if self.max_ecc==None else np.random.uniform(0, self.max_ecc, self.n)


    def masses(self):
        '''Generate the masses for the two stars.'''
        m_tot = np.random.choice(spec_data['Mass'], size=self.n)

        # Assuming a flat mass ratio distribution
        # You can customise this if you want, to test different mass ratio distributions!
        q = np.random.uniform(0.2, 1, self.n)
        self.m1 = m_tot/(1+q)
        self.m2 = q*self.m1
        

    def errors(self):
        '''Errors on the fake rv measurements drawn 
           from the real data distribution.'''
        err_rv1 = np.random.choice(rv_data['uncertRV'],size=self.n)
        err_rv2 = np.random.choice(rv_data['uncertRV'],size=self.n)


    # # Distance is relevant for visual binaries, not spectroscopic
    # def distance(self):
    #     '''Select the distance in kpc.'''
    #     if self.dshape == "flat":
    #         self.dist = np.random.uniform(self.dlim[0], self.dlim[1], size=self.n)
    #     elif self.dshape == "obs":
    #         dists = all_objects['dist']
    #         self.dist = np.abs(np.random.choice(dists, size=self.n))


    def phase_params(self):
        '''Select the inclination, orientation and mean anomaly of the system.'''
        self.inc = np.arcsin(np.random.uniform(size=self.n))
        self.phi = np.random.uniform(size=self.n)*2*np.pi
        self.meani = np.random.uniform(size=self.n)*2*np.pi


    def time_difference(self):
        '''Select the time difference between the two observations
           from the observed distribution.'''
        self.dtobs = np.random.choice(rv_data['deltaT'], size=self.n)
# --------------------------------------------------------------------------------------


def move_companion(a, ecc, inc, phi, meani, m1, m2, dtobs):
    '''Move the companion along its orbit to the time of the second observation.
       Args:
           a (float): Semi-major axis of the orbit.
           ecc (float): Eccentricity of the orbit.
           inc (float): Inclination of the orbit.
           phi (float): Longitude of the ascending node.
           meani (float): Mean anomaly at the first observation.
           m1 (float): Mass of the primary star.
           m2 (float): Mass of the secondary star.
           dtobs (float): Time difference between the two observations.

       Returns:
            rv2 (float): Radial velocity of the secondary star at the time 
                         of the second observation.
    '''
    # Equal chance of being either prograde or retrograde
    f = np.random.choice([-1,1])

    #Period (in days)
    T = np.sqrt(a**3/(m1 + m2))

    # Fraction of orbit moved
    dM = (dtobs/T)*2*np.pi
    meanf = meani + (f*dM)
    rv2 = rvbin.project(a, ecc, inc, phi, meanf, m1, m2)
    return rv2


def get_params(n, mean, std, ashape, max_ecc=None, dshape=None, alim=None):
    '''Get the orbital parameters of our sample.'''

    # Getting the orbital parameters of the system
    params = Params(n, mean, std, ashape, max_ecc, dshape, alim)
    params.semimajoraxis()   # in au
    params.eccentricity()    # dimensionless
    params.masses()          # in solar masses
    params.phase_params()    # in radians
    params.time_difference() # in days

    # Initialise empty rv arrays to store the radial velocities
    rv1 = np.zeros((n))
    rv2 = np.zeros((n))

    # Calculate radial velocities for each fake system
    # rv1 is the radial velocity at the first epoch
    # rv2 is the radial velocity after it's moved along its orbit
    # to the time of the second observation
    for i in range(n):
        rv1[i] = rvbin.project(params.a[i], 
                               params.ecc[i], 
                               params.inc[i],
                               params.phi[i], 
                               params.meani[i], 
                               params.m1[i], 
                               params.m2[i])

        # Getting the second epoch of observation
        rv2[i] = move_companion(params.a[i], 
                                params.ecc[i], 
                                params.inc[i], 
                                params.phi[i],
                                params.meani[i], 
                                params.m1[i], 
                                params.m2[i], 
                                params.dtobs[i])

    # Calculate the difference in rv and select an error
    drv = rv2 - rv1
    err_rv = drv * np.random.choice(rv_data['uncertRV']/rv_data['deltaRV'])
    # Return
    return drv, err_rv, params


def has_empty_rows(array):
    '''Check if any of the arrays have empty rows.'''
    for val in array:
        if type(val)==np.ma.core.MaskedConstant:
            raise ValueError("Array {} contains masked constants".format(array.name))
        

def selection_effects():
    '''Incorporate selection effects into the model.'''
    return

# Set random seed
# This makes the results reproducible
np.random.seed(10)

# Load in the observations using pandas, dropping all the NaNs
rvs_pd = pd.read_csv('TableRV_nonSB_SB1.csv').dropna()
MassAge_pd = pd.read_csv('TableSB_MassAge.csv').dropna()

# Convert into astropy table (because I like astropy tables)
rv_data = Table.from_pandas(rvs_pd)
spec_data = Table.from_pandas(MassAge_pd)

# Print the column names to check everything has loaded in correctly
print("Available data:")
print(spec_data.colnames)
print(rv_data.colnames)

# List to append the number of epochs each binary has
epochs = []

# List of which objects are SB1s
SB1 = []
for i in range(len(spec_data['obj'])):
    x = spec_data['obj'][i]
    binary = spec_data['SpecBin'][i]
    if binary =='YES':
        n_obs = list(sorted(rv_data['object'])).count(x)
        epochs.append(n_obs)
        SB1.append(x) if n_obs>0 else None

# Finding which of the sytems are binaries
mask = np.where((spec_data['SpecBin']=='YES') & (spec_data['SBtype']=='SB1'))[0]
binaries  = spec_data['obj'][mask]

# Calling the model
# model = Params(mean=250,
#                std=100,
#                ashape='lognormal',
#                dshape='obs',
#                theta=[250,100,10,10],
#                alim=[0,1],
#                n=1000)

drv, err_rv, params = get_params(n=1000,
                                 mean=100,
                                 std=50,
                                 ashape='lognormal',
                                 max_ecc=0.5)
                                #  ashape='flat',
                                #  alim=[2000,5000])

print(min(np.abs(rv_data['deltaRV'])))
print(np.shape(drv[drv > min(rv_data['deltaRV'])]))
# exit()

def compare_distributions(model_data, obs_data):
    '''Compare the distributions of the observed and model delta RVs.'''
    # Use a KS test to compare the two distributions
    ks_stat, p_value = ks_2samp(model_data, obs_data)
    print("KS statistic: {:.3f}, p-value: {:.3f}".format(ks_stat, p_value))
    return

# A low value of the KS statistic indicates that the two distributions are similar
compare_distributions(drv, rv_data['deltaRV'])
exit()

# Plot histograms of all the orbital parameters as a sanity check
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=[12,6])
ax[0,0].hist(params.a, bins=20)
ax[0,1].hist(params.ecc, bins=20)
ax[0,2].hist(params.inc, bins=20)
ax[1,0].hist(params.phi, bins=20)
ax[1,1].hist(params.meani, bins=20)
ax[1,2].hist(params.dtobs, bins=20)
# Labels
ax[0,0].set_xlabel("Semi-major axis (au)")
ax[0,1].set_xlabel("Eccentricity")
ax[0,2].set_xlabel("Inclination (radians)")
ax[1,0].set_xlabel("Orientation (radians)")
ax[1,1].set_xlabel("Mean anomaly (radians)")
ax[1,2].set_xlabel("Time between observations (days)")


# Now we want to compare the parameters of the fake binaries to the observed
# binaries. We can do this by plotting histograms of the semi-major axis (fake)
# and delta RV distributions (fake and observed).

# When the fake and observed delta RV distributions match, we can conclude that
# the semi-major axis distribution we assumed for the fake binaries is a good
# representation of the true semi-major axis distribution of the observed
# binaries.
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[9,4])
ax[0].hist(params.a, bins=20)
ax[0].set_xlabel("Semi-major axis (au)")

ax[1].hist(drv, 
           bins=np.arange(min(drv), max(drv) + 1, 1), 
           color="r", 
           label='model', 
           alpha=0.5)
ax[1].hist(rv_data['deltaRV'], 
           bins=np.arange(min(rv_data['deltaRV']), max(rv_data['deltaRV']) + 1, 1), 
           color="b", 
           label='observed', 
           alpha=0.3)
ax[0].tick_params(axis="both", direction="in")
ax[1].tick_params(axis="both", direction="in")
ax[0].set_ylabel("Number of binaries")
ax[1].set_xlabel("delta RV")
ax[1].set_xlim(-20,20)
ax[1].legend()
plt.show()