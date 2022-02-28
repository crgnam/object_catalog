"""
"""

from platform import node
import spiceypy as spice

import numpy as np
import pickle

class Catalog:
    def __init__(self, path_to_data = None):
        # Define physical constants
        AU = 149597870.6907 #(km) Astronomical Unit
        deg2rad = np.pi/180.0

        # If a previous catalog is given, load that:
        if path_to_data is not None:
            with open(path_to_data, 'rb') as handle:
                data = pickle.load(handle)

        # Otherwise, obtain data from HORIZONS Database:
        else:
            #TODO: Actually make this download from HORIZONS:
            self.mu = 132712440041.93938 #(km^3/s^2)
            self.e = 0.07850100198908602
            self.a = 2.766043062222408*AU #(km)
            self.i = 10.58769305845201*deg2rad #(rad)
            self.node = 80.26859547732911*deg2rad #(rad)
            self.peri = 73.63703979153577*deg2rad #(rad)
            self.M = 291.3755993017663*deg2rad #(rad)
            self.epoch = spice.str2et('jd {}'.format(2459600.5)) #(sec) Ephemeris Time (Barycentric Dynamical Time)
            return


    def save(self, output_file = 'data.pickle'):
        with open(output_file, 'wb') as handle:
            pickle.dump(self.data, handle, protocol=pickle.HIGHEST_PROTOCOL)


    def get_states(self, et=None):
        # Calculate the time since epoch for each object:
        if et is None:
            et = self.epoch
        tsince_epoch = et - self.epoch
        # print(tsince_epoch)

        # Convert to orbital elements:
        r_system = Catalog.kepler_to_states(self.mu, self.a, self.e, self.i, self.peri, self.node, self.M, tsince_epoch)

        # Apply parent system barycenter position:
        r_inertial = r_system

        return r_inertial


    @staticmethod
    def get_visible(data):
        return visible_dict


    @staticmethod
    def kepler_to_states(mu, a, e, i, peri, node, M0, tsince_epoch):
        # Mean Anomaly:
        n = np.sqrt(mu/(a*a*a))
        L = tsince_epoch.size
        theta = np.zeros((L,1))
        for idx, t in enumerate(tsince_epoch):
            M = M0 + n*t

            # Calculate eccentric anomaly:
            KeplerConverged = 0
            convergencePercentage = 0.05
            E = 1
            for iters in range(0,1000):
                # Set up Newton-Raphson method
                f = E - e*np.sin(E) - M
                df = 1 - e*np.cos(E)
                E_new = E - f / df

                # Check for convergence
                relativeDifference = abs(E_new - E) / E * 100
                if relativeDifference < convergencePercentage:
                    break
                E = E_new

            # Calculate the true anomaly:
            theta[idx] = np.arctan2(np.sqrt(1- e*e)*np.sin(E), np.cos(E) - e)

        # Calculate the orbital momentum and radial position:
        h = np.sqrt(mu*a*(1-e*e))
        r_mag = ((h*h)/mu)*(1/(1+(e*np.cos(theta))))
        
        # Calculate the states in peri-focal coordinates:
        R_pqw = np.zeros((L,3))
        R_pqw[:,0] = np.multiply(r_mag,np.cos(theta)).flatten()
        R_pqw[:,1] = np.multiply(r_mag,np.sin(theta)).flatten()

        # V_pqw = np.zeros((L,3))
        # V_pqw[:,0] = (mu/h)*-np.sin(theta)
        # V_pqw[:,1] = (mu/h)* (e+np.cos(theta))

        # Calculate vectorized rotations:
        a11 = np.cos(node)*np.cos(peri) -np.sin(node)*np.sin(peri)*np.cos(i)
        a12 = np.sin(node)*np.cos(peri) + np.cos(node)*np.sin(peri)*np.cos(i)
        a13 = np.sin(peri)*np.sin(i)
        a21 = -np.cos(node)*np.sin(peri) - np.sin(node)*np.cos(peri)*np.cos(i)
        a22 = -np.sin(node)*np.sin(peri) + np.cos(node)*np.cos(peri)*np.cos(i)
        a23 =  np.cos(peri) * np.sin(i)
        a31 =  np.sin(node)*np.sin(i)
        a32 = -np.cos(node)*np.sin(i)
        a33 =  np.cos(i)

        # Apply rotations to obtain position in inertial frame:
        r = np.vstack( ( (a11*R_pqw[:,0] + a12*R_pqw[:,1]),
                         (a21*R_pqw[:,0] + a22*R_pqw[:,1]),
                         (a31*R_pqw[:,0] + a32*R_pqw[:,1]) ) )
        # v = [(a1.*V_pqw[:,0])+(a2.*V_pqw[:,1]),(a4.*V_pqw[:,0]+a5.*V_pqw[:,1]),(a7.*V_pqw[:,0] + a8.*V_pqw[:,1])]
        return r
            