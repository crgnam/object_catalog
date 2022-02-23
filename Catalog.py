"""
"""

import spiceypy as spice
import numpy as np
import pickle

class Catalog:
    def __init__(self, path_to_dict = None):
        # If a previous catalog is given, load that:
        if path_to_dict is not None:
            with open(path_to_dict, 'rb') as handle:
                self.catalog_dict = pickle.load(handle)

        # Obtain data from HORIZONS Database:
        self.catalog_dict


    @classmethod
    def save(self,output_file = 'catalog_dict.pickle'):
        with open(output_file, 'wb') as handle:
            pickle.dump(self.catalog_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


    @classmethod
    def get_states(self, et=None):
        # Extract the orbital elements:
        mu = self.elements[:,0]
        a = self.elements[:,1]
        e = self.elements[:,2]
        i = self.elements[:,3]
        omega = self.elements[:,4]
        RAAN = self.elements[:,5]
        M0 = self.elements[:,6]

        # Calculate the time since epoch for each object:
        if et is None:
            et = self.epochs
        tsince_epoch = et - self.epochs

        # Convert to orbital elements:

        # Apply parent barycenter position:

        return states


    @staticmethod
    def kepler_to_states(mu, a, e, i, omega, RAAN, M0, tsince_epoch):

        # Mean Anomaly:
        n = np.sqrt(mu/(a^3))
        M = M0 + n*tsince_epoch

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
        theta = np.arctan2(np.sqrt(1- e*e)*np.sin(E), np.cos(E) - e)
        
        # Calculate the states in peri-focal coordinates:
        R_pqw(:,1) = r.*cos(theta)
        R_pqw(:,2) = r.*sin(theta)
        R_pqw(:,3) = 0
        V_pqw(:,1) = (mu./h).*-sin(theta)
        V_pqw(:,2) = (mu./h).* (e+cos(theta))
        V_pqw(:,3) = 0

        # Calculate the orbital momentum and radial position:
        h = np.sqrt(mu.*a.*(1-e.*e))
        r = ((h.*h)./mu).*(1./(1+(e.*cos(theta))))

        # Calculate vectorized rotations:
        a11 = np.cos(RAAN)*np.cos(w) -np.sin(RAAN)*np.sin(w)*np.cos(inc)
        a12 = np.sin(RAAN)*np.cos(w) + np.cos(RAAN)*np.sin(w)*np.cos(inc)
        a13 = np.sin(w)*np.sin(inc)
        a21 = -np.cos(RAAN)*np.sin(w) - np.sin(RAAN)*np.cos(w)*np.cos(inc)
        a22 = -np.sin(RAAN)*np.sin(w) + np.cos(RAAN)*np.cos(w)*np.cos(inc)
        a23 =  np.cos(w) * np.sin(inc)
        a31 =  np.sin(RAAN)*np.sin(inc)
        a32 = -np.cos(RAAN)*np.sin(inc)
        a33 =  np.cos(inc)

        # Apply rotations to obtain position in inertial frame:
        r = [(a1.*R_pqw(:,1))+(a2.*R_pqw(:,2)),(a4.*R_pqw(:,1)+a5.*R_pqw(:,2)),(a7.*R_pqw(:,1) + a8.*R_pqw(:,2))]
        v = [(a1.*V_pqw(:,1))+(a2.*V_pqw(:,2)),(a4.*V_pqw(:,1)+a5.*V_pqw(:,2)),(a7.*V_pqw(:,1) + a8.*V_pqw(:,2))]
        return r,v
            