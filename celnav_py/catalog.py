"""
"""
from platform import node
import spiceypy as spice
import requests
import os
from collections import namedtuple

import numpy as np
import pickle
import json

DataStruct = namedtuple("DataStruct", "mu spkid epoch a e i peri node M H G")

class Catalog:
    def __init__(self, path_to_data):
        # Define physical constants
        AU = 149597870.6907 #(km) Astronomical Unit
        deg2rad = np.pi/180.0

        # If a previous catalog is given, load that:
        if os.path.exists(path_to_data):
            print("Loading previously generated catalog...",flush=True,end='')
            with open(path_to_data, 'rb') as handle:
                data = pickle.load(handle)
            self.mu    = data.mu
            self.spkid = data.spkid
            self.epoch = data.epoch
            self.a     = data.a
            self.e     = data.e
            self.i     = data.i
            self.peri  = data.peri
            self.node  = data.node
            self.M     = data.M
            self.H     = data.H
            self.G     = data.G
            self.num_objects = self.mu.size
            print('DONE')

        # Otherwise, obtain data from HORIZONS Database:
        else:
            # Submit query:
            sbdb_api = "https://ssd-api.jpl.nasa.gov/sbdb_query.api?"
            query = "fields=spkid,epoch,a,e,i,w,om,ma,H,G&full-prec=1"
            url = sbdb_api + query
            print("Obtaining catalog using query: {} ...".format(query),flush=True,end='')
            retval = requests.get(url)
            print("DONE")
            print("Decoding json formatted response...",flush=True,end='')
            json_data = retval.text
            json_decoded = json.loads(json_data)
            data = json_decoded['data']
            print('DONE')

            # Preallocate numpy arrays:
            self.num_objects = len(data)
            self.mu = np.zeros((self.num_objects,1))
            self.spkid = np.zeros((self.num_objects,1))
            self.epoch = np.zeros((self.num_objects,1))
            self.a  = np.zeros((self.num_objects,1))
            self.e  = np.zeros((self.num_objects,1))
            self.i  = np.zeros((self.num_objects,1))
            self.peri = np.zeros((self.num_objects,1))
            self.node = np.zeros((self.num_objects,1))
            self.M = np.zeros((self.num_objects,1))
            self.H = np.zeros((self.num_objects,1))
            self.G = np.zeros((self.num_objects,1))

            # Loop through all loaded objects from catalog:
            print('Storing into class members...',flush=True,end='')
            remove_inds = []
            for idx,obj_data_str in enumerate(data):
                obj_data = np.array(obj_data_str, dtype=np.float)

                # Load in each value and convert units to radians and km
                self.mu[idx] = 132712440041.93938 #(km^3/s^2) TODO Aquire this for each object

                self.spkid[idx] = int(obj_data[0])
                self.epoch[idx] = spice.str2et('jd {}'.format(obj_data[1])) #(sec) Ephemeris Time (Barycentric Dynamical Time)

                self.a[idx] = obj_data[2]*AU 
                self.e[idx] = obj_data[3]
                self.i[idx] = obj_data[4]*deg2rad
                self.peri[idx] = obj_data[5]*deg2rad
                self.node[idx] = obj_data[6]*deg2rad
                self.M[idx] = obj_data[7]*deg2rad

                self.H[idx] = obj_data[8]
                self.G[idx] = obj_data[9]

                # Determine if it is a hyperbolic/parabolic trajectory:
                if (self.a[idx] <= 0 or self.e[idx] >= 1 or np.isnan(self.a[idx]) or np.isnan(self.e[idx]) or np.isnan(self.i[idx]) or
                    np.isnan(self.peri[idx]) or np.isnan(self.node[idx]) or np.isnan(self.M[idx]) or np.isnan(self.H[idx]) ):
                    remove_inds.append(idx)
            print('DONE')

            # Remove invalid objects
            print('Removing objects on hyperbolic/parabolic trajectory or missing elements...',flush=True,end='')
            self.mu = np.delete(self.mu, remove_inds)
            self.spkid = np.delete(self.spkid, remove_inds)
            self.epoch = np.delete(self.epoch, remove_inds)
            self.a = np.delete(self.a, remove_inds)
            self.e = np.delete(self.e, remove_inds)
            self.i = np.delete(self.i, remove_inds)
            self.peri = np.delete(self.peri, remove_inds)
            self.node = np.delete(self.node, remove_inds)
            self.M = np.delete(self.M, remove_inds)
            self.H = np.delete(self.H, remove_inds)
            self.G = np.delete(self.G, remove_inds)
            self.num_objects = self.num_objects - len(remove_inds)
            print('DONE (removed {} objects)'.format(len(remove_inds)))
            
            self.save(path_to_data)
        return


    def save(self, output_file):
        print("Saving data to .pickle file...",flush=True,end='')
        data = DataStruct(self.mu, self.spkid, self.epoch, self.a, self.e, self.i, self.peri, self.node, self.M, self.H, self.G)
        with open(output_file, 'wb') as handle:
            pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print('DONE')


    def get_states(self, et=None):
        print("Calculating position at given time...", flush=True, end='')
        # Calculate the time since epoch for each object:
        if et is None:
            tsince_epoch = np.zeros((self.num_objects, 0)) # Evaluate at epoch 
        elif type(et) is float:
            tsince_epoch = et - self.epoch
        #TODO Add error handling....

        # Convert to orbital elements:
        r_system = np.zeros((3, self.num_objects))
        r_system = Catalog.kepler_to_state(self.mu, self.a, self.e, self.i, self.peri, self.node, self.M, tsince_epoch)
        print("DONE")

        # Apply parent system barycenter position:
        print("Apply system barycenter to satellites...", flush=True, end='')
        r_inertial = r_system
        print("DONE")
        
        return r_inertial


    @staticmethod
    def write_catalog(objects, output_file):
        # Name     SPK a e i peri node m0 epoch GM J2
        # JUPITER  5   
        # Europa   502

        # Format the outputs:
        output_data = np.hstack((spk_ids, object_elements))

        # Write to a CSV file:
        np.savetxt(output_file, output_data, delimiter=",")
        return

    @staticmethod
    def get_visible(observer_position, object_positions, H, G, magnitude_maximum, sun_angle_minimum):
        # CITATION OF H-G MAGNITUDE LAW: 
        #    Karri Muinonen, Irina N. Belskaya, Alberto Cellino, Marco Delbò, Anny Chantal Levasseur-Regourd,
        #        et al.. A three-parameter magnitude phase function for asteroids. Icarus, Elsevier, 2010, 209 (2),
        #        pp.542-555. ff10.1016/j.icarus.2010.04.003ff. ffhal-00676207f
        #   url: https://hal.archives-ouvertes.fr/hal-00676207/document (pages 8 and 9 of document)

        # Pre-calculate relevant values:
        observer_position = observer_position.reshape((3,1))
        observer_to_sun_norm = np.linalg.norm(observer_position)
        observer_to_sun_unit = -observer_position/observer_to_sun_norm

        observer_to_objects = object_positions - observer_position
        observer_to_objects_norm = np.linalg.norm(observer_to_objects, axis=0)
        observer_to_objects_unit = observer_to_objects/observer_to_objects_norm

        sun_to_objects_norm = np.linalg.norm(object_positions, axis=0)
        sun_to_objects_unit = -object_positions/sun_to_objects_norm

        # Calculate objects that are too angularly close to the sun:
        sun_angle = np.sum(np.multiply(observer_to_objects_unit, observer_to_sun_unit), axis=0)

        # Calcualte the phase angle:
        alpha = np.sum(np.multiply(observer_to_objects_unit, sun_to_objects_unit), axis=0) # Phase angle

        # Calculate the apparent magnitude (using Approximated H-G Magnitude Law):
        PHI_1 = np.exp(-3.33*np.tan((1/2)*alpha)**(0.63)) # Equation 6 in referenced document
        PHI_2 = np.exp(-1.87*np.tan((1/2)*alpha)**(1.22)) # Equation 6 in referenced document
        magnitude = H - 2.5*np.log10( (1-G)*PHI_1 + G*PHI_2) # Equation 2 in referenced document

        # If any values are NaN, recalculate using simple magnitude only law:
        recalculate_indices = np.where(np.isnan(magnitude))[0].tolist()
        q_alpha = (2/3)*((1 - alpha[recalculate_indices]/np.pi)*np.cos(alpha[recalculate_indices]) + (1/np.pi)*np.sin(alpha[recalculate_indices]))
        d_bs = sun_to_objects_norm[recalculate_indices]
        d_bo = observer_to_objects_norm[recalculate_indices]
        d_o = 149597870.6907
        magnitude[recalculate_indices] = H[recalculate_indices] + 5*np.log10(d_bs*d_bo/(d_o**2)) - 2.5*np.log10(q_alpha)

        # Determine a boolean array of the potentially visible objects:
        visible_objects = (magnitude <= magnitude_maximum) & (sun_angle <= sun_angle_minimum)
        return visible_objects


    @staticmethod
    def kepler_to_state(mu, a, e, i, peri, node, M0, tsince_epoch):
        # Calculate Mean Anomaly:
        n = np.sqrt(mu/(a*a*a))
        M = M0 + n*tsince_epoch

        # Calculate eccentric anomaly:
        theta = np.zeros(M0.shape)
        for idx, ma in enumerate(M):
            En = ma
            Ens = En - (En-e[idx]*np.sin(En) - ma)/(1 - e[idx]*np.cos(En))
            for iters in range(0,1000):
                if np.abs(Ens - En) < 1e-12:
                    break
                En = Ens
                Ens = En - (En-e[idx]*np.sin(En) - ma)/(1 - e[idx]*np.cos(En))

            # Calculate the true anomaly:
            theta[idx] = np.arctan2(np.sqrt(1- e[idx]*e[idx])*np.sin(Ens), np.cos(Ens) - e[idx])

        # Calculate the orbital momentum and radial position:
        h = np.sqrt(mu*a*(1-e*e))
        r_mag = ((h*h)/mu)*(1/(1+(e*np.cos(theta))))
        
        # Calculate the states in peri-focal coordinates:
        R_pqw = np.zeros((3,theta.size))
        R_pqw[0,:] = np.multiply(r_mag,np.cos(theta)).flatten()
        R_pqw[1,:] = np.multiply(r_mag,np.sin(theta)).flatten()

        # Calculate vectorized rotations:
        a11 =  -np.sin(node)*np.cos(i)*np.sin(peri) + np.cos(node)*np.cos(peri)
        a12 =  -np.sin(node)*np.cos(i)*np.cos(peri) - np.cos(node)*np.sin(peri)

        a21 =  np.cos(node)*np.cos(i)*np.sin(peri) + np.sin(node)*np.cos(peri)
        a22 =  np.cos(node)*np.cos(i)*np.cos(peri) - np.sin(node)*np.sin(peri)

        a31 =  np.sin(peri)*np.sin(i)
        a32 =  np.cos(peri)*np.sin(i)

        # Apply rotations to obtain position in inertial frame:
        r = np.vstack( ( (a11*R_pqw[0,:] + a12*R_pqw[1,:]),
                         (a21*R_pqw[0,:] + a22*R_pqw[1,:]),
                         (a31*R_pqw[0,:] + a32*R_pqw[1,:]) ) )

        return r



            