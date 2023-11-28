import numpy as np
import matplotlib.pyplot as plt
import math 
from math import*

r_i = float(input("Add i position vector: ")) #position i vector
r_j = float(input("Add j position vector: ")) #position j vector
r_k = float(input("Add k position vector: ")) #position k vector

r = [r_i, r_j, r_k]

r_array = np.array(r) #position vector
r_mag = np.linalg.norm(r_array) #position magnitude

v_i = float(input("Add i velocity vector: ")) #velocity i vector
v_j = float(input("Add j velocity vector: ")) #velocity j vector
v_k = float(input("Add K velocity vector: ")) #velocity k vector

v = [v_i, v_j, v_k]
v_array = np.array(v) #velocity vector
v_mag = np.linalg.norm(v_array) #velocity magnitude

K = [0,0,1]

mu = 3.986004418 * 10**5

h = np.cross(r,v) #specific angular velocity vector
h_mag = np.linalg.norm(np.array(h)) #specific angular velocity magnitude

n = np.cross(K,h)
n_vec = np.array(n)
n_mag = np.linalg.norm(np.array(n))

r_coeff = (((v_mag)**2)-mu/r_mag)/mu #position coeff for e
r_new = r_coeff * r_array

v_coeff = np.dot(r,v)/mu #velocity coeff for e
v_new = v_coeff * v_array

E = float(((v_mag*v_mag))-(mu/r_mag)) #mechanical energy

e_vec = r_new - v_new #eccentricity vector
e_mag = np.linalg.norm(np.array(e_vec)) #eccentricity magnitude

p = (h_mag**2)/mu #momentum
a = p/(1-(e_mag**2)) #semi major axis calculation

i = math.acos(h[2]/h_mag) #inclination
omega = math.acos(abs(n_vec[0])/n_mag) #ascending node
w = math.acos((np.dot(n,e_vec))/(n_mag * e_mag)) #argument of periapsis
v = math.acos((np.dot(e_vec,r))/(e_mag * r_mag)) #true anomaly

i *= (180/pi)
omega *= (180/pi)
w *= (180/pi)
v *= (180/pi)

if n_vec[1]<0:
    omega += 180

if e_vec[2]<0:
    w += 180

if np.dot(r_array,v_array)<0:
    v +=180

print(f"a = Semi-major axis is {a:.3f} km")
print(f"e = Eccentricity is {e_mag:.3f}")

if math.isnan(i):
    print("i = Angle of inclination is UNDEFINED")

else:
    print(f"i = Angle of inclination is {i:.3f}°")


if math.isnan(omega):
    print(f"Ω = Longitude of ascending node is UNDEFINED")

else:
    print(f"Ω = Longitude of ascending node is {omega:.3f}°")

if math.isnan(w):
    print(f"ω = Argument of periapsis is UNDEFINED")

else:
    print(f"ω = Argument of periapsis is {w:.3f}°")

if math.isnan(v):
    print(f"v = True anomaly is UNDEFINED")

else:
    print(f"v = True anomaly is {v:.3f}°")


