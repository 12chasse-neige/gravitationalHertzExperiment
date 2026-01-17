import numpy as np
from scipy.integrate import tplquad
from scipy.optimize import approx_fprime
import math

#define the contants that will be used
H = 2   #height of the column (meter)
D = 5   #diameter of the column (meter)
d = 1   #diameter of the holes (meter)
s = 1.5 #distance from the center of the colum to the center of the holes (meter)
R = 2000    #distance from the source to the detector (meter)
rho = 1750  #the density of typical carbonfiber (kg/m^3)
G = 6.674E-11   #gravitational constant (m^3/kg*s^2)
c = 2.998E8 #speed of light (m/s)
cos = math.cos
sin = math.sin
pi = math.pi
omega = 20000 * 2 * pi    #the frequency of the rotation (rad/s)
def dirac_delta(i, j):
        if i==j:
            return 1
        else:
            return 0    #dirac delta function

#describe the rotation of the columns: we assume that at the time t=0, the center of one hole is at the direction of the x axis, and the center of the four holes are defined by x_k and y_k (k = 1,2,3,4)
def get_the_coordinate_of_the_hole(k, t):
    x_k = s * cos(omega * t + k * pi/2) 
    y_k = s * sin(omega * t + k * pi/2)
    return x_k, y_k

#calculate the components of the quadrupole tensor
#calculate the ij componet of the quadrupole tensor of one column
def calculate_one_column_tensor(i, j, D, rho, H):
    def y_lower(x):
        return -np.sqrt(D**2/4 - x**2)
    def y_upper(x):
        return np.sqrt(D**2/4 - x**2)
    def integrand(x, y, z):
        coordinates = [x, y, z]
        integrand = rho *(coordinates[i] * coordinates[j] - 1/3 * dirac_delta(i, j) * (x**2 + y**2 + z**2))
        return integrand
    component = tplquad(
            integrand,
            -D/2, D/2,                    
            y_lower, y_upper,         
            lambda x, y: -H/2,           
            lambda x, y: H/2  
        )[0]
    return component

#calculate the componets of quadrupole tensor of the holes relative to the center of the column
def hole_tensor_component(i, j, k, d, rho ,H ,t):
    component = calculate_one_column_tensor(i, j, d, rho, H)
    x = get_the_coordinate_of_the_hole(k, t)[0]
    y = get_the_coordinate_of_the_hole(k, t)[1]
    z = 0
    def relative_component(i, j, t):
        coordinates = [x, y, z]
        relativeComponent = rho * pi * d**2 / 4 * H * (coordinates[i] * coordinates[j] - 1/3 * dirac_delta(i, j) * (x**2 + y**2 + z**2))
        return relativeComponent
    component += relative_component(i, j, t)
    return component

#adding to get the whole quadrupole tensor I_ij
def calculate_the_whole_tensor(i, j, D, rho, H, t):
    big_column = calculate_one_column_tensor(i, j, D, rho, H)
    component = big_column
    for k in range(4):
        small_column = hole_tensor_component(i, j, k, d, -rho, H, t)
        component += small_column
    return component

#calculate the second-order derivative of the quadrupole tensor
def second_derivative_of_component(i, j, t):
    if i==2 or j==2:
        return 0
    else:
        derivative = - 4 * omega ** 2 * calculate_the_whole_tensor(i, j, D, rho, H, t)
        return derivative

#get the final metric tensor
t = float(input("current time:"))
t_rev = t - R/c
metricTensorComponent = [[0,0,0],[0,0,0],[0,0,0]]
for i in range(3):
    for j in range(3):
        metricTensorComponent[i][j] = 2 * G * second_derivative_of_component(i, j, t_rev) / (R * c**4)
        print(metricTensorComponent[i][j], end = " ")
    print("\t")

#project the mectric tensor to TT mode
def project_to_TT_mode(metricTensorComponent):
    projector = [[1,0,0],[0,1,0],[0,0,0]]
    metricTensorComponentTT = [[0,0,0],[0,0,0],[0,0,0]]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    metricTensorComponentTT[i][j] += (projector[i][k] * projector[j][l] - 0.5 * projector[i][j] * projector[k][l]) * metricTensorComponent[k][l]
    return metricTensorComponentTT

#output the final
metricTensorComponentTT = project_to_TT_mode(metricTensorComponent)
for i in range(3):
    for j in range(3):
        print(metricTensorComponentTT[i][j], end = " ")
    print("\t")

