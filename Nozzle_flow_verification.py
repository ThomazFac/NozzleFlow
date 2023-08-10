#!/usr/bin/env python
# coding: utf-8

# # Compressible Flow through Nozzle

import numpy as np
import sympy as smp
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from scipy.optimize import minimize


def CubicCurve(coef, x):
    Dcr = coef[0]
    Din = coef[1]
    Xm  = coef[2]
    L   = coef[3]
    
    D = 0
    
    if x/L < Xm:
        D = (Din - Dcr)*(1-(1/Xm**2)*(x/L)**3) + Dcr
    else:
        D = (Din - Dcr)*((1/(1 - Xm)**2)*(1-x/L)**3) + Dcr
    
    return D


def RingArea(R,r):
    
    A = np.pi*(R+r)*(R-r)
    
    return A


def AreaConv(x):
    Dwall_in = 0.160
    Dcent_in = 0.105
    Dwall_t  = 0.046
    Dcent_t  = 0.024
    Xm       = 0.45
    L        = 0.171
    
    Dwall_x = CubicCurve([Dwall_t, Dwall_in, Xm, L], x)
    Dcent_x = CubicCurve([Dcent_t, Dcent_in, Xm, L], x)
    
    Rwall_x = Dwall_x/2
    Rcent_x = Dcent_x/2
    
    A = RingArea(Rwall_x,Rcent_x)
    return A


def AreaDiv(x):
    Rwall_t  = 0.046/2
    Rcent_t  = 0.024/2
    x_t = 0.171
    
    Rwall_out = 0.046/2
    Rcent_out = 0.0086/2
    x_out = 0.496    
    
    slope = (Rcent_out - Rcent_t) / (x_out - x_t)
    
    y_intercept = Rcent_t - slope * x_t
    
    Rcent_x = x*slope + y_intercept
    
    Rwall_x = Rwall_out
    
    A = RingArea(Rwall_x,Rcent_x)
    
    return A


def Area_Mach(const, x):
    # substituir A por x
    gam   = const[0]
    R     = const[1]
    A_t   = const[2]
    rho_t = const[3]
    p_t   = const[4]
    T_t   = const[5]
    x_t   = const[6]

    if x < x_t:
        A = AreaConv(x)
        M0 = .001
    else:
        A = AreaDiv(x)
        M0 = 1.001
    
    rho, p, T, M = smp.symbols('\rho, p, T, M', real=True, positive=True)

    LeftSide = (A/A_t)**2
    RightSide = 1/M**2*(2/(gam+1)*(1+(gam-1)/2*M**2))**((gam+1)/(gam-1))

    area_mach = LeftSide-RightSide
    
    lam_f = smp.lambdify(M, area_mach**2)    
    
    res = minimize(lam_f, x0=M0, method='L-BFGS-B')
    M_sol = res.x[0]
    
    p   = p_0*(1+(gam-1)/2*M_sol**2)**(-gam/(gam-1))
    T   = T_0*(1+(gam-1)/2*M_sol**2)**-1
    rho = rho_0*(1+(gam-1)/2*M_sol**2)**(-1/(gam-1))
    
    return M_sol, p, T, rho


def ReadCFDResult(ResultFileName):
    file = open(ResultFileName,'r')

    lines = file.readlines()
    nLines = len(lines)
    file.close()

    xLists = [[]]
    yLists = [[]]

    i = 0
    flag = 1
    for j in range(nLines):
        line = lines[j]
        data = line.split(',')
        try:
            xLists[i].append(float(data[0]))
            yLists[i].append(float(data[1]))
            flag = 0
        except:
            if flag == 0:
                i+=1
                xLists.append([])
                yLists.append([])
                flag = 1
            pass
    return xLists, yLists


X = np.linspace(0.,0.496,50)

gamma_CH4 = 1.289
R_CH4 = 518.2 # J/kg.K

global p_0, T_0, rho_0

p_0  = 4e6    # Pa
T_0 = 300 # K
rho_0 = p_0/(R_CH4*T_0) # kg/m^3

A_t = 0.001209 # m^2
x_t = 0.171 # m

p_t = p_0/(1+(gamma_CH4-1)/2*(1)**2)**(gamma_CH4/(gamma_CH4-1))
T_t = T_0*(p_t/p_0)**((gamma_CH4-1)/gamma_CH4)
rho_t = p_t/(R_CH4*T_t)

InputList = [gamma_CH4, R_CH4, A_t, rho_t, p_t, T_t, x_t]

Sol_list = []

for x in X:
    Sol_list.append(Area_Mach(InputList, x))
    
Sol_list = np.array(Sol_list)
Sol_list = Sol_list.T

ResultFileNameM = "Mach_redlich.csv"
ResultFileNameT = "temp_redlich.csv"
ResultFileNameP = "pressure_redlich.csv"

xLists, MachLists = ReadCFDResult(ResultFileNameM)
xLists, TemperatureLists = ReadCFDResult(ResultFileNameT)
xLists, PressureLists = ReadCFDResult(ResultFileNameP)

xlist = np.array(xLists[0])
MachList = np.array(MachLists[0])
TemperatureList = np.array(TemperatureLists[0])
PressureList = np.array(PressureLists[0])

index = np.argmax(xlist > 0.496)


font = {'family' : 'Times New Roman',
        'size'   : 20}

mpl.rc('font', **font)


fig = plt.figure(figsize=(10,5))
plt.plot(xlist[:index],MachList[:index],'-',linewidth=3, label='CFD')
plt.plot(X, Sol_list[0],'--',linewidth=3, label='quasi-1D')
plt.plot([0.171,0.171], [min(Sol_list[0]), max(Sol_list[0])],'--',color='grey')
plt.text(0.18,0.3,'Throat', color='grey')

plt.text(0.21,0.7,r'$\left ( \frac{A}{A^*}\right )^{2} = \frac{1}{M^2} \left[\frac{2}{\gamma+1} \left(  1+\frac{\gamma-1}{2} M^2\right) \right]^{\frac{(\gamma+1)}{(\gamma -1)}}$', math_fontfamily='dejavuserif')

plt.grid()
plt.xlabel(r'x (m)')
plt.ylabel(r'Mach Number')
plt.legend()
plt.tight_layout()
plt.show()


fig = plt.figure(figsize=(10,5))
plt.plot(xlist[:index],TemperatureList[:index],'-',markersize=3, label='CFD')
plt.plot(X, Sol_list[2],'-',markersize=3, label='quasi-1D')
plt.plot([0.171,0.171], [210, max(Sol_list[2])],'--',color='grey')
plt.text(0.11,230,'Throat', color='grey')
plt.grid()
plt.xlabel(r'x (m)')
plt.ylabel(r'Static Temperature (K)')
plt.legend()
plt.tight_layout()
plt.show()


fig, ax1 = plt.subplots(figsize=(10,5))
ax2 = ax1.twinx()

ax1.plot(xlist[:index],PressureList[:index],'-',markersize=3, label='CFD', color = 'blue')
ax1.plot(X, Sol_list[1],'--',markersize=3, label='quasi-1D', color = 'blue')

ax2.plot(xlist[:index],MachList[:index],'-',markersize=3, label='CFD', color = 'green')
ax2.plot(X, Sol_list[0],'--',markersize=3, label='quasi-1D')

ax2.plot([0.171,0.171], [min(Sol_list[0]), max(Sol_list[0])],'--',color='grey')
ax2.text(0.11,0.13,'Throat', color='grey')

ax1.grid()
ax1.set_xlabel('x (m)')
ax1.set_ylabel('Static Pressure (Pa)', color = 'blue')
ax2.set_ylabel('Mach Number', color = 'green')
plt.legend(loc = 'center right')
plt.tight_layout()
plt.show()

