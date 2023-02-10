import numpy as np
import matplotlib.pyplot as plt

#This function solves the polynomial equation of the form: ax^4+bx^3+cx^2+dx+e=0 using Ferrari's method
def obtain_root(polynomial):
    a, b, c, d, e = polynomial
    p = (8*a*c-3*b*b)/(8*a*a)
    q = (b*b*b-4*a*b*c+8*a*a*d)/(8*a*a*a)
    delta0 = c*c-3*b*d+12*a*e
    delta1 = 2*c*c*c-9*b*c*d+27*b*b*e+27*a*d*d-72*a*c*e
    delta = (4*delta0*delta0*delta0-delta1*delta1)/27
    if delta > 0:
        phi = np.arccos(delta1/(2*np.sqrt(delta0*delta0*delta0)))
        S = 0.5*np.sqrt(-(2/3)*p+(2/(3*a))*np.sqrt(delta0)*np.cos(phi/3.0))
    else:
        Q = np.cbrt((delta1+np.sqrt(delta1*delta1-4*delta0*delta0*delta0))/2.0)
        S = 0.5*np.sqrt(-(2/3)*p+(1/(3*a))*(Q+delta0/Q))
    rh1 = 1/2*np.sqrt(-4*S*S-2*p+q/S)
    rh2 = 1/2*np.sqrt(-4*S*S-2*p-q/S)
    root = np.array([-b/(4*a)-S+rh1, -b/(4*a)-S-rh1, -b/(4*a)+S+rh2, -b/(4*a)+S-rh2])
    return root

#print(obtain_root([3, 6, -123, -126, 1080]))
#print(obtain_root([-20, 5, 17, -29, 87]))
#print(obtain_root([2, 4, 6, 8, 10]))



def get_beta(param,radius,Mass):
    beta = np.zeros(len(radius))
    g = 9.807
    alpha = param['alpha']
    density = param['density']
    for i in range(len(radius)):
        beta[i] = (alpha/radius[i])*(Mass+np.pi*density*radius[i]*radius[i])
    return beta

def gen_fabric_slope(param,radius, Mass, beta,root_number):
    g = 9.807
    slope = np.zeros(len(radius))
    for i in range(len(radius)):
        slope[i] = obtain_root([1,-2*beta[i],beta[i]*beta[i],-2*beta[i],beta[i]*beta[i]])[root_number-1]
    return slope

def read_fabric_slope(filename):
    slope = filename
    return slope

def fabric_height(param,radius,slope,z0=0.001):
    height = np.zeros(len(radius))
    height[0] = z0
    h = (param["rmax"]-param["rmin"])/(param["grid_number"])
    for i in range(len(radius)-1):
        k1 = (slope[i])*h
        k2 = (slope[i]+slope[i+1])*h*0.5
        k3 = (slope[i]+slope[i+1])*h*0.5
        k4 = slope[i+1]*h
        height[i+1] = height[i]+(k1+2*k2+2*k3+k4)/6.0
    return height

def Effective_Potential(param,radius,beta,slope,height,ball_mass,Angular_Momentum,Hamiltonian,Mass):
    g = 9.807
    V_eff = np.zeros(len(radius))
    potential_V = np.zeros(len(radius))
    potentialZ = np.zeros(len(radius))
    angular_mom = np.zeros(len(radius))
    for i in range(len(radius)):
        potentialZ[i]= (5.0/7.0)*ball_mass*9.807*height[i]
        angular_mom[i] = (0.5/ball_mass)*(Angular_Momentum/radius[i])*(Angular_Momentum/radius[i])
        potential_V[i] = angular_mom[i] + potentialZ[i]
        V_eff[i] = potential_V[i]*(1-2*(beta[i]/slope[i])+(beta[i]/slope[i])**2)+2*(Hamiltonian*beta[i]/slope[i])-Hamiltonian*(beta[i]/slope[i])**2
    #print(height)
    #startingpoint = 600
    #plt.axes(title='Effective Potential Curve',xlabel='radius',ylabel='energy')
    #plt.plot(radius[startingpoint:],V_eff[startingpoint:],'k', linewidth=3)
    #plt.plot(radius[startingpoint:],Hamiltonian*np.ones(len(radius)-startingpoint),'orange', linewidth=3)
    return V_eff


def solve_orbit(param,radius,V_eff,Hamiltonian,Angular_Momentum,ball_mass, N):
    location = np.where(V_eff < Hamiltonian)
    new_radius,new_V_eff  = radius[location],V_eff[location]
    diff_V_eff = Hamiltonian-new_V_eff
    BC = np.select([V_eff < Hamiltonian],[radius])
    drdphi = (new_radius*new_radius/Angular_Momentum)*np.sqrt(2*ball_mass*diff_V_eff)
    phi = np.zeros(len(diff_V_eff)) #in unit radian
    phi[0] = 0
    for i in range(len(diff_V_eff)-1):
        phi[i+1] = phi[i]+(new_radius[i+1]-new_radius[i])/drdphi[i]
    new_radius = np.hstack((new_radius,np.flip(new_radius)))
    new_radiuss = new_radius
    cycle_length = 2*len(phi)
    phi = np.hstack((phi, np.ones(len(phi))*phi[len(diff_V_eff)-1]+phi))
    phii = phi
    for i in range(N-1):
        new_radius = np.hstack((new_radius,new_radiuss))
        phi = np.hstack((phi, np.ones(cycle_length)*phi[-1] + phii))
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(phi, new_radius)
    ax.set_rmax(param['rmax'])
    ax.set_rticks([param['rmax']/2.0,param['rmax']])  # Less radial ticks
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)

    ax.set_title("A line plot on a polar axis", va='bottom')
    return phi, new_radius

param = {"density": 0.197, "elasticity": 112, "rmin": 0.0001, "rmax":1.0,"grid_number":20000, "alpha": 0.016532}
radius = np.linspace(param['rmin'],param['rmax'],param['grid_number'])
#print(radius)
Mass = 0.588 #in unit kilogram 
beta = get_beta(param,radius,Mass)
slope = gen_fabric_slope(param,radius,Mass,beta,3)
slope1 = gen_fabric_slope(param,radius,Mass,beta,4)
height = fabric_height(param,radius,slope,0.001)
ball_mass = 0.002113 #in kg
Hamiltonian = 0.0039#0.00328#-0.02576 #in Joule
Angular_Momentum= 0.000272#0.000284#0.00272 #in Kg m^2 s
V_eff = Effective_Potential(param,radius,beta,slope,height,ball_mass,Angular_Momentum,Hamiltonian,Mass)
N = 20
a = solve_orbit(param,radius,V_eff,Hamiltonian,Angular_Momentum,ball_mass, N)
plt.show()
