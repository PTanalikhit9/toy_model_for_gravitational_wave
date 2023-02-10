import numpy as np
import matplotlib.pyplot as plt

def obtain_root(polynomial):
    """
    This function calculates the roots of the polynomial equation with coefficients provided in 'polynomial' list
    The polynomial equation of the form: ax^4+bx^3+cx^2+dx+e=0 will be solved using Ferrari's method.
    
    Parameters:
    polynomial (tuple): A tuple of 5 coefficients (a, b, c, d, e) of the polynomial equation ax^4 + bx^3 + cx^2 + dx + e = 0.
    
    Returns:
    root (ndarray): An array of 4 roots of the polynomial equation.
    """
    
    a, b, c, d, e = polynomial # Coefficients of polynomial equation ax^4+bx^3+cx^2+dx+e=0
    
    # Intermediate variables for calculation
    p = (8*a*c-3*b*b)/(8*a*a)
    q = (b*b*b-4*a*b*c+8*a*a*d)/(8*a*a*a)
    delta0 = c*c-3*b*d+12*a*e
    delta1 = 2*c*c*c-9*b*c*d+27*b*b*e+27*a*d*d-72*a*c*e
    delta = (4*delta0*delta0*delta0-delta1*delta1)/27
    
    # Finding roots based on the value of delta
    if delta > 0:
        phi = np.arccos(delta1/(2*np.sqrt(delta0*delta0*delta0)))
        S = 0.5*np.sqrt(-(2/3)*p+(2/(3*a))*np.sqrt(delta0)*np.cos(phi/3.0))
    else:
        Q = np.cbrt((delta1+np.sqrt(delta1*delta1-4*delta0*delta0*delta0))/2.0)
        S = 0.5*np.sqrt(-(2/3)*p+(1/(3*a))*(Q+delta0/Q))
    rh1 = 1/2*np.sqrt(-4*S*S-2*p+q/S)
    rh2 = 1/2*np.sqrt(-4*S*S-2*p-q/S)
    
    # Calculation of roots
    root = np.array([-b/(4*a)-S+rh1, -b/(4*a)-S-rh1, -b/(4*a)+S+rh2, -b/(4*a)+S-rh2])
    return root

#This is to check whether the implemented function is working or not
###################################################
#print(obtain_root([3, 6, -123, -126, 1080]))
#print(obtain_root([-20, 5, 17, -29, 87]))
#print(obtain_root([2, 4, 6, 8, 10]))
###################################################


# This function calculates the value of beta for each radius
def get_beta(param,radius,Mass):
    """
    This function calculates the beta value for a given set of parameters, radius and Mass.
    
    Parameters:
    param (dict): A dictionary of parameters that includes the alpha and density value.
    radius (ndarray): An array of radii.
    Mass (float): The mass of the ball.
    
    Returns:
    beta (ndarray): An array of beta values corresponding to each radius.
    """
    
    # Initialize an array to store beta values
    beta = np.zeros(len(radius))
    # Gravitational constant
    g = 9.807
    # Get the value of alpha from the input parameters
    alpha = param['alpha']
    # Get the value of density from the input parameters
    density = param['density']
    
    # Loop through all the values of radius
    for i in range(len(radius)):
        # Calculate beta for each radius
        beta[i] = (alpha/radius[i])*(Mass+np.pi*density*radius[i]*radius[i])
        
    # Return the beta values
    return beta

# This function generates the fabric slope for each radius
def gen_fabric_slope(param,radius, Mass, beta,root_number):
    """
    This function generates the fabric slope for each radius.

    Parameters:
    param (dict): A dictionary of parameters. It should contain the following key-value pairs:
    alpha: The value of alpha.
    density: The value of density.
    radius (ndarray): An array of radii.
    Mass (float): The mass of the object.
    beta (ndarray): An array of beta values.
    root_number (int): An integer indicating which root of the polynomial to use in the calculation.
    
    Returns:
    slope (ndarray): An array of fabric slope values.
    """
    
    # Gravitational constant
    g = 9.807
    # Initialize an array to store the fabric slope values
    slope = np.zeros(len(radius))
    
    # Loop through all the values of radius
    for i in range(len(radius)):
        # Calculate the fabric slope for each radius
        slope[i] = obtain_root([1,-2*beta[i],beta[i]*beta[i],-2*beta[i],beta[i]*beta[i]])[root_number-1]
        
    # Return the fabric slope values
    return slope

# This function reads the fabric slope from a file
def read_fabric_slope(filename):
    """
    The function read_fabric_slope reads the fabric slope values from a file and returns the values.

    Parameters:
    filename (str): The name of the file containing the fabric slope values.
    
    Returns:
    slope (ndarray): An array of fabric slope values read from the file.
    """
    
    # Read the slope values from the file
    slope = filename
    # Return the slope values
    return slope

# This function calculates the fabric height for each radius
def fabric_height(param,radius,slope,z0=0.001):
    """
    This function calculates the fabric height for each radius by using a numerical integration method, 
    specifically the 4th-order Runge-Kutta method. 
    
    Parameters:
    param: a dictionary containing values such as rmax, rmin, and grid_number.
    radius: an array of radii values.
    slope: an array of slope values.
    z0: a default value for the first height, which is 0.001.
    
    Returns:
    height: an array of calculated height values.
    
    Note that the fabric_height function uses the Runge-Kutta method to update the height values for each radius.
    """
    
    # Initialize an array to store the fabric height values
    height = np.zeros(len(radius))
    # Set the first value of height to a default value of 0.001
    height[0] = z0
    # Calculate the step size for the radius values
    h = (param["rmax"]-param["rmin"])/(param["grid_number"])
    
    # Loop through all the values of radius
    for i in range(len(radius)-1):
        # Calculate the intermediate values for height
        """
        The function loops through all the values of radius (excluding the last value), 
        updating the height values using the 4th-order Runge-Kutta method. The method involves 
        calculating intermediate values of height (k1, k2, k3, k4) using the slope values, 
        and updating the height values using these intermediate values.
        """
        k1 = (slope[i])*h
        k2 = (slope[i]+slope[i+1])*h*0.5
        k3 = (slope[i]+slope[i+1])*h*0.5
        k4 = slope[i+1]*h
        # Update the height values using a Runge-Kutta method
        height[i+1] = height[i]+(k1+2*k2+2*k3+k4)/6.0
        
    # Return the fabric height values
    return height


def Effective_Potential(param,radius,beta,slope,height,ball_mass,Angular_Momentum,Hamiltonian,Mass):
    """
    Effective_Potential is a function that calculates the effective potential for each radius, given the parameters, 
    beta, slope, height, ball mass, angular momentum, Hamiltonian, and mass.

    Parameters:
    param: a dictionary containing values such as rmax and rmin.
    radius: an array of radii values.
    beta: an array of beta values.
    slope: an array of slope values.
    height: an array of height values.
    ball_mass: the mass of the ball.
    Angular_Momentum: the angular momentum of the ball.
    Hamiltonian: the Hamiltonian of the system.
    Mass: the mass of the system.

    Returns:
    V_eff: an array of calculated effective potential values.

    Note that the Effective_Potential function plots the effective potential curve 
    and the Hamiltonian as a line and displays it, and returns the effective potential values.
    """
    
    # Gravitational constant
    g = 9.807
    # Initialize an array to store the effective potential values
    V_eff = np.zeros(len(radius))
    # Initialize an array to store the potential values
    potential_V = np.zeros(len(radius))
    # Initialize an array to store the potential energy due to height
    potentialZ = np.zeros(len(radius))
    # Initialize an array to store the angular momentum values
    angular_mom = np.zeros(len(radius))
    
    # Loop through all the values of radius
    for i in range(len(radius)):
        # Calculate the potential energy due to height
        potentialZ[i]= (5.0/7.0)*ball_mass*9.807*height[i]
        # Calculate the angular momentum values
        angular_mom[i] = (0.5/ball_mass)*(Angular_Momentum/radius[i])*(Angular_Momentum/radius[i])
        # Calculate the potential values
        potential_V[i] = angular_mom[i] + potentialZ[i]
        # Calculate the effective potential values
        V_eff[i] = potential_V[i]*(1-2*(beta[i]/slope[i])+(beta[i]/slope[i])**2)+2*(Hamiltonian*beta[i]/slope[i])-Hamiltonian*(beta[i]/slope[i])**2
    
    
    #Display the graph of effective potential
    # Set the starting point for the plot
    startingpoint = 600
    # Plot the effective potential curve
    plt.axes(title='Effective Potential Curve',xlabel='radius',ylabel='energy')
    plt.plot(radius[startingpoint:],V_eff[startingpoint:],'k', linewidth=3)
    # Plot the Hamiltonian as a line
    plt.plot(radius[startingpoint:],Hamiltonian*np.ones(len(radius)-startingpoint),'orange', linewidth=3)
    # Return the effective potential values
    return V_eff



def solve_orbit(param,radius,V_eff,Hamiltonian,Angular_Momentum,ball_mass, N):
    """
    This function solves the orbit of a ball in a gravitational potential.
    
    Parameters:
    param (dict): A dictionary of parameters.
    radius (ndarray): An array of radii.
    V_eff (ndarray): An array of effective potentials.
    Hamiltonian (float): The Hamiltonian of the system.
    Angular_Momentum (float): The angular momentum of the system.
    ball_mass (float): The mass of the ball.
    N (int): The number of cycles.
    
    Returns:
    phi (ndarray): An array of angles.
    new_radius (ndarray): An array of updated radii.
    """
    
    # Get the location where the effective potential is less than the Hamiltonian
    location = np.where(V_eff < Hamiltonian)
    # Get the corresponding radius and effective potential
    new_radius,new_V_eff  = radius[location],V_eff[location]
    # Calculate the difference between the Hamiltonian and the effective potential
    diff_V_eff = Hamiltonian-new_V_eff
    # Select the radii corresponding to the location where the effective potential is less than the Hamiltonian
    BC = np.select([V_eff < Hamiltonian],[radius])
    # Calculate the change in radius with respect to the change in angle
    drdphi = (new_radius*new_radius/Angular_Momentum)*np.sqrt(2*ball_mass*diff_V_eff)
    # Initialize an array of angles
    phi = np.zeros(len(diff_V_eff)) #in unit radian
    
    # Set the first angle to 0
    phi[0] = 0
    
    
    # Loop through the array of difference in effective potential
    for i in range(len(diff_V_eff)-1):
        # Calculate the change in angle by integrating the change in radius with respect to the change in angle
        phi[i+1] = phi[i]+(new_radius[i+1]-new_radius[i])/drdphi[i]
    
    
    # Stack the original radius and its reverse to form a complete cycle
    new_radius = np.hstack((new_radius,np.flip(new_radius)))
    new_radiuss = new_radius
    # Calculate the length of a cycle
    cycle_length = 2*len(phi)
    # Stack the original angles and its shifted version to form a complete cycle
    phi = np.hstack((phi, np.ones(len(phi))*phi[len(diff_V_eff)-1]+phi))
    phii = phi
    
    
    # Repeat the cycle to form the desired number of cycles
    for i in range(N-1):
        new_radius = np.hstack((new_radius,new_radiuss))
        phi = np.hstack((phi, np.ones(cycle_length)*phi[-1] + phii))
        
        
    # Plot the orbit on a polar plot
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

#main
#Parameters
Mass = 0.588 #in unit kilogram 
beta = get_beta(param,radius,Mass)
slope = gen_fabric_slope(param,radius,Mass,beta,3)
slope1 = gen_fabric_slope(param,radius,Mass,beta,4)
height = fabric_height(param,radius,slope,0.001)
ball_mass = 0.002113 #in kg
Hamiltonian = 0.0039 #0.00328 #-0.02576 #in Joule
Angular_Momentum= 0.000272 #0.000284 #0.00272 #in Kg m^2 s
V_eff = Effective_Potential(param,radius,beta,slope,height,ball_mass,Angular_Momentum,Hamiltonian,Mass)
N = 20
a = solve_orbit(param,radius,V_eff,Hamiltonian,Angular_Momentum,ball_mass, N)

plt.show()
