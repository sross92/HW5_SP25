#Used ChatGPT for formatting, equation check, and debugging
#Used Dr. Smay Stem hw5c, and class lecture resources to build off of.
# region imports
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# endregion

# region functions
def ode_system(t, X, *params):
    """
    Defines the system of ODEs for a hydraulic valve/piston system.

    The state vector X is defined as follows:
      X[0] = x      : Position of the piston.
      X[1] = xdot   : Velocity of the piston.
      X[2] = p1     : Pressure on the right side of the piston.
      X[3] = p2     : Pressure on the left side of the piston.

    The equations governing the system are:
      1. dx/dt = xdot
      2. d(xdot)/dt = (p1 - p2) * A / m
      3. dp1/dt = { y*Kvalve*(ps - p1) - rho*A*xdot } * (beta/(rho*V))
      4. dp2/dt = -{ y*Kvalve*(p2 - pa) - rho*A*xdot } * (beta/(rho*V))

    Parameters
    ----------
    t : float
        Current time (required by solve_ivp, but not explicitly used in the equations).
    X : array-like, shape (4,)
        State vector [x, xdot, p1, p2].
    params : tuple
        Tuple of physical constants in the following order:
        (A, Cd, ps, pa, V, beta, rho, Kvalve, m, y)
          - A      : Cross-sectional area (m^2)
          - Cd     : Discharge coefficient (not explicitly used in the equations)
          - ps     : Supply pressure (Pa)
          - pa     : Ambient pressure (Pa)
          - V      : Effective fluid volume (m^3)
          - beta   : Bulk modulus (Pa)
          - rho    : Fluid density (kg/m^3)
          - Kvalve : Valve flow coefficient
          - m      : Mass of the piston (kg)
          - y      : Constant valve input

    Returns
    -------
    list
        Derivatives of the state variables [dx/dt, d^2x/dt^2, dp1/dt, dp2/dt].
    """
    # Unpack the parameters
    A, Cd, ps, pa, V, beta, rho, Kvalve, m, y = params

    # Unpack the state variables for clarity
    x    = X[0]
    xdot = X[1]
    p1   = X[2]
    p2   = X[3]

    # Compute the derivatives using the provided equations
    # Piston acceleration
    xddot = (p1 - p2) * A / m

    # Pressure rate on the right side of the piston
    p1dot = ( y * Kvalve * (ps - p1) - rho * A * xdot ) * (beta / (rho * V))

    # Pressure rate on the left side of the piston
    p2dot = - ( y * Kvalve * (p2 - pa) - rho * A * xdot ) * (beta / (rho * V))

    # Return the derivatives as a list in the same order as the state variables
    return [xdot, xddot, p1dot, p2dot]


def main():
    """
    Main function that solves the hydraulic valve ODE system and plots the results.

    Workflow:
      1. Define the time span for the simulation.
      2. Set the physical parameters and initial conditions:
           - x(0) = 0        (initial piston position)
           - xdot(0) = 0     (initial piston velocity)
           - p1(0) = pa      (initial pressure on the right side)
           - p2(0) = pa      (initial pressure on the left side)
      3. Solve the system of ODEs using solve_ivp.
      4. Extract and plot:
           a. The piston velocity xdot versus time.
           b. Both pressures p1 and p2 versus time on a separate subplot.
    """
    # Define the time vector for the simulation
    t = np.linspace(0, 0.02, 200)

    # Define the physical constants and parameters as a tuple:
    # (A, Cd, ps, pa, V, beta, rho, Kvalve, m, y)
    myargs = (4.909E-4, 0.6, 1.4E7, 1.0E5, 1.473E-4, 2.0E9, 850.0, 2.0E-5, 30.0, 0.002)

    # Extract ambient pressure pa from the parameters (4th element)
    pa = myargs[3]

    # Set initial conditions: x = 0, xdot = 0, p1 = pa, p2 = pa
    ic = [0.0, 0.0, pa, pa]

    # Solve the ODE system using solve_ivp
    sln = solve_ivp(
        fun=lambda t_, X_: ode_system(t_, X_, *myargs),
        t_span=(t[0], t[-1]),
        y0=ic,
        t_eval=t
    )

    # Unpack the solution
    xvals = sln.y[0]
    xdot  = sln.y[1]
    p1    = sln.y[2]
    p2    = sln.y[3]

    # Plotting results
    # Top subplot: Piston position and velocity
    plt.subplot(2, 1, 1)
    plt.title("Piston Position and Velocity vs Time")
    plt.plot(t, xvals, 'r-', label='$x$')
    plt.ylabel('Piston Position, $x$ (m)')
    plt.legend(loc='upper left')

    # Create a twin axis for velocity on the same subplot
    ax2 = plt.twinx()
    ax2.plot(t, xdot, 'b-', label='$\\dot{x}$')
    ax2.set_ylabel('Piston Velocity, $\\dot{x}$ (m/s)')
    ax2.legend(loc='lower right')

    # Bottom subplot: Pressures p1 and p2
    plt.subplot(2, 1, 2)
    plt.title("Pressures $p_1$ and $p_2$ vs Time")
    plt.plot(t, p1, 'b-', label='$p_1$')
    plt.plot(t, p2, 'r-', label='$p_2$')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (Pa)')
    plt.legend(loc='lower right')

    plt.tight_layout()
    plt.show()

# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion
