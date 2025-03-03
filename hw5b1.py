#Used ChatGPT for formatting, equation check, and debugging
#Used Dr. Smay Stem hw5b, and class lecture resources to build off of.
# region imports
import hw5a1 as pta
import random as rnd
from matplotlib import pyplot as plt
import numpy as np
# endregion

# region functions
def ffPoint(Re, rr):
    """
    Calculate the friction factor for a given Reynolds number and relative roughness.

    The function uses different approaches based on the flow regime:
      - For Re >= 4000 (turbulent flow): it uses the Colebrook equation.
      - For Re <= 2000 (laminar flow): it uses the analytical formula (f = 64/Re).
      - For transitional flow (2000 < Re < 4000): it performs a linear interpolation
        between the laminar and turbulent predictions and adds a probabilistic component.

    Parameters:
        Re (float): Reynolds number.
        rr (float): Relative roughness (ε/D).

    Returns:
        float: The computed friction factor.
    """
    if Re >= 4000:
        return pta.ff(Re, rr, CBEQN=True)
    if Re <= 2000:
        return pta.ff(Re, rr)
    CBff = pta.ff(Re, rr, CBEQN=True)  # Turbulent prediction
    Lamff = pta.ff(Re, rr)            # Laminar prediction
    mean = Lamff + (CBff - Lamff) * ((Re - 2000) / 2000)  # Linear interpolation
    sig = 0.2 * mean  # Standard deviation as 20% of mean
    return rnd.normalvariate(mean, sig)


def computeReynolds(D, epsilon, Q):
    """
    Compute the Reynolds number and related parameters based on user input.

    The function converts inputs into appropriate units, calculates velocity, and then
    computes the Reynolds number for the given pipe and flow parameters.

    Parameters:
        D (float): Pipe diameter in inches.
        epsilon (float): Pipe roughness in micro-inches.
        Q (float): Flow rate in gallons per minute (GPM).

    Returns:
        tuple: A tuple containing:
            - Re (float): Reynolds number.
            - rr (float): Relative roughness (ε/D).
            - V (float): Flow velocity in ft/s.
            - D_ft (float): Pipe diameter in feet.
    """
    D_ft = D / 12  # Convert inches to feet
    epsilon_in = epsilon * 1e-6  # Convert micro-inches to inches
    rr = epsilon_in / D_ft  # Relative roughness (ε/D)
    V = (Q / 448.831) / (np.pi * (D_ft / 2) ** 2)  # Convert GPM to ft/s
    nu = 1.08e-5  # Kinematic viscosity of water in ft²/s
    Re = (V * D_ft) / nu
    return Re, rr, V, D_ft


def computeHeadLoss(f, V, D):
    """
    Compute head loss per foot using the Darcy-Weisbach equation.

    Parameters:
        f (float): Friction factor.
        V (float): Flow velocity in ft/s.
        D (float): Pipe diameter in feet.

    Returns:
        float: Head loss per foot (ft/ft).
    """
    g = 32.2  # ft/s² (acceleration due to gravity)
    return f * (V ** 2 / (2 * g * D))


def main():
    """
    Main function for collecting user input, computing hydraulic parameters, and plotting results.

    Workflow:
      1. Repeatedly prompt the user for pipe diameter, roughness, and flow rate.
      2. For each input set, compute the Reynolds number, friction factor, and head loss.
      3. Store the computed (Re, friction factor, transition flag) in a list.
      4. After input is complete, draw the Moody chart using pta.plotMoody().
      5. Overlay the stored data points onto the Moody chart.
      6. Display the final combined figure.

    Note:
      - This function temporarily overrides plt.show() to prevent pta.plotMoody() from
        displaying its figure prematurely, allowing us to overlay our data points.
    """
    data_points = []  # List to store (Re, friction factor, transition flag) tuples

    while True:
        try:
            D = float(input("Enter pipe diameter (in inches): "))
            epsilon = float(input("Enter pipe roughness (in micro-inches): "))
            Q = float(input("Enter flow rate (in gallons per minute, GPM): "))
        except ValueError:
            print("Invalid input. Please enter numerical values.")
            continue

        Re, rr, V, D_ft = computeReynolds(D, epsilon, Q)
        f = ffPoint(Re, rr)
        hf_per_L = computeHeadLoss(f, V, D_ft)

        print(f"Reynolds Number: {Re:.2f}")
        print(f"Friction Factor: {f:.5f}")
        print(f"Head Loss per Foot: {hf_per_L:.5f} ft/ft")

        # Determine if the flow is in the transition region
        transition = 2000 < Re < 4000
        data_points.append((Re, f, transition))

        cont = input("Would you like to enter another set of parameters? (y/n): ")
        if cont.lower() != 'y':
            break

    # Temporarily override plt.show() to prevent premature display by pta.plotMoody()
    original_show = plt.show
    plt.show = lambda: None

    # Draw the Moody chart using the pta module
    pta.plotMoody(plotPoint=False)

    # Restore the original plt.show() function
    plt.show = original_show

    # Grab the current axis (used by pta.plotMoody)
    ax = plt.gca()

    # Overlay all stored data points onto the Moody chart
    for Re, f, transition in data_points:
        marker = '^' if transition else 'o'
        ax.scatter(Re, f, marker=marker, color='r', s=100, edgecolors='black', zorder=3)

    # Display the final combined figure
    plt.show()

# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion

