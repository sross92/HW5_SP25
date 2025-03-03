#Used ChatGPT for formatting, equation check, and debugging
#Used Dr. Smay Stem hw5a, and class lecture resources to build off of.
# region imports
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


# endregion

# region functions
def ff(Re, rr, CBEQN=False):
    """
    Calculate the Darcy-Weisbach friction factor for given flow conditions.

    The friction factor (f) is used in the Darcy-Weisbach equation to compute head loss in a pipe:
        h_f = f * (L/D) * (V^2 / 2g)
    where:
        - L = pipe length
        - D = pipe diameter
        - V = velocity of the fluid
        - g = acceleration due to gravity

    The function determines f based on the flow regime:
    - **Laminar flow (Re < 2000)**: Uses the analytical formula:
        f = 64 / Re
    - **Turbulent flow (Re > 4000)**: Solves the implicit **Colebrook equation**:
        1 / sqrt(f) = -2.0 * log10((ϵ/d) / 3.7 + 2.51 / (Re * sqrt(f)))
      Since this equation cannot be solved explicitly, numerical methods are required.
    - **Transition flow (2000 < Re < 4000)**: Not explicitly defined due to flow instability.

    :param Re: Reynolds number (dimensionless, characterizing flow type).
    :param rr: Relative roughness (pipe wall roughness / diameter ratio).
    :param CBEQN: Boolean flag indicating whether to use the Colebrook equation (True for turbulent flow).
    :return: The computed friction factor (dimensionless).
    """
    if CBEQN:
        # Initial guess for f to improve numerical stability in fsolve.
        initial_guess = 0.25 / (np.log10(rr / 3.7 + 5.74 / Re ** 0.9)) ** 2

        # Colebrook equation, implicit form, requires numerical solving.
        cb = lambda f: 1 / np.sqrt(f) + 2.0 * np.log10(rr / 3.7 + 2.51 / (Re * np.sqrt(f)))

        # Solve for f using fsolve.
        f_solution, _, flag, _ = fsolve(cb, initial_guess, full_output=True)

        # Check for convergence issues.
        if flag != 1:
            print(f"Warning: fsolve did not converge for Re={Re}, rr={rr}. Initial guess: {initial_guess}.")

        return f_solution[0]
    else:
        return 64 / Re


def plotMoody(plotPoint=False, Re=0, rr=0):
    """
    Generate and display the Moody Diagram, showing the friction factor versus Reynolds number.

    The Moody diagram illustrates how the friction factor varies across different flow regimes:
    - **Laminar Flow** (f = 64/Re) is shown as a solid black line.
    - **Transition Flow** (2000 < Re < 4000) is shown as a dashed black line.
    - **Turbulent Flow** is computed numerically using the **Colebrook equation** and plotted for various pipe roughness values.

    :param plotPoint: Boolean flag to indicate whether to highlight a specific point.
    :param Re: Reynolds number of the point to plot (if plotPoint=True).
    :param rr: Relative roughness for the point to plot (if plotPoint=True).
    """
    # Step 1: Define Reynolds number ranges for different flow regimes.
    ReValsCB = np.logspace(np.log10(4000), 8, 100)  # Turbulent range
    ReValsL = np.logspace(np.log10(600.0), np.log10(2000.0), 100)  # Laminar range
    ReValsTrans = np.logspace(np.log10(2000), np.log10(4000), 100)  # Transition range

    # Step 2: Define roughness values (ϵ/d) for turbulent flow curves.
    rrVals = np.array([0, 1E-6, 5E-6, 1E-5, 5E-5, 1E-4, 2E-4, 4E-4, 6E-4, 8E-4,
                       1E-3, 2E-3, 4E-3, 6E-3, 8E-3, 1E-2, 1.5E-2, 2E-2, 3E-2, 4E-2, 5E-2])

    # Step 3: Compute friction factor in laminar and transition regions.
    ffLam = np.array([64 / Re for Re in ReValsL])
    ffTrans = np.array([64 / Re for Re in ReValsTrans])

    # Step 4: Compute friction factor values for turbulent range using Colebrook equation.
    ffCB = np.array([[ff(Re, relRough, CBEQN=True) for Re in ReValsCB] for relRough in rrVals])

    # Step 5: Plot the Moody diagram.
    plt.loglog(ReValsL, ffLam, 'k', label='Laminar Flow')  # Solid black line
    plt.loglog(ReValsTrans, ffTrans, 'k--', label='Transition Flow')  # Dashed black line

    for nRelR, relRough in enumerate(rrVals):
        plt.loglog(ReValsCB, ffCB[nRelR], 'k')  # Black lines for turbulent flow curves
        plt.text(ReValsCB[-1], ffCB[nRelR][-1], f'rr={rrVals[nRelR]}', fontsize=10, va='center', ha='left',
                 color='black')

    # Axis settings.
    plt.xlim([600, 1E8])
    plt.ylim([0.008, 0.10])
    plt.xlabel(r"Reynolds number Re", fontsize=16)
    plt.ylabel(r"Friction factor f", fontsize=16)
    plt.text(2.5E8, 0.02, r"Relative roughness $\frac{\epsilon}{d}$", rotation=90, fontsize=16)
    plt.legend()

    # Grid and axis formatting.
    ax = plt.gca()
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=12)
    ax.tick_params(axis='both', grid_linewidth=1, grid_linestyle='solid', grid_alpha=0.5)
    ax.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
    plt.grid(which='both')
    plt.title("Moody Diagram", fontsize=16)

    # Highlight a specific point if requested.
    if plotPoint:
        f = ff(Re, rr, Re > 4000)
        plt.plot(Re, f, 'ro', markersize=12, label=f'Re={Re}, rr={rr}')
        plt.legend()

    plt.show()


def main():
    """
    Main function to generate and display the Moody diagram.
    """
    plotMoody()


# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion
