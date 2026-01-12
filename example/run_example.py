'''Script to run an example of iobs_ngc to verify correct installation.'''

import os
import matplotlib.pyplot as plt
import numpy as np
from iobs_ngc import IOBSParameters, IOBSCalculator
import json


def main():
    # Load example fit data
    fit_file = os.path.join(os.path.dirname(__file__), 'example_fit.xy')
    fit_data = np.loadtxt(fit_file)
    s_values = fit_data[:, 0]
    intensity = fit_data[:, 1]

    params_file = os.path.join(os.path.dirname(__file__), 'example_params.json')
    with open(params_file, 'r') as f:
        params_dict = json.load(f)

    # Initialize parameters
    params = IOBSParameters(params_dict)

    # Create IOBS calculator
    iobs_calculator = IOBSCalculator()

    # calculate IOBS result
    iobs_result = iobs_calculator.calculate(params, s_values)

    # Check if results match the example calculation
    max_diff = np.max(np.abs(iobs_result - intensity))
    rel_diff = np.max(np.abs((iobs_result - intensity) / intensity))

    print(f"Maximum absolute difference: {max_diff:.6e}")
    print(f"Maximum relative difference: {rel_diff:.6e}")

    if np.allclose(iobs_result, intensity, rtol=1e-5, atol=1e-8):
        print("✓ Results match the example calculation!")
    else:
        print("✗ Results differ from the example calculation.")

    # Create overlay plot
    plt.figure(figsize=(10, 6))
    plt.plot(s_values, intensity, 'b-', label='Example calculated result', linewidth=2)
    plt.plot(s_values, iobs_result, 'r--', label='Calculated result', linewidth=1.5, alpha=0.7)
    plt.xlabel('s values')
    plt.ylabel('Intensity')
    plt.title('Comparison of Results')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Save plot to example directory
    plot_file = os.path.join(os.path.dirname(__file__), 'comparison_plot.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {plot_file}")

    plt.show()


if __name__ == '__main__':
    main()

