'''Script to benchmark the execution time of iobs_ngc calculation.'''

import os
import numpy as np
from iobs_ngc import IOBSParameters, IOBSCalculator
import json
import time


def main():
    # Load example fit data
    fit_file = os.path.join(os.path.dirname(__file__), 'example_fit.xy')
    fit_data = np.loadtxt(fit_file)
    s_values = fit_data[:, 0]

    params_file = os.path.join(os.path.dirname(__file__), 'example_params.json')
    with open(params_file, 'r') as f:
        params_dict = json.load(f)

    # Initialize parameters
    params = IOBSParameters(params_dict)

    # Create IOBS calculator
    iobs_calculator = IOBSCalculator()

    # Warm-up run
    _ = iobs_calculator.calculate(params, s_values)

    # Benchmark the calculation
    num_runs = 50
    times = []
    
    print(f"Running {num_runs} iterations...")
    for i in range(num_runs):
        start_time = time.perf_counter()
        iobs_result = iobs_calculator.calculate(params, s_values)
        end_time = time.perf_counter()
        times.append(end_time - start_time)

    times = np.array(times)
    
    print(f"\n{'='*50}")
    print(f"Benchmark Results ({num_runs} runs)")
    print(f"{'='*50}")
    print(f"Mean execution time:   {times.mean()*1000:.3f} ms")
    print(f"Median execution time: {np.median(times)*1000:.3f} ms")
    print(f"Std deviation:         {times.std()*1000:.3f} ms")
    print(f"Min execution time:    {times.min()*1000:.3f} ms")
    print(f"Max execution time:    {times.max()*1000:.3f} ms")
    print(f"{'='*50}\n")


if __name__ == '__main__':
    main()
