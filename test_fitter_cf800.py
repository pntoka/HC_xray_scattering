#!/usr/bin/env python3
"""
Test script for IOBSFitter using CF800 experimental data and predicted parameters.
Runs 10 iterations of fitting to evaluate performance.
"""

import numpy as np
import json
import os
import pickle
from scipy.interpolate import interp1d
from iobs_ngc import IOBSFitter, IOBSParameters

# ============================================================================
# CONFIGURATION - Choose which test to run
# ============================================================================
USE_BASE_CASE = True  # Set to True for base case, False for CF800 experimental data

# Paths for base case test
BASE_CASE_PATH = '/home/ptoka/projects/xrd/data/calc_result_base_case_4.pkl'
START_PARAMS_RANDOM_PATH = '/home/ptoka/projects/iobs_model_dev/data/resnet_mlp_pred/start_params_random.json'

# Paths for CF800 experimental test
SPECTRA_PATH = '/home/ptoka/projects/xrd/data/test_data'
PRED_PARAMS_PATH = '/home/ptoka/projects/iobs_model_dev/data/resnet_mlp_pred/CF800_model_pred_params.json'
EXP_DATA_FILE = 'CF800.xy'


def load_base_case(path):
    """Load base case data from pickle file"""
    with open(path, 'rb') as f:
        base_data = pickle.load(f)
    return base_data['test_s'], base_data['calc_result']


def open_xy(path, file_name):
    """Load XY data file"""
    with open(os.path.join(path, file_name), 'r') as f:
        data = f.readlines()
        data_clean = [line.replace('\t', ' ').split(' ') for line in data]
        data_array = np.array(data_clean).astype(float)
        sorted_data = data_array[data_array[:, 0].argsort()]
    return sorted_data


def two_theta_to_s(x, wavelength=1.5418):
    """Convert 2-theta to scattering vector s"""
    return 2/wavelength * np.sin(x/2*np.pi/180)


def interp_spectra(path, file_name):
    """Load and interpolate experimental spectra to 1000 points"""
    data = open_xy(path, file_name)
    x = two_theta_to_s(data[:, 0])
    test_s = np.linspace(x.min(), x.max(), 1000)
    f_interp = interp1d(x, data[:, 1])
    y_interp = f_interp(test_s)
    return test_s, y_interp


def main():
    print("="*70)
    if USE_BASE_CASE:
        print("IOBSFitter Test - Base Case 4 with Random Start Params")
    else:
        print("IOBSFitter Test - CF800 Experimental Data")
    print("="*70)
    
    # Load data based on configuration
    print("\n[1/4] Loading data...")
    if USE_BASE_CASE:
        s_values, intensities = load_base_case(BASE_CASE_PATH)
        print(f"  Source: Base case pickle file")
    else:
        s_values, intensities = interp_spectra(SPECTRA_PATH, EXP_DATA_FILE)
        print(f"  Source: CF800 experimental data")
    
    print(f"  Data points: {len(s_values)}")
    print(f"  s range: [{s_values.min():.4f}, {s_values.max():.4f}]")
    print(f"  Intensity range: [{intensities.min():.2f}, {intensities.max():.2f}]")
    
    # Load initial parameters based on configuration
    print("\n[2/4] Loading initial parameters...")
    if USE_BASE_CASE:
        params_path = START_PARAMS_RANDOM_PATH
    else:
        params_path = PRED_PARAMS_PATH
    
    with open(params_path, 'r') as f:
        params_dict = json.load(f)
    initial_params = IOBSParameters(params_dict)
    print(f"  Parameters loaded from: {params_path}")
    print(f"  Initial mu: {initial_params.mu:.4f}")
    print(f"  Initial beta: {initial_params.beta:.4f}")
    print(f"  Initial a3: {initial_params.a3:.4f}")
    
    # Create fitter
    print("\n[3/4] Configuring fitter...")
    fitter = IOBSFitter()
    fitter.set_max_iterations(10)  # 10 iterations as requested
    fitter.set_num_threads(8)
    fitter.set_verbose(True)
    fitter.set_function_tolerance(1e-15)
    fitter.set_parameter_tolerance(1e-15)
    fitter.set_gradient_tolerance(1e-15)
    print("  Max iterations: 10")
    print("  Threads: 8")
    print("  Verbose: True")
    
    # Run fitting
    print("\n[4/4] Running fit...")
    print("-"*70)
    result = fitter.fit(s_values.tolist(), intensities.tolist(), initial_params)
    print("-"*70)
    
    # Display results
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    print(f"Success: {result.success}")
    print(f"Convergence: {result.convergence_message}")
    print(f"\nPerformance:")
    print(f"  Iterations: {result.num_iterations}")
    print(f"  Function evaluations: {result.num_function_evaluations}")
    print(f"  Execution time: {result.execution_time_seconds:.2f} seconds")
    print(f"\nGoodness of fit:")
    print(f"  Final cost: {result.final_cost:.6e}")
    print(f"  RSS: {result.rss:.6e}")
    print(f"  R²: {result.r_squared:.6f}")
    
    print(f"\nFinal parameters (unscaled):")
    print(f"  mu:     {result.final_params.mu:.6f}")
    print(f"  beta:   {result.final_params.beta:.6f}")
    print(f"  a3:     {result.final_params.a3:.6f}")
    print(f"  da3:    {result.final_params.da3:.6f}")
    print(f"  sig3:   {result.final_params.sig3:.6f}")
    print(f"  eta:    {result.final_params.eta:.6f}")
    print(f"  alpha:  {result.final_params.alpha:.6f}")
    print(f"  sig1:   {result.final_params.sig1:.6f}")
    print(f"  lcc:    {result.final_params.lcc:.6f}")
    print(f"  q:      {result.final_params.q:.6f}")
    print(f"  g:      {result.final_params.g:.6f}")
    print(f"  u3:     {result.final_params.u3:.6f}")
    print(f"  const1: {result.final_params.const1:.6f}")
    print(f"  const2: {result.final_params.const2:.6f}")
    print(f"  k:      {result.final_params.k:.6f}")
    
    print(f"\nFinal parameters (scaled [0,1]):")
    param_names = ['mu', 'beta', 'a3', 'da3', 'sig3', 'eta', 'alpha', 
                   'sig1', 'lcc', 'q', 'g', 'u3', 'const1', 'const2', 'k']
    for i, name in enumerate(param_names):
        print(f"  {name:7s}: {result.final_params_scaled[i]:.6f}")
    
    print("\n" + "="*70)
    print("Test completed successfully!")
    print("="*70)


if __name__ == "__main__":
    main()
