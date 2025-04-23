# Running dadi
# I tested these models in an interactive notebook, then run them in full & in parallel in this script
# Script will iterate over a dictionary of models and parameters, run the models in 10x parallel (running 10 optimisation runs with different start points for each model)
# Each model is started using a small grid space (40) then largere grids are used for final likelihood calculations.
# To-Do - save the models as a separate file and import as a library to enable greater reusability
# Tristan Dennis 2025/04/21


# Set up env
import dadi
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from tabulate import tabulate
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import os
import logging
from datetime import datetime
import pickle
import warnings
warnings.filterwarnings("ignore")

# Create output directory for results
output_dir = f"dadi_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
os.makedirs(output_dir, exist_ok=True)

# Set up logging
log_file = os.path.join(output_dir, "model_fitting.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

# Toggle to enable GPU computing
#dadi.cuda_enabled(True)

# Import our own SFS
sfs_array = np.load('rift_joint.sfs.npy')
print(f"Loaded SFS with shape: {sfs_array.shape}")

# Convert to dadi Spectrum object
# You need to provide the sample sizes
sample_sizes = (sfs_array.shape[0] - 1, sfs_array.shape[1] - 1)
sfs = dadi.Spectrum(sfs_array)

# Mask the (0,0) entry which represents invariant sites
sfs.mask[0, 0] = True

# Models, TODO - save as separate file and import as library for better reuse. 
# Pop1 = West, Pop2 = East.
def iso_model(params, ns, pts):
    """
    Isolation model with no migration
    params = (nu1, nu2, T)
    - nu1: Final size of west population
    - nu2: Final size of east population
    """
    nu1, nu2, T = params

    # Convert pts to integer if it's a list
    if isinstance(pts, list):
        pts = pts[0]

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def constant_migration(params, ns, pts):
    """
    Constant migration model
    params = (nu1, nu2, m12, m21, T)
    - nu1: Final size of west population
    - nu2: Final size of east population
    - m12: Migration from west to east
    - m21: Migration from east to wast
    - T: Time since split
    """
    nu1, nu2, m12, m21, T = params

    # Convert pts to integer if it's a list
    if isinstance(pts, list):
        pts = pts[0]

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def founder_effect(params, ns, pts):
    """
    East population founded from West with potential founder effect.

    params = [nuW, nuE, T]
    - nuW: Final size of west population
    - nuE: Final size of east population (founder effect if nuE < 1)
    - T: Time since split
    """
    nuW, nuE, T = params

    # Convert pts to integer if it's a list
    if isinstance(pts, list):
        pts = pts[0]

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuW, nu2=nuE)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def split_with_selection(params, ns, pts):
    """
    Simple split model with selection in one or both populations.

    params = [nu1, nu2, T, gamma1, gamma2]
    - nu1: Size of population 1 (West)
    - nu2: Size of population 2 (East)
    - T: Time since split
    - gamma1: Selection coefficient in population 1 (West)
    - gamma2: Selection coefficient in population 2 (East)
    """
    nu1, nu2, T, gamma1, gamma2 = params

    # Convert pts to integer if it's a list
    if isinstance(pts, list):
        pts = pts[0]

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Integrate with selection
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nu1, nu2=nu2, gamma1=gamma1, gamma2=gamma2)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def founder_with_migration(params, ns, pts):
    """
    East population founded from West with migration after split.

    params = [nuW, nuE, m12, m21, T]
    - nuW: Final size of west population
    - nuE: Final size of east population (founder effect if nuE < 1)
    - m12: Migration from east to west
    - m21: Migration from west to east
    - T: Time since split
    """
    nuW, nuE, m12, m21, T = params

    # Convert pts to integer if it's a list
    if isinstance(pts, list):
        pts = pts[0]

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuW, nu2=nuE, m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def split_asymmetric_size_change2(params, ns, pts):
    """
    Split with asymmetric size change where pop2 contracts
    params = (nu1, nu2_initial, nu2_final, T1, T2, m12, m21)

    nu1 = size of pop1 (constant)
    nu2_initial = initial size of pop2 after split
    nu2_final = final size of pop2 after contraction
    T1 = time between split and start of size change
    T2 = time of size change to present
    m12, m21 = migration rates
    """
    nu1, nu2_initial, nu2_final, T1, T2, m12, m21 = params

    # Convert pts to integer if it's a list
    if isinstance(pts, list):
        pts = pts[0]

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Phase 1: Initial split with migration
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1, nu2_initial, m12=m12, m21=m21)

    # Phase 2: Pop2 changes size
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1, nu2_final, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def two_phase_size_change(params, ns, pts):
    """
    Both populations experience size changes but with different magnitudes
    params = (nu1_initial, nu1_final, nu2_initial, nu2_final, T1, T2, m12, m21)

    nu1_initial = initial size of pop1 after split
    nu1_final = final size of pop1 after size change
    nu2_initial = initial size of pop2 after split
    nu2_final = final size of pop2 after size change
    T1 = time between split and start of size change
    T2 = time of size change to present
    m12, m21 = migration rates
    """
    nu1_initial, nu1_final, nu2_initial, nu2_final, T1, T2, m12, m21 = params

    # Convert pts to integer if it's a list
    if isinstance(pts, list):
        pts = pts[0]

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Phase 1: Initial split with migration
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1_initial, nu2_initial, m12=m12, m21=m21)

    # Phase 2: Both populations change size
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1_final, nu2_final, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def compare_models_aic(model_results, output_dir="."):
    """
    Compare models using AIC
    model_results: list of tuples (model_name, log_likelihood, num_params)
    """
    # Remove any models with None log-likelihood (failed optimization)
    valid_results = [(name, ll, params) for name, ll, params in model_results if ll != -float('inf')]

    # Check if we have any valid models
    if not valid_results:
        msg = "No valid models to compare!"
        print(msg)
        logging.warning(msg)
        return []

    results = []

    # Calculate AIC for each model
    for model_name, ll, num_params in valid_results:
        aic = 2 * num_params - 2 * ll
        results.append((model_name, ll, num_params, aic))

    # Sort by AIC (lower is better)
    results.sort(key=lambda x: x[3])

    # Calculate delta AIC and AIC weights
    min_aic = results[0][3]
    total = 0

    for i in range(len(results)):
        model_name, ll, num_params, aic = results[i]
        delta_aic = aic - min_aic
        results[i] = (model_name, ll, num_params, aic, delta_aic)
        total += np.exp(-0.5 * delta_aic)

    # Add AIC weights
    table = []
    for i in range(len(results)):
        model_name, ll, num_params, aic, delta_aic = results[i]
        weight = np.exp(-0.5 * delta_aic) / total
        table.append([model_name, ll, num_params, aic, delta_aic, weight])

    # Print and log results table
    headers = ["Model", "Log-likelihood", "Parameters", "AIC", "ΔAIC", "AIC weight"]
    table_str = tabulate(table, headers=headers, floatfmt=".4f")
    print(table_str)
    logging.info("Model comparison results:")
    logging.info("\n" + table_str)

    # Save table to file
    with open(os.path.join(output_dir, "model_comparison.txt"), "w") as f:
        f.write(table_str)

    return results

def assess_model_fit(model_func, params, data, pts_l, model_name=None, output_dir="."):
    """Check goodness-of-fit for the best model"""
    if params is None:
        logging.warning("No valid parameters to assess fit!")
        return None

    # Generate the expected SFS
    model_sfs = model_func(params, data.sample_sizes, pts_l[0])

    # Calculate residuals
    resid = data - model_sfs

    # Plot comparison
    fig = plt.figure(figsize=(16, 6))

    # Data
    ax1 = fig.add_subplot(1, 3, 1)
    dadi.Plotting.plot_single_2d_sfs(data, vmin=1, ax=ax1)
    ax1.set_title('Observed SFS')

    # Model
    ax2 = fig.add_subplot(1, 3, 2, sharex=ax1, sharey=ax1)
    dadi.Plotting.plot_single_2d_sfs(model_sfs, vmin=1, ax=ax2)
    ax2.set_title('Model SFS')

    # Residuals
    ax3 = fig.add_subplot(1, 3, 3, sharex=ax1, sharey=ax1)
    dadi.Plotting.plot_2d_resid(resid, resid_range=3, ax=ax3)
    ax3.set_title('Residuals')

    plt.tight_layout()

    # Save figure
    if model_name:
        fig_path = os.path.join(output_dir, f"{model_name}_fit.png")
        plt.savefig(fig_path, dpi=300)
        logging.info(f"Saved model fit plot to {fig_path}")

    # Compute Anscombe residuals
    resid = dadi.Inference.Anscombe_Poisson_residual(model_sfs, data)

    # Log summary statistics
    ll = dadi.Inference.ll_multinom(model_sfs, data)
    theta = dadi.Inference.optimal_sfs_scaling(model_sfs, data)
    ss_resid = np.sum(resid**2)
    avg_resid = np.mean(np.abs(resid))

    logging.info(f"Maximum log-likelihood: {ll}")
    logging.info(f"Theta: {theta}")
    logging.info(f"Sum of squared residuals: {ss_resid}")
    logging.info(f"Average absolute residual: {avg_resid}")

    print(f"Maximum log-likelihood: {ll}")
    print(f"Theta: {theta}")
    print(f"Sum of squared residuals: {ss_resid}")
    print(f"Average absolute residual: {avg_resid}")

    return resid

def run_single_optimization(args):
    """Single optimization run for multiprocessing"""
    model_name, model_func, p0, data, pts_l, lower_bound, upper_bound, run_id = args

    try:
        popt = dadi.Inference.optimize_log(p0, data, model_func, pts_l[0],
                                        lower_bound=lower_bound,
                                        upper_bound=upper_bound,
                                        verbose=False, maxiter=100)
        model = model_func(popt, data.sample_sizes, pts_l[0])
        ll = dadi.Inference.ll_multinom(model, data)

        return model_name, run_id, popt, ll, True
    except Exception as e:
        return model_name, run_id, None, -float('inf'), False

def run_model_optimization(model_name, model_func, param_names, lower_bounds, upper_bounds, data, pts_l, nstarts, output_dir):
    """Run optimization for a single model with parallel runs"""
    msg = f"\n{'-'*50}\nRunning model: {model_name}\n{'-'*50}"
    print(msg)
    logging.info(msg)

    # Prepare all optimization tasks for this model
    all_tasks = []

    for i in range(nstarts):
        # Generate random starting parameters between bounds
        p0 = [np.random.uniform(lower_bounds[j], upper_bounds[j]) for j in range(len(lower_bounds))]
        all_tasks.append((model_name, model_func, p0, data, pts_l, lower_bounds, upper_bounds, i))

    # Determine number of processes
    num_cpus = multiprocessing.cpu_count()
    max_workers = min(40, num_cpus - 1)  # Use up to 40 threads but leave 1 CPU free
    max_workers = max(1, max_workers)  # Ensure at least 1 worker

    msg = f"Running {nstarts} optimization runs using {max_workers} parallel processes..."
    print(msg)
    logging.info(msg)

    # Run optimizations in parallel for this model
    model_runs = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_single_optimization, task) for task in all_tasks]

        for future in futures:
            try:
                run_model_name, run_id, popt, ll, success = future.result()
                if success:
                    model_runs.append((popt, ll))
                    msg = f"{model_name} - Run {run_id+1}: Log-likelihood = {ll:.4f}"
                    print(msg)
                    logging.info(msg)
                else:
                    msg = f"{model_name} - Run {run_id+1} failed"
                    print(msg)
                    logging.warning(msg)
            except Exception as e:
                msg = f"Task failed with error: {e}"
                print(msg)
                logging.error(msg)

    # Save all runs for this model
    all_runs_file = os.path.join(output_dir, f"{model_name}_all_runs.txt")
    with open(all_runs_file, 'w') as f:
        f.write("Run\tLogLikelihood\tParameters\n")
        for i, (params, ll) in enumerate(model_runs):
            f.write(f"{i+1}\t{ll:.4f}\t{dict(zip(param_names, params))}\n")

    # Process results for this model
    if model_runs:
        # Find best result for this model
        best_run = max(model_runs, key=lambda x: x[1])
        best_params, best_ll = best_run

        param_dict = dict(zip(param_names, best_params))

        msg = f"\n{model_name} - Best log-likelihood: {best_ll:.4f}"
        print(msg)
        logging.info(msg)

        msg = f"{model_name} - Best parameters: {param_dict}"
        print(msg)
        logging.info(msg)

        # Save best parameters to file
        with open(os.path.join(output_dir, f"{model_name}_best_params.txt"), "w") as f:
            f.write(f"Log-likelihood: {best_ll:.4f}\n")
            for param, value in param_dict.items():
                f.write(f"{param}: {value:.6f}\n")

        # Assess fit for this model
        msg = f"\nAssessing fit of {model_name}:"
        print(msg)
        logging.info(msg)
        assess_model_fit(model_func, best_params, data, pts_l, model_name, output_dir)

        # Save model SFS
        model_sfs = model_func(best_params, data.sample_sizes, pts_l[0])
        output_path = os.path.join(output_dir, f"{model_name}_sfs.pkl")
        with open(output_path, 'wb') as f:
            pickle.dump(model_sfs, f)

        # Create a comprehensive results dictionary
        full_model_results = {
            "model_name": model_name,
            "model_func": model_func,
            "best_params": best_params,
            "param_dict": param_dict,
            "log_likelihood": best_ll,
            "model_sfs": model_sfs,
            "observed_sfs": data,
            "residuals": data - model_sfs,
            "anscombe_residuals": dadi.Inference.Anscombe_Poisson_residual(model_sfs, data),
            "theta": dadi.Inference.optimal_sfs_scaling(model_sfs, data)
        }

        # Save the complete model results object
        model_results_path = os.path.join(output_dir, f"{model_name}_full_results.pkl")
        with open(model_results_path, 'wb') as f:
            pickle.dump(full_model_results, f)
        logging.info(f"Saved full model results to {model_results_path}")

        return {
            "params": param_dict,
            "ll": best_ll
        }
    else:
        msg = f"\n{model_name} - All runs failed"
        print(msg)
        logging.warning(msg)
        return {
            "params": None,
            "ll": -float('inf')
        }

def sequential_model_comparison(data, output_dir="."):
    """Run models sequentially, with parallel optimization runs for each model"""
    # Set grid points - using smaller grid for faster processing
    pts_l = [70, 80, 90]  # Larger grid for real data
    #pts_l = [10, 20, 30]  # Smaller grid for testing

    logging.info("Starting sequential model comparison")
    logging.info(f"Data shape: {data.shape}")

    model_definitions = {
        "Isolation": (
            iso_model,
            ["nu1", "nu2", "T"],
            [0.1, 0.1, 0.1],
            [10, 10, 5]
        ),
        "Constant_Migration": (
            constant_migration,
            ["nu1", "nu2", "m12", "m21", "T"],
            [0.1, 0.1, 0.1, 0.1, 0.1],
            [10, 10, 10, 10, 5]
        ),
        "Founder_Effect": (
            founder_effect,
            ["nuW", "nuE", "T"],
            [1.0, 0.3, 0.1],  # Eastern pop smaller (founder effect)
            [10.0, 10.0, 5.0]
        ),
        "Split_With_Selection": (
            split_with_selection,
            ["nu1", "nu2", "T", "gamma1", "gamma2"],
            [1.0, 0.5, 0.1, 1.0, 5.0],  # Starting values: stronger selection in east
            [10.0, 10.0, 5.0, 20.0, 20.0]  # Upper bounds
        ),
        "Founder_With_Migration": (
        founder_with_migration,
        ["nuW", "nuE", "m12", "m21", "T"],
        [1.0, 0.3, 0.5, 2.0, 0.1],  # More migration from west to east
        [10.0, 10.0, 10.0, 10.0, 5.0]
    ),

        "Two_Phase_Size_Change": (
            two_phase_size_change,
            ["nu1_initial", "nu1_final", "nu2_initial", "nu2_final", "T1", "T2", "m12", "m21"],
            [0.1, 0.5, 0.1, 0.05, 0.1, 0.1, 0.1, 0.1],
            [10, 20, 10, 5, 3, 3, 10, 10]
        ),
        "Pop2_Contraction": (
            split_asymmetric_size_change2,
            ["nu1", "nu2_initial", "nu2_final", "T1", "T2", "m12", "m21"],
            [0.1, 0.5, 0.01, 0.1, 0.1, 0.1, 0.1],
            [10, 10, 5, 3, 3, 10, 10]
        ),


    }

    logging.info(f"Testing {len(model_definitions)} demographic models")

    nstarts = 10  # Number of optimization runs per model
    model_results = {}
    compare_results = []

    # Run each model sequentially
    for model_name, (model_func, param_names, lower, upper) in model_definitions.items():
        # Run optimization for this model
        result = run_model_optimization(
            model_name=model_name,
            model_func=model_func,
            param_names=param_names,
            lower_bounds=lower,
            upper_bounds=upper,
            data=data,
            pts_l=pts_l,
            nstarts=nstarts,
            output_dir=output_dir
        )

        # Store results
        model_results[model_name] = result
        compare_results.append((model_name, result["ll"], len(param_names)))

    # Save observed SFS to pickle file
    output_path = os.path.join(output_dir, "observed_sfs.pkl")
    with open(output_path, 'wb') as f:
        pickle.dump(data, f)

    # Compare all models
    print("\nModel Comparison:")
    logging.info("Model Comparison:")
    aic_results = compare_models_aic(compare_results, output_dir)

    # Check if we have any valid models
    if not aic_results:
        msg = "No models fitted successfully!"
        print(msg)
        logging.error(msg)
        return model_results, None

    # Get best model
    best_model_name = aic_results[0][0]  # First model after sorting by AIC

    return model_results, best_model_name

# Run the sequential model comparison
model_results, best_model = sequential_model_comparison(sfs, output_dir)

# Print and log detailed results for best model
if best_model is not None:
    msg = f"\n{'='*50}\nBest fitting model: {best_model}\n{'='*50}"
    print(msg)
    logging.info(msg)

    logging.info("Best model parameters:")
    print("Best model parameters:")
    for param, value in model_results[best_model]["params"].items():
        msg = f"  {param}: {value:.6f}"
        print(msg)
        logging.info(msg)

    # Save final summary
    with open(os.path.join(output_dir, "final_results.txt"), "w") as f:
        f.write(f"Best fitting model: {best_model}\n")
        f.write("Parameters:\n")
        for param, value in model_results[best_model]["params"].items():
            f.write(f"  {param}: {value:.6f}\n")
else:
    msg = "\nNo models fitted successfully. Please check your data or model setup."
    print(msg)
    logging.error(msg)

print(f"\nAll results saved to directory: {output_dir}")
logging.info(f"Analysis complete. All results saved to: {output_dir}")
