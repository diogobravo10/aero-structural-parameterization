from nastran.model import nastran_model
import nastran.utility as nastran_utils
import nastran.results as nastran_results
from optimparallel import minimize_parallel
from multiprocessing import Manager
import shutil as sh
import matplotlib.pyplot as plt
import os
import datetime
import math
import numpy as np
# from dataclasses import dataclass
import _prmtr_aero as aero_parameterization
import _prmtr_geomapping as geomapping
import _prmtr_prop as prop_parameterization
import _prmtr_struct as geo_parameterization
import _opt_setup as opt_setup
import _opt_run_utils as _opt_run_utils

msc_nastran_cmd = 'C:/Program Files/MSC.Software/MSC_Nastran/2024.2/bin/nastran.exe'


def stick_bdf_parameterization(sa, di, ar, k_geo, k_mass_wing, k_mass_fuse, k_E, k_G, pid):

    aero_in_bdf = nastran_model.from_file("aero/wing-dlm-coarser.bdf")
    aero_parameterized_bdf = nastran_model.from_file("aero/wing-dlm-coarser.bdf")
    stick_geo_out_bdf = nastran_model.from_file("struct/wing-stick-geo.bdf")
    stick_mass_out_bdf = nastran_model.from_file("mass/wing-stick.bdf")
    fuse_mass_out_bdf = nastran_model.from_file("mass/fuse.bdf")
    stick_spline_out_bdf = nastran_model.from_file("aero/wing-spline-stick.bdf")
    stick_winglet_out_bdf = nastran_model.from_file("aero/nacelle-pylon-winglet-plotting-stick.bdf")
    # winglet_trim_mass_out_bdf = nastran_model.from_file("mass/winglet-trim-m=0.5kg-stick.bdf")
    winglet_trim_mass_out_bdf = nastran_model.from_file("mass/wingtip-trim-m=0.1kg-stick.bdf")
    stick_prop_parameterized_bdf = nastran_model.from_file("struct/wing-stick-prop-opt.bdf")

    aero_parameterization.sweep_angle(aero_parameterized_bdf, sa, n_wc = 0.25)
    aero_parameterization.dihedral_angle(aero_parameterized_bdf, di, n_wc = 0.25)
    aero_parameterization.aspect_ratio(aero_parameterized_bdf, ar, n_wc = 0.25)
    geomapping.mapping(aero_in_bdf, aero_parameterized_bdf, stick_geo_out_bdf, stick_mass_out_bdf, stick_spline_out_bdf, stick_winglet_out_bdf, winglet_trim_mass_out_bdf)
    prop_parameterization.parameterization(stick_prop_parameterized_bdf, k_E, k_G)
    geo_parameterization.geo_parameterization(stick_geo_out_bdf, stick_mass_out_bdf, fuse_mass_out_bdf, aero_parameterized_bdf, k_geo, k_mass_wing, k_mass_fuse)

    # os.makedirs(f"prmtr", exist_ok=True)  # Create directory if it doesn't exist

    # Save parameterized aero/prop bdf files and Save mapped struct bdf files
    aero_parameterized_bdf.to_file(f"prmtr-out/tmp-{pid}/wing-dlm-coarser-prmtr.bdf")
    stick_geo_out_bdf.to_file(f"prmtr-out/tmp-{pid}/wing-stick-geo-prmtr.bdf")
    stick_mass_out_bdf.to_file(f"prmtr-out/tmp-{pid}/wing-stick-mass-prmtr.bdf")
    fuse_mass_out_bdf.to_file(f"prmtr-out/tmp-{pid}/fuse-mass-prmtr.bdf")
    stick_spline_out_bdf.to_file(f"prmtr-out/tmp-{pid}/wing-spline-prmtr.bdf")
    stick_winglet_out_bdf.to_file(f"prmtr-out/tmp-{pid}/n-p-w-plotting-prmtr.bdf")
    # winglet_trim_mass_out_bdf.to_file(f"prmtr/{p_name}-{p_calc:.2f}/winglet-trim-m=0.5kg-stick-mapped.bdf")
    winglet_trim_mass_out_bdf.to_file(f"prmtr-out/tmp-{pid}/wingtip-trim-prmtr.bdf")
    stick_prop_parameterized_bdf.to_file(f"prmtr-out/tmp-{pid}/wing-stick-prop-prmtr.bdf")


def create_deck_prmtr_stick(pid, read_stick_deck_tmp_pid_filename, write_stick_deck_tmp_pid_filename):
    """ Create a temporary STICK modal deck with the current process id.
    """

    with open(read_stick_deck_tmp_pid_filename, 'r') as file:
        stick_deck = file.read()
        stick_deck = stick_deck.replace(f'{{pid}}', str(pid))

    with open(write_stick_deck_tmp_pid_filename, 'w') as file_out:
        file_out.write(stick_deck)



def obj_fun(design_vars: list[float], velocity_point, callback_boolean, shared, lock):

    pid = os.getpid()   # Get the current process ID
    os.mkdir(f"prmtr-out/tmp-{pid}")  # Create a temporary directory for the current process ID
    
    #### Flutter tmp directories #####
    flutter_hybrid_out_tmp_pid_filepath = f'prmtr-out/tmp-{pid}/flutter-deck-hybrid-p-{pid}'    
    flutter_hybrid_pid_filename = "flutter-deck-hybrid-p-pid.bdf"
    flutter_hybrid_tmp_pid_filename = f"flutter-deck-hybrid-p-{pid}.bdf"

    sa, di, ar, k_geo, k_mass_wing, k_mass_fuse, k_E, k_G = design_vars
    stick_bdf_parameterization(sa, di, ar, k_geo, k_mass_wing, k_mass_fuse, k_E, k_G, pid) # Update the BDF model with the design variables
    create_deck_prmtr_stick(pid, flutter_hybrid_pid_filename, flutter_hybrid_tmp_pid_filename) # Create a temporary STICK modal deck with the current process id.

    nastran_utils.run_nastran_analysis(flutter_hybrid_tmp_pid_filename, 
                                       flutter_hybrid_out_tmp_pid_filepath, 
                                       nastran_cmd=msc_nastran_cmd, 
                                       nastran_args=["old=no", "nlines=9999999", "memory=0.2GB", "memorymaximum=0.2GB"],
                                       verbose=False
                                       )
    flutter_results = nastran_results.read_f06(flutter_hybrid_out_tmp_pid_filepath + ".f06")
    flutter_sorted_points = nastran_results.sort_flutter_points(flutter_results.flutter_points, nastran_results.flutter_modes_sort_mode.velocity, reverse=False)
    critical_flutter_results = nastran_results.nastran_results(flutter_points=nastran_results.find_critical_flutter_points(flutter_results.flutter_points))


    # damping_pitch = flutter_sorted_points[velocity_point].modes[target_mode].damping_ratio # positive damping -> unstable -> flutter
    # maximize_positive_damping = -damping_pitch

    lambda_pitch = flutter_sorted_points[velocity_point].modes[6].eigenvalue # positive damping -> unstable -> flutter
    lambda_1elastic = flutter_sorted_points[velocity_point].modes[7].eigenvalue # positive damping -> unstable -> flutter

    coalescing_modes = abs(lambda_pitch - lambda_1elastic)**2 

    os.remove(flutter_hybrid_tmp_pid_filename)
    sh.rmtree(f'prmtr-out/tmp-{pid}')

    print(f"Coalescing Value: {coalescing_modes}")

    with lock:  # protect read-then-write
        if coalescing_modes < shared.analog:
            shared.analog = coalescing_modes

    if callback_boolean:
        shared.lambda_pitch = lambda_pitch
        shared.lambda_1elastic = lambda_1elastic
        shared.coalescing_modes = coalescing_modes
        flutter_point = critical_flutter_results.flutter_points[0]
        if flutter_point.modes:
            shared.flutter_mode_number = next(iter(flutter_point.modes))

    # print(maximize_positive_damping)
    return coalescing_modes


# Callback function definition
def callback(xk, obj_fun, velocity_point, ax, log_filename, callback_boolean, shared, lock) -> None:

    callback.iteration += 1

    nit = callback.iteration

    shared.flutter_mode_number = None
    callback_boolean = True
    obj_fun_value = obj_fun(xk, velocity_point, callback_boolean, shared, lock)
    callback_boolean = False

    ax.plot(nit, obj_fun_value, 'ro')
    plt.draw()
    plt.pause(0.1)

    time_stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print(f"Current design variable values: {xk}")
    print(f"Current objective function value: {obj_fun_value}")
    print(f"Current iteration: {nit}")
    print(f"Current time: {time_stamp}")

    print("Analog =", shared.analog)

    with open(log_filename, 'a') as file:
        file.write(f"{nit}, {time_stamp}, {obj_fun_value}, {shared.analog}, {shared.lambda_pitch}, {shared.lambda_1elastic}, {shared.coalescing_modes}, {shared.flutter_mode_number}")
        for x in xk:
            file.write(f", {x}")
        file.write("\n")


def summary(log_filename, velocity_point, velocity, design_vars, shared, lock):

    log_filename = "log_velocity_points.csv"

    with open(log_filename, 'a') as file:
        # file.write(f"\n{velocity_point}, {velocity:.2f}, {target_mode}, {damping:.5f}")
        file.write(f"\n{velocity_point}, {velocity:.2f}, {shared.lambda_pitch}, {shared.lambda_1elastic}, {shared.coalescing_modes}, {shared.flutter_mode_number}")
        for i, x in enumerate(design_vars):
            file.write(f", {x:.5f}")

    return

def main():
    
    ###################################### Create Directories ######################################
    os.makedirs('prmtr-out', exist_ok=True)
    # os.makedirs('struct/tmp', exist_ok=True)
        
    ###################################### Inital guess for design variables ######################################
    design_vars_0 = [32.84, 3.36, 15.0, 0, 0, 0, 1, 100]

    lb = [15, 1, 5, -1, -1, -1, 1e-3, 1e-3]
    ub = [45, 15, 20, 1, 1, 1, 1e3, 1e3]
    lbub = list(zip(lb, ub))

    ###################################### Velocities to verify flutter ######################################
    aero_bdf = nastran_model.from_file("aero/flutter-h=500m-coarser.bdf")
    velocities = np.array(aero_bdf.flfact_cards[3].factors)
    # velocities = [10., 16.41026, 22.82051, 29.23077, 35.64103, 42.05128, 48.46154,
    #               54.87179, 61.28205, 67.69231, 74.10256, 80.51282, 86.92308, 93.33333, 99.74359,
    #               106.1538, 112.5641, 118.9744, 125.3846, 131.7949]

    log_summary = "log_velocity_points.csv"
    with open(log_summary, 'w') as file:
        file.write("Velocity Point, Velocity, Pitch Eigenvalue, 1st Elastic Mode Eigenvalue, Coalescing Value, Flutter Mode Number")

    for velocity_point, velocity in enumerate(velocities):

        ###################################### Log file setup ######################################
        log_filename = "log.csv"

        time_stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        with open(log_filename, 'w') as file:
            file.write("Iteration, Time, Objective Function Value, Analog, Pitch Eigenvalue, 1st Elastic Mode Eigenvalue, Coalescing Value, Flutter Mode Number")
            for i, x in enumerate(design_vars_0):
                file.write(f", x{i}")
            file.write("\n")
            
        ###################################### Initial design point X0 ######################################   
        with open(log_filename, 'a') as file:    
            file.write(f"0, {time_stamp}, 0, 0")
            for x in design_vars_0:
                file.write(f", {x}")
            file.write("\n")

        ###################################### Objective function plot ######################################
        plt.ion()
        fig, ax = plt.subplots()
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Objective Function Value')
        ax.set_title(f'Damping of Flutter Mode at {velocity}ms-1')

        callback.iteration = 0
        callback_boolean = False

        with Manager() as mgr:
            shared = mgr.Namespace()
            shared.analog = math.inf   
            shared.lambda_pitch = None
            shared.lambda_1elastic = None
            shared.coalescing_modes = None
            shared.flutter_mode_number = None
            lock = mgr.Lock()

            # stages = [
            #     {"rel_eps": 1e-2, "ftol": 1e-3, "gtol": 1e-3, "maxiter": 800},
            #     {"rel_eps": 2e-3, "ftol": 5e-4, "gtol": 5e-4, "maxiter": 1500},
            #     {"rel_eps": 5e-4, "ftol": 1e-5, "gtol": 1e-5, "maxiter": 3000},
            # ]
            ##################################### Parallel Optimization ######################################
            result = minimize_parallel(fun = obj_fun, 
                                    x0 = design_vars_0, 
                                    args = (velocity_point, callback_boolean, shared, lock), 
                                    bounds = lbub,
                                    options = { "disp": None, "maxcor": 10, "ftol": 1e-03, "gtol": 1e-03, "eps": 1e-2, "maxfun": 15000, "maxiter": 15000, "iprint": -1, "maxls": 20 }, 
                                    callback=lambda x: callback(x, obj_fun, velocity_point, ax, log_filename, callback_boolean, shared, lock), 
                                    parallel={'max_workers': 20})
            
            result = minimize_parallel(fun = obj_fun, 
                                    x0 = result.x, 
                                    args = (velocity_point, callback_boolean, shared, lock), 
                                    bounds = lbub,
                                    options = { "disp": None, "maxcor": 10, "ftol": 1e-05, "gtol": 1e-05, "eps": 1e-5, "maxfun": 15000, "maxiter": 15000, "iprint": -1, "maxls": 20 }, 
                                    callback=lambda x: callback(x, obj_fun, velocity_point, ax, log_filename, callback_boolean, shared, lock), 
                                    parallel={'max_workers': 20})
            
            summary(log_summary, velocity_point, velocity, result.x, shared, lock)

        ###################################### Clear pycache ######################################
        sh.rmtree('__pycache__', ignore_errors=True)
        sh.rmtree('nastran/__pycache__', ignore_errors=True)
        sh.rmtree('nastran/cards/__pycache__', ignore_errors=True)

        print("Optimization completed successfully.")
        os.rename("log.csv", "log-" + f"{velocity_point}" + "-" + f"{velocity:.2f}" + ".csv")
        # _opt_run_utils.update_stick_bdf('struct/wing-stick-prop-init.bdf', result.x, 'struct/wing-stick-prop-opt-' + analysis_title + '.bdf')

        plt.close(fig)

        # plt.show(block = True)


if __name__ == '__main__':
    main()
