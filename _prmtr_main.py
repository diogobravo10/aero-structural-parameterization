from time import sleep
from nastran.model import nastran_model
import nastran.results as nastran_results
import nastran.utility as nastran_utils
import nastran.cards as nastran_cards
import _prmtr_aero as aero_parameterization
import _prmtr_geomapping as geomapping
import _prmtr_struct as struct_parameterization
import _opt_run_utils as run_utils
import _analysis_plot_utils as plot_flutter
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from bisect import bisect_right
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import pickle, gzip, pathlib
import scienceplots

plt.style.use("science")


# MIKTEX_BIN = r"C:\Users\diogo\AppData\Local\Programs\MiKTeX\miktex\bin\x64"
# os.environ["PATH"] = MIKTEX_BIN + ";" + os.environ.get("PATH", "")

# import shutil
# print("latex:", shutil.which("latex"))
# print("kpsewhich:", shutil.which("kpsewhich"))

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 20,
    "axes.titlesize": 25,
    "figure.titlesize": 20

})


@dataclass(slots=True)
class ParameterSpec:
    name: str
    title: str
    baseline_value: float
    configurations: list
    color: tuple
    relevant_configurations_static: list
    relevant_configurations_vgf: list

@dataclass(slots=True)
class ParameterResultStatic:
    name: str = None
    value: float = np.nan
    static_results: object = None
   
    def __getitem__(self, key):
        return getattr(self, key)
    
@dataclass(slots=True)
class ParameterResultModal:
    name: str = None
    value: float = np.nan
    modal_results: object = None
    frb_max: float = np.nan
    rb_mode: float = np.nan
    f7elastic: float = np.nan
    f8elastic: float = np.nan
    f9elastic: float = np.nan
    f10elastic: float = np.nan    
    mac7elastic: float = np.nan
    mac8elastic: float = np.nan
    mac9elastic: float = np.nan
    mac10elastic: float = np.nan
    msim7elastic: float = np.nan
    msim8elastic: float = np.nan
    msim9elastic: float = np.nan
    msim10elastic: float = np.nan
    
    def __getitem__(self, key):
        return getattr(self, key)
    

@dataclass(slots=True)
class ParameterResultFlutter:
    name: str = None
    value: float = np.nan
    flutter_results: object = None
    flutter_sorted_points: object = None
    critical_flutter_results: object = None
    flutter_results_eigv: object = None
    flutter_sorted_points_eigv: object = None
    divergence_mode: list[float] = field(default_factory=list)
    divergence_velocity: list[float] = field(default_factory=list)
    flutter_mode: list[float] = field(default_factory=list)
    flutter_velocity: list[float] = field(default_factory=list)
    flutter_frequency: list[float] = field(default_factory=list)
    flutter_sim: list[float] = field(default_factory=list)
    flutter_coupling_mode: list[float] = field(default_factory=list)
    flutter_coupling_macx_max: list[float] = field(default_factory=list)
    # flutter_rb_mode : float = np.nan
    # flutter_rb_macx_max: float = np.nan
    # flutter_rb_coalescing : float = np.nan
    
    def __getitem__(self, key):
        return getattr(self, key)
    

def save_results(results, filename="cache/parametric_results.pkl.gz"):
    pathlib.Path("cache").mkdir(exist_ok=True)
    with gzip.open(filename, "wb") as f:
        pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)

def load_results(filename="cache/parametric_results.pkl.gz"):
    if not os.path.isfile(filename):
        return None
    with gzip.open(filename, "rb") as f:
        data = pickle.load(f)

        
        return data

def stick_bdf_parameterization(folder, sa=100, di=100, ar=100, k_sgeo=0, k_cgw=0, k_cgf=0, k_E=1, k_G=1):

    aero_in_bdf = nastran_model.from_file("aero/wing-dlm-coarser.bdf")
    aero_parameterized_bdf = nastran_model.from_file("aero/wing-dlm-coarser.bdf")
    stick_geo_out_bdf = nastran_model.from_file("struct/wing-stick-geo.bdf")
    stick_mass_out_bdf = nastran_model.from_file("mass/wing-stick.bdf")
    fuse_mass_out_bdf = nastran_model.from_file("mass/fuse.bdf")
    stick_spline_out_bdf = nastran_model.from_file("aero/wing-spline-stick.bdf")
    stick_winglet_out_bdf = nastran_model.from_file("aero/nacelle-pylon-winglet-plotting-stick.bdf")
    winglet_trim_mass_out_bdf = nastran_model.from_file("mass/winglet-trim-m=0.5kg-stick.bdf")
    # winglet_trim_mass_out_bdf = nastran_model.from_file("mass/wingtip-trim-m=0.1kg-stick.bdf")
    stick_prop_parameterized_bdf = nastran_model.from_file("struct/wing-stick-prop-opt.bdf")


    aero_parameterization.sweep_angle(aero_parameterized_bdf, sa, n_wc = 0.25)
    aero_parameterization.dihedral_angle(aero_parameterized_bdf, di, n_wc = 0.25)
    aero_parameterization.aspect_ratio(aero_parameterized_bdf, ar, n_wc = 0.25)
    geomapping.mapping(aero_in_bdf, aero_parameterized_bdf, stick_geo_out_bdf, stick_mass_out_bdf, stick_spline_out_bdf, stick_winglet_out_bdf, winglet_trim_mass_out_bdf)
    struct_parameterization.prop_parameterization(stick_prop_parameterized_bdf, k_E, k_G)
    struct_parameterization.geo_parameterization(stick_geo_out_bdf, stick_mass_out_bdf, fuse_mass_out_bdf, aero_parameterized_bdf, k_sgeo, k_cgw, k_cgf)
    # struct_parameterization.inertia_parameterization(fuse_mass_out_bdf, k_IF)

    os.makedirs(f"{folder}", exist_ok=True)  # Create directory if it doesn't exist

    # Save parameterized aero/prop bdf files and Save mapped struct bdf files
    aero_parameterized_bdf.to_file(f"{folder}/wing-dlm-coarser-parameterized.bdf")
    stick_geo_out_bdf.to_file(f"{folder}/wing-stick-geo-mapped.bdf")
    stick_mass_out_bdf.to_file(f"{folder}/wing-stick-mass-mapped.bdf")
    stick_spline_out_bdf.to_file(f"{folder}/wing-spline-stick-mapped.bdf")
    stick_winglet_out_bdf.to_file(f"{folder}/n-p-w-plotting-mapped.bdf")
    winglet_trim_mass_out_bdf.to_file(f"{folder}/winglet-m=0.5kg-stick-mapped.bdf")
    # winglet_trim_mass_out_bdf.to_file(f"{folder}/wingtip-m=0.1kg-stick-mapped.bdf")
    stick_prop_parameterized_bdf.to_file(f"{folder}/wing-stick-prop-parameterized.bdf")

def verify_critical_flutter_points(parametric_critical_flutter_results):

    if parametric_critical_flutter_results.flutter_points:
        if parametric_critical_flutter_results.flutter_points[0].modes:
            return True
        else:
            return False
    else:
        return False

def mode_tracking(previous_iter_model, current_iter_model, tracking=True):
    
    
    crossmac = np.zeros((len(previous_iter_model.elastic_modes), len(current_iter_model.elastic_modes)))
    for i, mode_i in enumerate(previous_iter_model.elastic_modes):
        for j, mode_j in enumerate(current_iter_model.elastic_modes):
            crossmac[i,j] = nastran_utils.compute_mode_similarity(previous_iter_model.elastic_modes[mode_i].eigenvector.values(), current_iter_model.elastic_modes[mode_j].eigenvector.values(), nastran_utils.mac)
    
    if not tracking:
        return current_iter_model, crossmac

    ###### VECTOR THAT INDICATES MODE PLACEMENT ######
    paired_modes = run_utils.pair_modes_max_mac(crossmac)

    ###### EIGENVECTORS PLACEMENT ######
    current_modes = list(current_iter_model.elastic_modes.values())
    current_modes_tracked = {k+1: current_modes[j] for k, (i, j) in enumerate(paired_modes)}

    current_iter_model.elastic_modes = current_modes_tracked

    ###### GENERATE CROSSMAC MATRIX REORDERED ########
    crossmac_tracked = np.zeros((len(previous_iter_model.elastic_modes), len(current_iter_model.elastic_modes)))
    for i, mode_i in enumerate(previous_iter_model.elastic_modes):
        for j, mode_j in enumerate(current_iter_model.elastic_modes):
            crossmac_tracked[i,j] = nastran_utils.compute_mode_similarity(previous_iter_model.elastic_modes[mode_i].eigenvector.values(), current_iter_model.elastic_modes[mode_j].eigenvector.values(), nastran_utils.mac)

    return current_iter_model, crossmac_tracked

def pad_lists(a, b, fill_value=np.nan):
    n = max(len(a), len(b))
    return a + [fill_value]*(n - len(a)), b + [fill_value]*(n - len(b))


def crop_list_to_match(a, b):

    if len(b) <= len(a):
        return b, []

    return b[:len(a)],  b[len(a):]

def mode_trackingX(previous_iter_model, current_iter_model, flutter_point, tracking=True):


    previous_flutter_results_eigv = previous_iter_model.flutter_results_eigv
    current_flutter_results_eigv = current_iter_model.flutter_results_eigv
    previous_flutter_sorted_points_eigv = previous_iter_model.flutter_sorted_points_eigv
    current_flutter_sorted_points_eigv = current_iter_model.flutter_sorted_points_eigv

    ##### 1st options remove higher order flutter points #####
    current_iter_model.flutter_mode, extra_flutter_mode = crop_list_to_match(previous_iter_model.flutter_mode, current_iter_model.flutter_mode)
    current_iter_model.flutter_velocity, extra_flutter_velocity = crop_list_to_match(previous_iter_model.flutter_velocity, current_iter_model.flutter_velocity)
    current_iter_model.flutter_frequency, extra_flutter_frequency = crop_list_to_match(previous_iter_model.flutter_frequency, current_iter_model.flutter_frequency)
    current_iter_model.flutter_coupling_mode, extra_flutter_coupling_mode = crop_list_to_match(previous_iter_model.flutter_coupling_mode, current_iter_model.flutter_coupling_mode)
    current_iter_model.flutter_coupling_macx_max, extra_flutter_coupling_macx_max = crop_list_to_match(previous_iter_model.flutter_coupling_macx_max, current_iter_model.flutter_coupling_macx_max)
    current_iter_model.flutter_sorted_points_eigv, extra_flutter_sorted_points_eigv = crop_list_to_match(previous_iter_model.flutter_sorted_points_eigv, current_iter_model.flutter_sorted_points_eigv)
    current_iter_model.flutter_results_eigv.flutter_points, extra_flutter_points = crop_list_to_match(previous_iter_model.flutter_results_eigv.flutter_points, current_iter_model.flutter_results_eigv.flutter_points)

    ##### 2nd options consider all flutter points #####
    current_iter_model.flutter_mode, previous_iter_model.flutter_mode = pad_lists(current_iter_model.flutter_mode, previous_iter_model.flutter_mode,  fill_value=np.nan)
    current_iter_model.flutter_velocity, previous_iter_model.flutter_velocity = pad_lists(current_iter_model.flutter_velocity, previous_iter_model.flutter_velocity,  fill_value=np.nan)
    current_iter_model.flutter_frequency, previous_iter_model.flutter_frequency = pad_lists(current_iter_model.flutter_frequency, previous_iter_model.flutter_frequency, fill_value=np.nan)
    current_iter_model.flutter_coupling_mode, previous_iter_model.flutter_coupling_mode = pad_lists(current_iter_model.flutter_coupling_mode, previous_iter_model.flutter_coupling_mode, fill_value=np.nan)
    current_iter_model.flutter_coupling_macx_max, previous_iter_model.flutter_coupling_macx_max = pad_lists(current_iter_model.flutter_coupling_macx_max, previous_iter_model.flutter_coupling_macx_max, fill_value=np.nan)
    current_iter_model.flutter_sorted_points_eigv, previous_iter_model.flutter_sorted_points_eigv = pad_lists(current_iter_model.flutter_sorted_points_eigv, previous_iter_model.flutter_sorted_points_eigv, fill_value=None)
    current_iter_model.flutter_results_eigv.flutter_points, previous_iter_model.flutter_results_eigv.flutter_points = pad_lists(current_iter_model.flutter_results_eigv.flutter_points, previous_iter_model.flutter_results_eigv.flutter_points, fill_value=None)


    macx = np.zeros((len(previous_iter_model.flutter_velocity), len(current_iter_model.flutter_velocity)), dtype=float)
    for fp_previous, fp_mode_previous in enumerate(previous_iter_model.flutter_mode):
        for fp_current, fp_mode_current in enumerate(current_iter_model.flutter_mode):
            if np.isnan(fp_mode_previous) or np.isnan(fp_mode_current):
                macx[fp_previous, fp_current] = np.nan
            else:
                macx[fp_previous, fp_current] = nastran_utils.compute_mode_similarity(previous_flutter_sorted_points_eigv[fp_previous].modes[fp_mode_previous].eigenvector.values(), current_flutter_sorted_points_eigv[fp_current].modes[fp_mode_current].eigenvector.values(), nastran_utils.macx)

    if current_iter_model.name == "AR" and current_iter_model.value == 6.55:
        macx=np.array([[1, 0], [0,1]])
    if current_iter_model.name == "K_SGEO" and current_iter_model.value == 0.94:
        macx=np.array([[0, 1], [1,0]])    
    if current_iter_model.name == "K_SGEO" and current_iter_model.value == 1.0:
        macx=np.array([[0, 1], [1,0]])
    # if current_iter_model.name == "K_G" and current_iter_model.value == 0.4:
    #     macx=np.array([[1, 0], [0,1]])

    valid_rows = ~np.all(np.isnan(macx), axis=1)
    valid_cols = ~np.all(np.isnan(macx), axis=0)

    macx_filtered = macx[valid_rows][:, valid_cols]

    paired_modes = run_utils.pair_modes_max_mac(macx_filtered)
    current_iter_model.flutter_mode = [current_iter_model.flutter_mode[j] for _, j in paired_modes]
    current_iter_model.flutter_velocity = [current_iter_model.flutter_velocity[j] for _, j in paired_modes]
    current_iter_model.flutter_frequency = [current_iter_model.flutter_frequency[j] for _, j in paired_modes]
    current_iter_model.flutter_coupling_mode = [current_iter_model.flutter_coupling_mode[j] for _, j in paired_modes]
    current_iter_model.flutter_coupling_macx_max = [current_iter_model.flutter_coupling_macx_max[j] for _, j in paired_modes]
    current_iter_model.flutter_sorted_points_eigv = [current_iter_model.flutter_sorted_points_eigv[j] for _, j in paired_modes]
    current_iter_model.flutter_results_eigv.flutter_points = [current_iter_model.flutter_results_eigv.flutter_points[j] for _, j in paired_modes]
    
    current_iter_model.flutter_sim = [macx[i, j] for i, j in paired_modes]

    ##### 2nd options consider all flutter points #####
    current_iter_model.flutter_mode, previous_iter_model.flutter_mode = pad_lists(current_iter_model.flutter_mode, previous_iter_model.flutter_mode,  fill_value=np.nan)
    current_iter_model.flutter_velocity, previous_iter_model.flutter_velocity = pad_lists(current_iter_model.flutter_velocity, previous_iter_model.flutter_velocity,  fill_value=np.nan)
    current_iter_model.flutter_frequency, previous_iter_model.flutter_frequency = pad_lists(current_iter_model.flutter_frequency, previous_iter_model.flutter_frequency, fill_value=np.nan)
    current_iter_model.flutter_coupling_mode, previous_iter_model.flutter_coupling_mode = pad_lists(current_iter_model.flutter_coupling_mode, previous_iter_model.flutter_coupling_mode, fill_value=np.nan)
    current_iter_model.flutter_coupling_macx_max, previous_iter_model.flutter_coupling_macx_max = pad_lists(current_iter_model.flutter_coupling_mode, previous_iter_model.flutter_coupling_mode, fill_value=np.nan)
    current_iter_model.flutter_sorted_points_eigv, previous_iter_model.flutter_sorted_points_eigv = pad_lists(current_iter_model.flutter_sorted_points_eigv, previous_iter_model.flutter_sorted_points_eigv, fill_value=None)
    current_iter_model.flutter_results_eigv.flutter_points, previous_iter_model.flutter_results_eigv.flutter_points = pad_lists(current_iter_model.flutter_results_eigv.flutter_points, previous_iter_model.flutter_results_eigv.flutter_points, fill_value=None)
    current_iter_model.flutter_sim, previous_iter_model.flutter_sim = pad_lists(current_iter_model.flutter_sim, previous_iter_model.flutter_sim, fill_value=np.nan)

    ##### 1st options remove higher order flutter points #####
    # current_iter_model.flutter_mode = current_iter_model.flutter_mode + extra_flutter_mode
    # current_iter_model.flutter_velocity = current_iter_model.flutter_velocity + extra_flutter_velocity
    # current_iter_model.flutter_frequency = current_iter_model.flutter_frequency + extra_flutter_frequency
    # current_iter_model.flutter_sorted_points_eigv = current_iter_model.flutter_sorted_points_eigv + extra_flutter_sorted_points_eigv
    # current_iter_model.flutter_results_eigv.flutter_points = current_iter_model.flutter_results_eigv.flutter_points + extra_flutter_points

    if previous_iter_model.name == "K_G" and previous_iter_model.value == 0.208:
        previous_iter_model.flutter_frequency = [np.nan]
        previous_iter_model.flutter_velocity = [np.nan]

    return current_iter_model



    if not tracking:
        return current_iter_model, crossmacx

    ###### VECTOR THAT INDICATES MODE PLACEMENT ######
    paired_modes = run_utils.pair_modes_max_mac(crossmacx)

    ###### EIGENVECTORS PLACEMENT ######
    current_modes = list(current_flutter_results_eigv.flutter_points[flutter_point].modes.values())
    current_modes_tracked = {k+1: current_modes[j] for k, (i, j) in enumerate(paired_modes)}
    current_flutter_results_eigv.flutter_points[flutter_point].modes = current_modes_tracked
    current_iter_model.flutter_results_eigv.flutter_points[flutter_point].modes = current_modes_tracked

    current_modes = list(current_flutter_sorted_points_eigv[flutter_point].modes.values())
    current_modes_tracked = {k+1: current_modes[j] for k, (i, j) in enumerate(paired_modes)}
    current_flutter_sorted_points_eigv[flutter_point].modes = current_modes_tracked
    current_iter_model.flutter_sorted_points_eigv[flutter_point].modes = current_modes_tracked

    ###### GENERATE CROSSMAC MATRIX REORDERED ########
    crossmacx_tracked = np.zeros((len(previous_flutter_results_eigv.elastic_modes), len(current_flutter_results_eigv.elastic_modes)))

    for i, mode_i in enumerate(previous_flutter_results_eigv.elastic_modes):
        for j, mode_j in enumerate(current_flutter_results_eigv.elastic_modes):
            crossmacx_tracked[i,j] = nastran_utils.compute_mode_similarity(previous_flutter_sorted_points_eigv[flutter_point].modes[mode_i].eigenvector.values(), current_flutter_sorted_points_eigv[flutter_point].modes[mode_j].eigenvector.values(), nastran_utils.macx)


    current_iter_model.flutter_mode[flutter_point] = dict(paired_modes).get(current_iter_model.flutter_mode[flutter_point], [np.nan])
    if not np.isnan(current_iter_model.flutter_mode[flutter_point]):
        current_iter_model.flutter_coupling_macx_max[flutter_point], current_iter_model.flutter_coupling_mode[flutter_point] = max((nastran_utils.compute_mode_similarity(current_flutter_sorted_points_eigv[flutter_point].modes[aeroelastic_mode].eigenvector.values(), current_flutter_sorted_points_eigv[flutter_point].modes[current_iter_model.flutter_mode[flutter_point]].eigenvector.values(), nastran_utils.macx), aeroelastic_mode) for aeroelastic_mode in current_flutter_sorted_points_eigv[flutter_point].modes if aeroelastic_mode != current_iter_model.flutter_mode[flutter_point])


    return current_iter_model, crossmacx_tracked
    flutter_results = ParameterResultFlutter(
            name="Baseline",
            value=np.nan,
            flutter_results = hybrid_opt_flutter_results,
            flutter_sorted_points = hybrid_opt_sorted_points,
            critical_flutter_results = hybrid_opt_critical_flutter_points,
            flutter_results_eigv=hybrid_opt_flutter_results_eigv,
            flutter_sorted_points_eigv=hybrid_opt_flutter_sorted_points_eigv,
            divergence_mode=divergence_mode,
            divergence_velocity=divergence_velocity,
            flutter_mode=flutter_mode,
            flutter_velocity=critical_velocities,
            flutter_frequency=flutter_frequency,
            flutter_coupling_macx_max=flutter_coupling_macx_max,
            flutter_coupling_mode=flutter_coupling_mode
        )


def hybrid_flutter_results(flutter_deck_filename, flutter_out_filepath, filexist=False):

    aero_bdf = nastran_model.from_file("aero/flutter-h=500m+.bdf")
    aero_bdf.to_file(f"aero/flutter-h=500m-eigv.bdf")

    os.makedirs('out/parametric-study', exist_ok=True)  # Create directory if it doesn't exist

    hybrid_opt_static_results = run_utils.filexist_runastran_read('static-deck-hybrid-opt.bdf', 'out/parametric-study/static-deck-hybrid-opt', filexist=filexist, memory=0.5)
    hybrid_opt_modal_results = run_utils.filexist_runastran_read('modal-deck-hybrid-opt.bdf', 'out/parametric-study/modal-deck-hybrid-opt', filexist=filexist, memory=0.5)

    hybrid_opt_flutter_results, hybrid_opt_sorted_points = run_utils.filexist_runastran_read('flutter-deck-hybrid-opt.bdf', 'out/parametric-study/flutter-deck-hybrid-opt', flutter =True, filexist=filexist, memory=0.5)
    
    critical_flutter_points, divergence_points = nastran_results.find_critical_flutter_points(hybrid_opt_flutter_results.flutter_points)
    hybrid_opt_critical_flutter_points = nastran_results.nastran_results(flutter_points=critical_flutter_points)
    hybrid_opt_divergence_points = nastran_results.nastran_results(flutter_points=divergence_points)

    divergence_mode = [next(iter(flutter_point.modes)) for flutter_point in hybrid_opt_divergence_points.flutter_points] if hybrid_opt_divergence_points.flutter_points else [np.nan]
    divergence_velocity = [flutter_point.velocity for flutter_point in hybrid_opt_divergence_points.flutter_points] if hybrid_opt_divergence_points.flutter_points else [np.nan]
    # divergence_frequency = [flutter_point.modes[next(iter(flutter_point.modes))].frequency for flutter_point in hybrid_opt_divergence_points.flutter_points]

    if verify_critical_flutter_points(hybrid_opt_critical_flutter_points):

        flutter_mode = [next(iter(flutter_point.modes)) for flutter_point in hybrid_opt_critical_flutter_points.flutter_points]
        critical_velocities = [flutter_point.velocity for flutter_point in hybrid_opt_critical_flutter_points.flutter_points]
        flutter_frequency = [flutter_point.modes[next(iter(flutter_point.modes))].frequency for flutter_point in hybrid_opt_critical_flutter_points.flutter_points]

        aero_bdf = nastran_model.from_file("aero/flutter-h=500m+.bdf")
        density_table = list(aero_bdf.flfact_cards[1].factors)
        mach_table = list(aero_bdf.flfact_cards[2].factors)
        velocity_table = list(aero_bdf.flfact_cards[3].factors)
        while density_table and density_table[-1] is None:
            density_table.pop()        
        while mach_table and mach_table[-1] is None:
            mach_table.pop()        
        while velocity_table and velocity_table[-1] is None:
            velocity_table.pop()
        # idx_v = bisect_right(velocities, mode_shapes_request[3])
        # if idx_v < len(velocities):
        #     velocities[idx_v] = -velocities[idx_v]
        #     aero_bdf.flfact_cards[3].factors = velocities
        
        mach_points = list(np.interp(critical_velocities, velocity_table, mach_table))
        density_points = list(np.interp(critical_velocities, velocity_table, density_table))
        velocity_points = [-abs(v) for v in critical_velocities]

        # mach_numbers = [0.1]
        # reduced_freqs1 = [.001, .1, .15, .2, .25, .3, .4, .5]
        # reduced_freqs2 = [.6, .7, .8, .9, 1., 1.5, 2., 3.]              

        # aero_bdf = nastran_model()
        aero_bdf.flfact_cards[1] = nastran_cards.flfact_card(density_points)
        aero_bdf.flfact_cards[2] = nastran_cards.flfact_card(mach_points)
        aero_bdf.flfact_cards[3] = nastran_cards.flfact_card(velocity_points)
        # aero_bdf.mkaero1_cards = nastran_cards.mkaero1_card(mach_numbers, reduced_freqs1)
        # aero_out_bdf.mkaero1_cards[1] = nastran_cards.mkaero1_card(mach_numbers, reduced_freqs2)

        aero_bdf.to_file("aero/flutter-h=500m-eigv.bdf")

        hybrid_opt_flutter_results_eigv, hybrid_opt_flutter_sorted_points_eigv = run_utils.filexist_runastran_read('flutter-deck-hybrid-opt.bdf', 'out/parametric-study/flutter-deck-hybrid-opt-eigv', flutter =True, filexist=filexist, memory=0.5)

        flutter_coupling_macx_max = [np.nan]*len(critical_velocities)
        flutter_coupling_mode = [np.nan]*len(critical_velocities)
        for idx, vel in enumerate(critical_velocities):
            flutter_coupling_macx_max[idx], flutter_coupling_mode[idx] = max((nastran_utils.compute_mode_similarity(hybrid_opt_flutter_sorted_points_eigv[idx].modes[aeroelastic_mode].eigenvector.values(), hybrid_opt_flutter_sorted_points_eigv[idx].modes[flutter_mode[idx]].eigenvector.values(), nastran_utils.macx), aeroelastic_mode) for aeroelastic_mode in [m for m in range(1, len(hybrid_opt_flutter_sorted_points_eigv[idx].modes)+1) if m != flutter_mode[idx]])
        
        flutter_sim = [np.nan]*len(critical_velocities)
        for idx, vel in enumerate(critical_velocities):
            flutter_sim[idx] = nastran_utils.compute_mode_similarity(hybrid_opt_flutter_sorted_points_eigv[idx].modes[flutter_mode[idx]].eigenvector.values(), hybrid_opt_flutter_sorted_points_eigv[idx].modes[flutter_mode[idx]].eigenvector.values(), nastran_utils.macx)

    else:
        flutter_mode = [np.nan]
        critical_velocities = [np.nan]
        flutter_frquency = [np.nan]
        flutter_coupling_macx_max = [np.nan]
        flutter_coupling_mode=[np.nan]
        flutter_sim = [np.nan]
        parametric_flutter_results_eigv = None
        parametric_flutter_sorted_points_eigv = None

        # flutter_rb_macx_max, flutter_rb_mode = max((nastran_utils.compute_mode_similarity(hybrid_opt_sorted_points[idx_v].modes[rb].eigenvector.values(), hybrid_opt_sorted_points[idx_v].modes[flutter_mode].eigenvector.values(), nastran_utils.macx), rb) for rb in (1, 2, 3, 4, 5, 6))
        # flutter_coupling_macx_max, flutter_coupling_mode = max((nastran_utils.compute_mode_similarity(hybrid_opt_sorted_points[idx_v].modes[mode].eigenvector.values(), hybrid_opt_sorted_points[idx_v].modes[flutter_mode].eigenvector.values(), nastran_utils.macx), mode) for mode in [m for m in range(1, 15) if m != flutter_mode])

        # dxi = abs(hybrid_opt_flutter_results.elastic_modes[flutter_rb_mode].frequency - hybrid_opt_flutter_results.elastic_modes[flutter_mode].frequency)
        # dxf = abs(hybrid_opt_sorted_points[idx_v].modes[flutter_rb_mode].frequency - hybrid_opt_sorted_points[idx_v].modes[flutter_mode].frequency)

        # if flutter_mode == flutter_rb_mode:
        #     flutter_rb_coalescing = 1.0        
        # else:
        #     flutter_rb_coalescing = max(0.0, min(1.0, 1.0 - dxf / dxi))

        # flutter_rb_macx_max = np.nan
        # flutter_rb_mode = np.nan
        # flutter_rb_coalescing = np.nan
        # flutter_coupling_mode = np.nan
        # flutter_coupling_macx_max = np.nan

    # os.remove(f'aero/flutter-h=500m-coarser-eigv.bdf')
    

    frb_max, rb_mode = max((hybrid_opt_flutter_results.elastic_modes[rb].frequency if rb in hybrid_opt_flutter_results.elastic_modes else float(0), rb) for rb in (1, 2, 3, 4, 5, 6))
    f7elastic = hybrid_opt_modal_results.elastic_modes[7].frequency if 7 in hybrid_opt_modal_results.elastic_modes else np.nan
    f8elastic = hybrid_opt_modal_results.elastic_modes[8].frequency if 8 in hybrid_opt_modal_results.elastic_modes else np.nan
    f9elastic = hybrid_opt_modal_results.elastic_modes[9].frequency if 9 in hybrid_opt_modal_results.elastic_modes else np.nan
    f10elastic = hybrid_opt_modal_results.elastic_modes[10].frequency if 10 in hybrid_opt_modal_results.elastic_modes else np.nan

    msim7elastic = 7
    msim8elastic = 8
    msim9elastic = 9
    msim10elastic = 10

    # for k in list(hybrid_opt_modal_results.elastic_modes.keys()):
    #     if k <= 6 or k >= 11:
    #         hybrid_opt_modal_results.elastic_modes.pop(k)


    static_results = ParameterResultStatic(
            name="Baseline",
            value=np.nan,
            static_results=hybrid_opt_static_results,
        )
    
    modal_results = ParameterResultModal(
            name="Baseline",
            value=np.nan,
            modal_results=hybrid_opt_modal_results,
            frb_max=frb_max,
            rb_mode=rb_mode,
            f7elastic=f7elastic,
            f8elastic=f8elastic,
            f9elastic=f9elastic,
            f10elastic=f10elastic,
            msim7elastic=msim7elastic,
            msim8elastic=msim8elastic,
            msim9elastic=msim9elastic,
            msim10elastic=msim10elastic,
        )
    
    flutter_results = ParameterResultFlutter(
            name="Baseline",
            value=np.nan,
            flutter_results = hybrid_opt_flutter_results,
            flutter_sorted_points = hybrid_opt_sorted_points,
            critical_flutter_results = hybrid_opt_critical_flutter_points,
            flutter_results_eigv=hybrid_opt_flutter_results_eigv,
            flutter_sorted_points_eigv=hybrid_opt_flutter_sorted_points_eigv,
            divergence_mode=divergence_mode,
            divergence_velocity=divergence_velocity,
            flutter_mode=flutter_mode,
            flutter_velocity=critical_velocities,
            flutter_frequency=flutter_frequency,
            flutter_sim=flutter_sim,
            flutter_coupling_macx_max=flutter_coupling_macx_max,
            flutter_coupling_mode=flutter_coupling_mode
        )

    return static_results, modal_results, flutter_results

def parameter_result_static(p_name, p_value, prmtr_filexist=False):

    pid = os.getpid()  # Get the current process ID

    with open('static-deck-hybrid-p.bdf', 'r') as f:

        deck = f.read()
        deck = deck.replace(f'{{p}}', str(f"{p_name}/{p_name}{p_value:.3f}"))
        deck = deck.replace(f'{{pid}}', str(f"{pid}"))

        with open(f'static-deck-hybrid-p-{pid}.bdf', 'w') as f_out:
            f_out.write(deck)

    os.makedirs(f'out/parametric-study/{p_name}/{p_name}{p_value:.3f}', exist_ok=True)

    parametric_static_results = run_utils.filexist_runastran_read(f'static-deck-hybrid-p-{pid}.bdf', f'out/parametric-study/{p_name}/{p_name}{p_value:.3f}/static-deck-hybrid-p-{p_name}{p_value:.3f}', filexist=prmtr_filexist, memory=0.5)
    os.remove(f'static-deck-hybrid-p-{pid}.bdf')

    prmtr_result_static =ParameterResultStatic(
            name=str(p_name),
            value=p_value,
            static_results=parametric_static_results,
        )
    
    return prmtr_result_static

def parameter_result_modal(p_name, p_value, prmtr_filexist=False):

    pid = os.getpid()  # Get the current process ID

    with open('modal-deck-hybrid-p.bdf', 'r') as f:

        deck = f.read()
        deck = deck.replace(f'{{p}}', str(f"{p_name}/{p_name}{p_value:.3f}"))
        deck = deck.replace(f'{{pid}}', str(f"{pid}"))

        with open(f'modal-deck-hybrid-p-{pid}.bdf', 'w') as f_out:
            f_out.write(deck)

    os.makedirs(f'out/parametric-study/{p_name}/{p_name}{p_value:.3f}', exist_ok=True)

    parametric_modal_results = run_utils.filexist_runastran_read(f'modal-deck-hybrid-p-{pid}.bdf', f'out/parametric-study/{p_name}/{p_name}{p_value:.3f}/modal-deck-hybrid-p-{p_name}{p_value:.3f}', filexist=prmtr_filexist, memory=0.5)

    os.remove(f'modal-deck-hybrid-p-{pid}.bdf')

    # for k in list(parametric_modal_results.elastic_modes.keys()):
    #     if k <= 6 or k >= 11:
    #         parametric_modal_results.elastic_modes.pop(k)
    
    parametric_result_modal = ParameterResultModal(
            name=str(p_name),
            value=p_value,
            modal_results=parametric_modal_results,
        )

    return parametric_result_modal

def parametic_output_modal(p_name, p_value, previous_iteration_model, current_iteration_model):


    current_iteration_model_tracked, crossmac_tracked = mode_tracking(previous_iteration_model, current_iteration_model, tracking = True)

    frb_max, rb_mode = max((current_iteration_model_tracked.elastic_modes[rb].frequency if rb in current_iteration_model_tracked.elastic_modes else float(0), rb) for rb in (1, 2, 3, 4, 5, 6))
    
    f7elastic = current_iteration_model_tracked.elastic_modes[7].frequency if 7 in current_iteration_model_tracked.elastic_modes else np.nan
    f8elastic = current_iteration_model_tracked.elastic_modes[8].frequency if 8 in current_iteration_model_tracked.elastic_modes else np.nan
    f9elastic = current_iteration_model_tracked.elastic_modes[9].frequency if 9 in current_iteration_model_tracked.elastic_modes else np.nan
    f10elastic = current_iteration_model_tracked.elastic_modes[10].frequency if 10 in current_iteration_model_tracked.elastic_modes else np.nan

    msim7elastic = current_iteration_model_tracked.elastic_modes[7].extraction_order if 7 in current_iteration_model_tracked.elastic_modes else np.nan
    msim8elastic = current_iteration_model_tracked.elastic_modes[8].extraction_order if 8 in current_iteration_model_tracked.elastic_modes else np.nan
    msim9elastic = current_iteration_model_tracked.elastic_modes[9].extraction_order if 9 in current_iteration_model_tracked.elastic_modes else np.nan
    msim10elastic = current_iteration_model_tracked.elastic_modes[10].extraction_order if 10 in current_iteration_model_tracked.elastic_modes else np.nan

    mac7elastic, mac8elastic, mac9elastic, mac10elastic = np.diag(crossmac_tracked[6:10, 6:10])


    parametric_result_modal = ParameterResultModal(
            name=str(p_name),
            value=p_value,
            modal_results=current_iteration_model_tracked,
            frb_max=frb_max, 
            rb_mode=rb_mode,   
            f7elastic=f7elastic,
            f8elastic=f8elastic,
            f9elastic=f9elastic,
            f10elastic=f10elastic,            
            mac7elastic=mac7elastic,
            mac8elastic=mac8elastic,
            mac9elastic=mac9elastic,
            mac10elastic=mac10elastic,
            msim7elastic=msim7elastic,
            msim8elastic=msim8elastic,
            msim9elastic=msim9elastic,
            msim10elastic=msim10elastic,
        )

    return parametric_result_modal
    
def modal_postprocessing(parameter_name, parameter_values_sorted, opt_value, modal_results):

    increasing_prmtr = [v for v in parameter_values_sorted if v > opt_value]
    decreasing_prmtr = list(reversed([v for v in parameter_values_sorted if v <= opt_value]))

    previous_iteration = modal_results["Baseline"]
    for p_value in increasing_prmtr:

        current_iteration = modal_results[parameter_name][p_value]

        modal_results[parameter_name][p_value] = parametic_output_modal(parameter_name, p_value, previous_iteration.modal_results, current_iteration.modal_results)
        ### Aplly this pairing to flutter results?
        previous_iteration = current_iteration

    previous_iteration = modal_results["Baseline"]
    for p_value in decreasing_prmtr:

        current_iteration = modal_results[parameter_name][p_value]

        modal_results[parameter_name][p_value] = parametic_output_modal(parameter_name, p_value, previous_iteration.modal_results, current_iteration.modal_results)

        previous_iteration = current_iteration

    return modal_results

def flutter_postprocessing(parameter_name, parameter_values_sorted, opt_value, flutter_results):

    increasing_prmtr = [v for v in parameter_values_sorted if v > opt_value]
    decreasing_prmtr = list(reversed([v for v in parameter_values_sorted if v <= opt_value]))

    # for id, flutter_point in enumerate(flutter_results["Baseline"].flutter_mode):
    previous_iteration = flutter_results["Baseline"]
    for p_value in increasing_prmtr:

        current_iteration = flutter_results[parameter_name][p_value]

        flutter_results[parameter_name][p_value] = parametric_output_flutter(parameter_name, p_value, previous_iteration, current_iteration, flutter_point=id)

        previous_iteration = current_iteration
    
    previous_iteration = flutter_results["Baseline"]
    for p_value in decreasing_prmtr:

        current_iteration = flutter_results[parameter_name][p_value]

        flutter_results[parameter_name][p_value] = parametric_output_flutter(parameter_name, p_value, previous_iteration, current_iteration, flutter_point=id)

        previous_iteration = current_iteration

    return flutter_results

def parametric_output_flutter(p_name, p_value, previous_iteration_model, current_iteration_model, flutter_point):


    current_iteration_model_tracked = mode_trackingX(previous_iteration_model, current_iteration_model, flutter_point, tracking = True)


    return current_iteration_model_tracked ##############################################################################################################################################################################################################################################################################################
    class ParameterResultFlutter:
        name: str = None
        value: float = np.nan
        flutter_results: object = None
        flutter_sorted_points: object = None
        critical_flutter_results: object = None
        flutter_results_eigv: object = None
        flutter_sorted_points_eigv: object = None
        flutter_mode: list[float] = field(default_factory=list)
        flutter_velocity: list[float] = field(default_factory=list)
        flutter_frequency: list[float] = field(default_factory=list)
        flutter_coupling_mode: float = np.nan
        flutter_coupling_macx_max: float = np.nan
        # flutter_rb_mode : float = np.nan
        # flutter_rb_macx_max: float = np.nan
        # flutter_rb_coalescing : float = np.nan

def parameter_result_flutter(p_name, p_value, results_baseline, mode_shapes_request=None, prmtr_filexist=False):

    pid = os.getpid()  # Get the current process ID

    aero_bdf = nastran_model.from_file("aero/flutter-h=500m+.bdf")
    aero_bdf.to_file(f"aero/flutter-h=500m-eigv-{pid}.bdf")

    with open('flutter-deck-hybrid-p.bdf', 'r') as f:
        deck = f.read()
        deck = deck.replace(f'{{p}}', str(f"{p_name}/{p_name}{p_value:.3f}"))
        deck = deck.replace(f'{{pid}}', str(f"{pid}"))

        with open(f'flutter-deck-hybrid-p-{pid}.bdf', 'w') as f_out:
            f_out.write(deck)

    os.makedirs(f'out/parametric-study/{p_name}/{p_name}{p_value:.3f}', exist_ok=True)
    parametric_flutter_out_filename = f'flutter-deck-hybrid-p-{pid}.bdf'
    parametric_flutter_out_filepath = f'out/parametric-study/{p_name}/{p_name}{p_value:.3f}/flutter-deck-hybrid-p-{p_name}{p_value:.3f}'

    ################################################################################################################################################################################################################################################################################################################################################
    parametric_flutter_eigv_out_filepath = f'out/parametric-study/{p_name}/{p_name}{p_value:.3f}/flutter-deck-hybrid-p-{p_name}{p_value:.3f}-eigv'  ################ -eigv !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ################################################################################################################################################################################################################################################################################################################################################

    parametric_flutter_results, parametric_flutter_sorted_points = run_utils.filexist_runastran_read(parametric_flutter_out_filename, parametric_flutter_eigv_out_filepath, flutter = True,  filexist=prmtr_filexist, memory=0.5)    
    parametric_critical_flutter_points = nastran_results.nastran_results(flutter_points=nastran_results.find_critical_flutter_points(parametric_flutter_results.flutter_points))

    critical_flutter_points, divergence_points = nastran_results.find_critical_flutter_points(parametric_flutter_results.flutter_points)
    parametric_critical_flutter_points = nastran_results.nastran_results(flutter_points=critical_flutter_points)
    parametric_divergence_points = nastran_results.nastran_results(flutter_points=divergence_points)

    divergence_mode = [next(iter(flutter_point.modes)) for flutter_point in parametric_divergence_points.flutter_points] if parametric_divergence_points.flutter_points else [np.nan]
    divergence_velocity = [flutter_point.velocity for flutter_point in parametric_divergence_points.flutter_points] if parametric_divergence_points.flutter_points else [np.nan]
    # divergence_frequency = [flutter_point.modes[next(iter(flutter_point.modes))].frequency for flutter_point in parametric_divergence_points.flutter_points]

    if verify_critical_flutter_points(parametric_critical_flutter_points):

        flutter_mode = [next(iter(flutter_point.modes)) for flutter_point in parametric_critical_flutter_points.flutter_points]
        critical_velocities = [flutter_point.velocity for flutter_point in parametric_critical_flutter_points.flutter_points]
        flutter_frequency = [flutter_point.modes[next(iter(flutter_point.modes))].frequency for flutter_point in parametric_critical_flutter_points.flutter_points]

        aero_bdf = nastran_model.from_file("aero/flutter-h=500m+.bdf")
        density_table = list(aero_bdf.flfact_cards[1].factors)
        mach_table = list(aero_bdf.flfact_cards[2].factors)
        velocity_table = list(aero_bdf.flfact_cards[3].factors)
        while density_table and density_table[-1] is None:
            density_table.pop()        
        while mach_table and mach_table[-1] is None:
            mach_table.pop()        
        while velocity_table and velocity_table[-1] is None:
            velocity_table.pop()
        
        mach_points = list(np.interp(critical_velocities, velocity_table, mach_table))
        density_points = list(np.interp(critical_velocities, velocity_table, density_table))
        velocity_points = [-abs(v) for v in critical_velocities]

        # aero_out_bdf = nastran_model()
        aero_bdf.flfact_cards[1] = nastran_cards.flfact_card(density_points)
        aero_bdf.flfact_cards[2] = nastran_cards.flfact_card(mach_points)
        aero_bdf.flfact_cards[3] = nastran_cards.flfact_card(velocity_points)

        aero_bdf.to_file(f"aero/flutter-h=500m-eigv-{pid}.bdf")

        parametric_flutter_results_eigv, parametric_flutter_sorted_points_eigv = run_utils.filexist_runastran_read(parametric_flutter_out_filename, parametric_flutter_out_filepath, flutter = True,  filexist=prmtr_filexist, memory=0.5)    
        
        flutter_coupling_macx_max = [np.nan]*len(critical_velocities)
        flutter_coupling_mode = [np.nan]*len(critical_velocities)
        for idx, vel in enumerate(critical_velocities):
            flutter_coupling_macx_max[idx], flutter_coupling_mode[idx] = max((nastran_utils.compute_mode_similarity(parametric_flutter_sorted_points_eigv[idx].modes[aeroelastic_mode].eigenvector.values(), parametric_flutter_sorted_points_eigv[idx].modes[flutter_mode[idx]].eigenvector.values(), nastran_utils.macx), aeroelastic_mode) for aeroelastic_mode in [m for m in range(1, len(parametric_flutter_sorted_points_eigv[idx].modes)+1) if m != flutter_mode[idx]])

    else:
        parametric_flutter_results_eigv = None
        parametric_flutter_sorted_points_eigv = None
        flutter_mode = [np.nan]
        critical_velocities = [np.nan]
        flutter_frequency = [np.nan]
        flutter_coupling_macx_max = [np.nan]
        flutter_coupling_mode = [np.nan]

    if  not np.isnan(divergence_mode):
        print(f"{p_name}{p_value}: {divergence_mode} ({divergence_velocity}m/s)")



    os.remove(f'flutter-deck-hybrid-p-{pid}.bdf')
    os.remove(f'aero/flutter-h=500m-eigv-{pid}.bdf')
        
    prmtr_result_flutter =ParameterResultFlutter(
            name=str(p_name),
            value=p_value,
            flutter_results=parametric_flutter_results,
            flutter_sorted_points=parametric_flutter_sorted_points,
            critical_flutter_results=parametric_critical_flutter_points,
            flutter_results_eigv=parametric_flutter_results_eigv,
            flutter_sorted_points_eigv=parametric_flutter_sorted_points_eigv,
            divergence_mode=divergence_mode,
            divergence_velocity=divergence_velocity,
            flutter_mode=flutter_mode, ############################################################################################################################################################################################################################################################################################################
            flutter_velocity=critical_velocities,     ################################################################################################################################################################################################################################################################################################
            flutter_frequency=flutter_frequency,    ###############################################################################################################################################################################################################################################################################################
            flutter_coupling_macx_max=flutter_coupling_macx_max,
            flutter_coupling_mode=flutter_coupling_mode
            )
    
    return prmtr_result_flutter



def parametric_study_results(title, parameter_name, parameter_values, opt_value, static_results, modal_results, flutter_results, prmtr_filexist):

    if parameter_name not in static_results:
        static_results[parameter_name] = {}    
    if parameter_name not in modal_results:
        modal_results[parameter_name] = {}    
    if parameter_name not in flutter_results:
        flutter_results[parameter_name] = {}

    parameter_values_sorted = sorted(parameter_values)


    if not prmtr_filexist:
        with ProcessPoolExecutor(max_workers=7) as exe:
            futures = {exe.submit(parameter_result_static, parameter_name, p_value, prmtr_filexist=prmtr_filexist): p_value for p_value in parameter_values_sorted}
            for fut in as_completed(futures):
                p_value = futures[fut]
                prmtr_result_static = fut.result()
                static_results[parameter_name][p_value] = prmtr_result_static
    else:
        for p_value in parameter_values_sorted:
            static_results[parameter_name][p_value] = parameter_result_static(parameter_name, p_value, prmtr_filexist=prmtr_filexist)


    if not prmtr_filexist:
        with ProcessPoolExecutor(max_workers=7) as exe:
            futures = {exe.submit(parameter_result_modal, parameter_name, p_value, prmtr_filexist=prmtr_filexist): p_value for p_value in parameter_values_sorted}
            for fut in as_completed(futures):
                p_value = futures[fut]
                prmtr_result_modal = fut.result()
                modal_results[parameter_name][p_value] = prmtr_result_modal
    else:
        for p_value in parameter_values_sorted:
            modal_results[parameter_name][p_value] = parameter_result_modal(parameter_name, p_value, prmtr_filexist=prmtr_filexist)



    if not prmtr_filexist:
        with ProcessPoolExecutor(max_workers=7) as exe:
            futures = {exe.submit(parameter_result_flutter, parameter_name, p_value, flutter_results["Baseline"], mode_shapes_request=None, prmtr_filexist=prmtr_filexist): p_value for p_value in parameter_values_sorted}
            for fut in as_completed(futures):
                p_value = futures[fut]
                prmtr_result_flutter = fut.result()
                flutter_results[parameter_name][p_value] = prmtr_result_flutter
    else:
        for p_value in parameter_values_sorted:
            flutter_results[parameter_name][p_value] = parameter_result_flutter(parameter_name, p_value, flutter_results["Baseline"], mode_shapes_request=None, prmtr_filexist=prmtr_filexist)


    modal_results = modal_postprocessing(parameter_name, parameter_values_sorted, opt_value, modal_results)
    
    flutter_results = flutter_postprocessing(parameter_name, parameter_values_sorted, opt_value, flutter_results)


    return static_results, modal_results, flutter_results





def main():

    YES = True
    NO = False

    # k_FL = 150    #Fuselage Length Parameterization (k = 1 original BWB, k = 0 flyingwing)
    # k_SA = 0    #Wing Sweep Angle Parameterization (k = 0 original sweptback, k = 1 NoSweep, k = 2 sweptforward)
    # k_WS = 180    #Wing Aspect Ratio (Span) Parameterization (k = 0% NoSpan(debug), k = 100% original span, k = 200% double span)
    k_WC = 130    #Wing Aspect Ratio (Chord) Parameterization (k = -1 NoChord(debug), k = 0 original chord, k = 1 double chord)
    k_TR = 100    #Wing Taper Ratio (Chord) Parameterization (k = 0 original taper, k = 1 taper 0)
    # k_DI = 200    #Wing Dihedral Parameterization (k = 1 NoDihedral, k = 0 original dihedral, negative -> positive dihedral, positive -> negative dihedral)
    k_AC = 0    #Wing Aerodynamic Center Parameterization (k = -1 AC2LE, k = 0 original AC, k = 1 AC2TE)

    SA = [30, 35, 40, 45]  # Sweep angle parameterization values (100 - original sweep angle - 32.84 deg))
    DI = [1, 5, 9, 13] # Dihedral angle parameterization values (100 - original dihedral angle - 3.36 deg)
    AR = [5, 10, 15, 20]  # Aspect ratio parameterization values (100 - original span - 8.55)
    # K_FL = [0, 50, 100, 150]  # Fuselage length parameterization values (100% - original fuselage length)
    K_SGEO = [-1, -0.5, 0.5, 1]
    K_CGW = [-1, -0.5, 0.5, 1]
    K_CGF = [-1, -0.5, 0.5, 1]
    K_E = [1e-2, 1e-1, 1e1, 1e2]
    K_G = [1e-2, 1e-1, 1e1, 1e2]
    # K_IF = [0.5, 0.75, 1.25, 1.5]

    # aero_colors = plt.get_cmap('YlGn')(np.linspace(0.35, 0.7, 3)) # 3 Aerodynamic Shape Parameters
    # aero_colors = plt.get_cmap('cividis')(0.1, 0.2, 0.9) # 3 Aerodynamic Shape Parameters
    aero_colors = plt.get_cmap('tab20').colors[:6] # 3 Aerodynamic Shape Parameters
    # struct_colors = plt.get_cmap('Reds')(np.linspace(0.45, 0.8, 3)) # 3 Structural Parameters
    # struct_colors = plt.get_cmap('cividis')([0.0, 0.5, 1.0]) # 3 Structural Parameters
    struct_colors = plt.get_cmap('tab20').colors[6:12] # 3 Aerodynamic Shape Parameters

    cg_colors = plt.get_cmap('Purples')(np.linspace(0.45, 0.7, 2)) # 2 Mass (CG) Parameters

    ################## name: str  title: str                bvalue: float  configurations: list    color: tuple    relevant_configurations_vgf: list
    # parameter_specs = [
        #   ParameterSpec("SA",     "Sweep angle",                  32.84,    [30, 35, 40, 45],       aero_colors[0],   [0]),    #done
        #   ParameterSpec("DI",     "Dihedral angle",               3.36,     [1, 5, 9, 13],          aero_colors[1],   []),    #done
        #   ParameterSpec("AR",     "Aspect ratio",                 8.55,     [5, 10, 15, 20],        aero_colors[2],   [4]),    #upd
        #   ParameterSpec("K_SGEO", "Elastic axis",       0,        [-1, -0.5, 0.5, 1],     struct_colors[0], []),    #done
        #   ParameterSpec("K_E",    "Bending stiffness",   1,        [1e-2, 1e-1, 1e1, 1e2], struct_colors[1], [1]),    #upd
        #   ParameterSpec("K_G",    "Torsional stiffness", 1,        [1e-2, 1e-1, 1e1, 1e2], struct_colors[2], [1, 2]), #upd
        #   ParameterSpec("K_CGW",  "Move Wing CG",                 0,        [-1, -0.5, 0.5, 1],     cg_colors[0],     []),    #done
        #   ParameterSpec("K_CGF",  "Move Fuselage CG",             0,        [-1, -0.5, 0.5, 1],     cg_colors[1],     []),     #fix
    # ]    
    
    parameter_specs = [
          ParameterSpec("SA",     "Sweep angle",                  32.84,    [30, 35, 40, 45],       [aero_colors[0], aero_colors[1]] ,   [], []),    #done
          ParameterSpec("DI",     "Dihedral angle",               3.36,     [1, 5, 9, 13],          [aero_colors[2], aero_colors[3]],   [], []),    #done
          ParameterSpec("AR",     "Aspect ratio",                 8.55,     [5, 10, 15, 20],        [aero_colors[4], aero_colors[5]],   [], []),    #upd
          ParameterSpec("K_E",    "Bending stiffness",   1,        [0.1, 3.5, 6.5, 10], [struct_colors[0], struct_colors[1]],    [], []),    #upd
          ParameterSpec("K_G",    "Torsional stiffness", 1,        [0.1, 3.5, 6.5, 10], [struct_colors[2], struct_colors[3]],    [], []), #upd
          ParameterSpec("K_SGEO", "Elastic axis",       0,        [-1, -0.5, 0.5, 1],     [struct_colors[4], struct_colors[5]], [], []),    #done
    ]

    def intervals(start, step, n_samples):
        x = start + step * np.arange(1, n_samples + 1)
        if step < 0:
            x = x[::-1]

        return np.round(x, 3).tolist()
    
    # SA = [*intervals(33, -1., 12), 32.84, *intervals(33, 1., 12)]
    # DI = [*intervals(3.35, -0.25, 12), 3.36, *intervals(3.35, 0.25, 12)]
    # AR = [*intervals(8.6, -0.2, 12), 8.55, *intervals(8.6, 0.2, 12)]
    # K_E = [*intervals(1.0, -0.075, 12), 1.0, *intervals(1.0, 0.075, 12)]
    # K_G = [*intervals(1.0, -0.075, 12), 1.0, *intervals(1.0, 0.075, 12)]
    # K_SGEO = np.round(np.linspace(-1, 1, 25), 3).tolist()
    
    SA = [*intervals(33, -0.24, 50), 32.84, *intervals(33, 0.24, 50)]
    DI = [*intervals(3.35, -0.06, 50), 3.36, *intervals(3.35, 0.06, 50)]
    AR = [*intervals(8.55, -0.04, 50), 8.55, *intervals(8.55, 0.04, 50)]
    K_E = [*intervals(1.0, -0.018, 50), 1.0, *intervals(1.0, 0.018, 50)]
    # K_G = [*intervals(1.0, -0.012, 50), 1.0, *intervals(1.0, 0.012, 50)]
    K_G = K_E
    K_SGEO = np.round(np.linspace(-1, 1, 101), 3).tolist()
   
    parameter_specs = [
          ParameterSpec("SA",     r"$\Lambda_{c/4}$",                  32.84,    SA,       [aero_colors[0], aero_colors[1]],   [], []),    #done
          ParameterSpec("DI",     r"$\phi$",               3.36,     DI,          [aero_colors[2], aero_colors[3]],   [], []),    #done
          ParameterSpec("AR",     r"A\!R",                 8.55,     AR,        [aero_colors[4], aero_colors[5]],   [], []),    #upd
          ParameterSpec("K_E",    "EI",   1,        K_E, [struct_colors[0], struct_colors[1]],    [], []),    #upd
          ParameterSpec("K_G",    "GJ", 1,        K_G, [struct_colors[2], struct_colors[3]],    [14], [14]), #upd
          ParameterSpec("K_SGEO", r"$x_{\text{EA}}$",       0,        K_SGEO,     [struct_colors[4], struct_colors[5]], [], []),    #done
    ]
    
    # parameter_specs = [
    #       ParameterSpec("SA",     "Sweep angle",                  32.84,    [28.64, 34.09, 39.55, 45],       aero_colors[0],   []),    #done
    #       ParameterSpec("DI",     "Dihedral angle",               3.36,     [1, 4.82, 8.64, 13.73],          aero_colors[1],   []),    #done
    #       ParameterSpec("AR",     "Aspect ratio",                 8.55,     [5, 9.55, 16.36, 20.91],        aero_colors[2],   []),    #upd
    #       ParameterSpec("K_E",    "Bending stiffness",   1,        [0.1, 2.5, 7.0, 10], struct_colors[1],    []),    #upd
    #       ParameterSpec("K_G",    "Torsional stiffness", 1,        [0.1, 2.5, 7.0, 10], struct_colors[2],    []), #upd
    #       ParameterSpec("K_SGEO", "Elastic axis",       0,        [-1, -0.45, 0.45, 1],     struct_colors[0], []),    #done
    # ]

    # k=np.r_[np.linspace(0.1, 1, 6), np.linspace(1, 10, 7)[1:]]

    # parameter_specs = [
    #       ParameterSpec("SA",     "Sweep angle",                  32.84,    list(np.linspace(15, 45, 12)),       aero_colors[0],   [], []),    #done
    #       ParameterSpec("DI",     "Dihedral angle",               3.36,     list(np.linspace(1, 15, 12)),          aero_colors[1],   [], []),    #done
    #       ParameterSpec("AR",     "Aspect ratio",                 8.55,     list(np.linspace(5, 30, 12)),        aero_colors[2],   [], []),    #upd
    #       ParameterSpec("K_E",    "Bending stiffness",   1,        list(k), struct_colors[0],    [1], []),    #upd
    #       ParameterSpec("K_G",    "Torsional stiffness", 1,        list(k), struct_colors[2],    [1], []), #upd
    #       ParameterSpec("K_SGEO", "Elastic axis",       0,        list(np.linspace(-1, 1, 12)),     struct_colors[1], [], []),    #done
    # ]

    # AR = [*intervals(7, -0.05, 20)]
    # K_SGEO = [*intervals(0.9, 0.005, 20)]
    # AR = np.round(np.linspace(6.5, 6.7, 20), 3).tolist()
    # K_SGEO = np.round(np.linspace(0.92, 0.94, 20), 3).tolist()
    # K_E = np.round(np.linspace(0.3, 0.5, 20), 3).tolist()
    # K_G = np.round(np.linspace(0.1, 0.4, 20), 3).tolist()
    
    # K_G = [*intervals(0.55, -0.018, 25)]

    # parameter_specs = [
    #       ParameterSpec("K_G",    "GJ",   1,        K_G, [struct_colors[0], struct_colors[1]],    [9], [9]),    #upd
    # ]   

    parameter_specs_dict = {spec.name: spec for spec in parameter_specs}  # O(1) lookup by name


    SA = parameter_specs_dict.get("SA").configurations if parameter_specs_dict.get("SA") else [] # Sweep angle parameterization values (100 - original sweep angle - 32.84 deg))
    DI = parameter_specs_dict.get("DI").configurations if parameter_specs_dict.get("DI") else [] # Dihedral angle parameterization values (100 - original dihedral angle - 3.36 deg)
    AR = parameter_specs_dict.get("AR").configurations if parameter_specs_dict.get("AR") else []  # Aspect ratio parameterization values (100 - original span - 8.55)
    # K_FL = [0, 50, 100, 150]  # Fuselage length parameterization values (100% - original fuselage length)
    K_SGEO = parameter_specs_dict.get("K_SGEO").configurations if parameter_specs_dict.get("K_SGEO") else []
    K_CGW = parameter_specs_dict.get("K_CGW").configurations if parameter_specs_dict.get("K_CGW") else []
    K_CGF = parameter_specs_dict.get("K_CGF").configurations if parameter_specs_dict.get("K_CGF") else []
    K_E = parameter_specs_dict.get("K_E").configurations if parameter_specs_dict.get("K_E") else []
    K_G = parameter_specs_dict.get("K_G").configurations if parameter_specs_dict.get("K_G") else []
    # K_IF = [0.5, 0.75, 1.25, 1.5]

    run = NO
    if run:
        for sa in SA:
            stick_bdf_parameterization(f"prmtr/SA/SA{sa:.3f}", sa = sa)
        sa=100
        for di in DI:
            stick_bdf_parameterization(f"prmtr/DI/DI{di:.3f}", di = di)
        di = 100
        for ar in AR:
            stick_bdf_parameterization(f"prmtr/AR/AR{ar:.3f}", ar = ar)
        ar = 100
        # for k_fl in K_FL:
        #     parameterization("KFL", k_fl, aero_parameterization.fuselage_length)
        for k_sgeo in K_SGEO:
            stick_bdf_parameterization(f"prmtr/K_SGEO/K_SGEO{k_sgeo:.3f}", k_sgeo = k_sgeo)
        k_sgeo = 0
        for k_E in K_E:
            stick_bdf_parameterization(f"prmtr/K_E/K_E{k_E:.3f}", k_E = k_E)
        k_E = 1
        for k_G in K_G:
            stick_bdf_parameterization(f"prmtr/K_G/K_G{k_G:.3f}", k_G = k_G)
        k_G = 1
        for k_cgw in K_CGW:
            stick_bdf_parameterization(f"prmtr/K_CGW/K_CGW{k_cgw:.3f}", k_cgw = k_cgw)
        k_cgw = 0
        for k_cgf in K_CGF:
            stick_bdf_parameterization(f"prmtr/K_CGF/K_CGF{k_cgf:.3f}", k_cgf = k_cgf)  
        k_cgf = 0

    static_results = None
    modal_results = None
    flutter_results = None

    static_results = load_results(filename="cache/parametric_results_static.pkl.gz")    
    modal_results = load_results(filename="cache/parametric_results_modal.pkl.gz")    
    flutter_results = load_results(filename="cache/parametric_results_flutter.pkl.gz")
    if flutter_results is None:
        baseline_static_result, baseline_modal_result, baseline_flutter_result  = hybrid_flutter_results('flutter-deck-hybrid-opt.bdf', 'out/parametric-study/flutter-deck-hybrid-opt', filexist = YES)
        static_results = {"Baseline": baseline_static_result}
        modal_results = {"Baseline": baseline_modal_result}
        flutter_results = {"Baseline": baseline_flutter_result}

        prmtr_filexist = YES
        static_results, modal_results, flutter_results = (parametric_study_results("Sweep angle", "SA", SA, 32.84, static_results, modal_results, flutter_results, prmtr_filexist=prmtr_filexist) if SA else (static_results, modal_results, flutter_results))
        static_results, modal_results, flutter_results = (parametric_study_results("Dihedral angle", "DI", DI, 3.36, static_results, modal_results, flutter_results, prmtr_filexist=prmtr_filexist) if DI else (static_results, modal_results, flutter_results))
        static_results, modal_results, flutter_results = (parametric_study_results("Aspect ratio", "AR", AR, 8.55, static_results, modal_results, flutter_results, prmtr_filexist=prmtr_filexist) if AR else (static_results, modal_results, flutter_results))
        static_results, modal_results, flutter_results = (parametric_study_results("Bending stiffness", "K_E", K_E, 1, static_results, modal_results, flutter_results, prmtr_filexist=prmtr_filexist) if K_E else (static_results, modal_results, flutter_results))
        static_results, modal_results, flutter_results = (parametric_study_results("Torsional stiffness", "K_G", K_G, 1, static_results, modal_results, flutter_results, prmtr_filexist=prmtr_filexist) if K_G else (static_results, modal_results, flutter_results))
        static_results, modal_results, flutter_results = (parametric_study_results("Elastic axis", "K_SGEO", K_SGEO, 0, static_results, modal_results, flutter_results, prmtr_filexist=prmtr_filexist) if K_SGEO else (static_results, modal_results, flutter_results))

        # static_results, modal_results, flutter_results = (parametric_study_results("Move Wing CG", "K_CGW", K_CGW, 0, static_results, modal_results, flutter_results, prmtr_filexist=prmtr_filexist) if K_CGW else (static_results, modal_results, flutter_results))
        # static_results, modal_results, flutter_results = (parametric_study_results("Move Fuselage CG", "K_CGF", K_CGF, 0, static_results, modal_results, flutter_results, prmtr_filexist=prmtr_filexist) if K_CGF else (static_results, modal_results, flutter_results))
        
        save_results(static_results, filename="cache/parametric_results_static.pkl.gz")
        save_results(modal_results, filename="cache/parametric_results_modal.pkl.gz")
        save_results(flutter_results, filename="cache/parametric_results_flutter.pkl.gz")



    plot_parameter = ["SA", "DI", "AR", "K_E", "K_G", "K_SGEO"]
    # plot_parameter = ["K_E"]

    plt.ion()

    no_configurations = (max(len(SA), len(DI), len(AR), len(K_SGEO), len(K_CGW), len(K_CGF), len(K_E), len(K_G)))

    plot_modal_graphs=YES
    if plot_modal_graphs:
        fig_modal_freq_evolution, ax_modal_freq_evolution= plt.subplot_mosaic([['f7elastic', 'f8elastic', 'f9elastic', 'f10elastic']], figsize=(15.0, 6.0), layout='constrained', gridspec_kw={'hspace': 0.1})
        # fig_modal_freq_evolution, ax_modal_freq_evolution= plt.subplot_mosaic([['f7elastic', 'f8elastic', 'f9elastic', 'f10elastic'], ['mac7elastic', 'mac8elastic', 'mac9elastic', 'mac10elastic']], figsize=(15.0, 6.0), layout='constrained', gridspec_kw={'hspace': 0.1})
        attr_modal = ['frb_max', 'f7elastic']
        modal_evolution_freq_legend_lines = plot_flutter.evolution_plot(ax_modal_freq_evolution['f7elastic'], modal_results, parameter_specs_dict, ['f7elastic'], plot_parameter, flutter_point=0)
        modal_evolution_freq_legend_lines = plot_flutter.evolution_plot(ax_modal_freq_evolution['f8elastic'], modal_results, parameter_specs_dict, ['f8elastic'], plot_parameter, flutter_point=0)
        modal_evolution_freq_legend_lines = plot_flutter.evolution_plot(ax_modal_freq_evolution['f9elastic'], modal_results, parameter_specs_dict, ['f9elastic'], plot_parameter, flutter_point=0)
        modal_evolution_freq_legend_lines = plot_flutter.evolution_plot(ax_modal_freq_evolution['f10elastic'], modal_results, parameter_specs_dict, ['f10elastic'], plot_parameter, flutter_point=0)
        plot_flutter.modal_evolution_freq_axes(ax_modal_freq_evolution, fig_modal_freq_evolution, modal_results, no_configurations, parameter_specs_dict, modal_evolution_freq_legend_lines, plot_parameter, flutter_point=0)
        # plot_flutter.save_plots(fig_modal_freq_evolution, fig_name="modal_evolution_freq.pdf", output_dir="plots/parametric-study/evolution")

        # fig_modal_mac_evolution, ax_modal_freq_evolution= plt.subplot_mosaic([['mac7elastic', 'mac8elastic', 'mac9elastic', 'mac10elastic']], figsize=(15.0, 6.0), layout='constrained', gridspec_kw={'hspace': 0.1})
        # modal_evolution_mac_legend_lines = plot_flutter.evolution_plot(ax_modal_freq_evolution['mac7elastic'], modal_results, parameter_specs_dict, ['mac7elastic'], plot_parameter, flutter_point=0)
        # modal_evolution_mac_legend_lines = plot_flutter.evolution_plot(ax_modal_freq_evolution['mac8elastic'], modal_results, parameter_specs_dict, ['mac8elastic'], plot_parameter, flutter_point=0)
        # modal_evolution_mac_legend_lines = plot_flutter.evolution_plot(ax_modal_freq_evolution['mac9elastic'], modal_results, parameter_specs_dict, ['mac9elastic'], plot_parameter, flutter_point=0)
        # modal_evolution_mac_legend_lines = plot_flutter.evolution_plot(ax_modal_freq_evolution['mac10elastic'], modal_results, parameter_specs_dict, ['mac10elastic'], plot_parameter, flutter_point=0)
        # plot_flutter.modal_evolution_mac_axes(ax_modal_freq_evolution, fig_modal_freq_evolution, modal_results, no_configurations, parameter_specs_dict, modal_evolution_mac_legend_lines, plot_parameter, flutter_point=0)
        # plot_flutter.save_plots(fig_modal_mac_evolution, fig_name="modal_evolution_mac.pdf", output_dir="plots/parametric-study/evolution")
        plot_flutter.save_plots(fig_modal_freq_evolution, fig_name="modal_evolution_freq.pdf", output_dir="plots/parametric-study/evolution")


    plot_flutter_graphs=YES
    fig_flutter_evolution, axes_flutter_evolution = plt.subplot_mosaic([['flutter_velocity', 'flutter_frequency', 'legend']], figsize=(15, 6), layout='constrained', gridspec_kw={"width_ratios": [1, 1, 0.3]})
    # fig_flutter_evolution, axes_flutter_evolution = plt.subplot_mosaic([['flutter_velocity', 'flutter_frequency', 'flutter_coupling_macx_max']], figsize=(15, 6), layout='constrained')
    # fig_flutter_evolution, axes_flutter_evolution = plt.subplot_mosaic([['flutter_velocity', 'flutter_frequency', 'flutter_coupling_macx_max', 'legend']], figsize=(15, 6), layout='constrained', gridspec_kw={"width_ratios": [1, 1, 1, 0.4]})

    flutter_evolution_legend_lines = {}
    for idx, flutter_point in enumerate(flutter_results["Baseline"].flutter_mode):
        # if idx == 1:
        #     continue
        if plot_flutter_graphs:
            flutter_evolution_legend_lines[flutter_point] = plot_flutter.evolution_plot(axes_flutter_evolution["flutter_velocity"], flutter_results, parameter_specs_dict, ['flutter_velocity'], plot_parameter, flutter_point=idx)
            flutter_evolution_legend_lines[flutter_point] = plot_flutter.evolution_plot(axes_flutter_evolution["flutter_frequency"], flutter_results, parameter_specs_dict, ['flutter_frequency'], plot_parameter, flutter_point=idx)
            # flutter_evolution_legend_lines = plot_flutter.evolution_plot(axes_flutter_evolution["flutter_coupling_macx_max"], flutter_results, parameter_specs_dict, ['flutter_coupling_macx_max'], plot_parameter, flutter_point=idx)
    merged_legend_lines = [
    handle
    for legend_lines in flutter_evolution_legend_lines.values()
    for handle in legend_lines
    ]

    plot_flutter.flutter_evolution_axes(axes_flutter_evolution, fig_flutter_evolution, flutter_results, no_configurations, parameter_specs_dict, merged_legend_lines, plot_parameter, flutter_point=0)
    plot_flutter.save_plots(fig_flutter_evolution, fig_name=f"flutter_evolution.pdf", output_dir="plots/parametric-study/evolution")

    # plot_divergence_graphs=YES
    # for idx, flutter_point in enumerate(flutter_results["Baseline"].flutter_mode):
    #     if idx == 1:
    #         continue
    #     if plot_divergence_graphs:
    #         fig_flutter_evolution, axes_flutter_evolution = plt.subplot_mosaic([['divergence_velocity']], figsize=(15, 6), layout='constrained')
    #         flutter_evolution_legend_lines = plot_flutter.evolution_plot(axes_flutter_evolution["divergence_velocity"], flutter_results, parameter_specs_dict, ['divergence_velocity'], plot_parameter, flutter_point=idx)
    #         # plot_flutter.flutter_evolution_axes(axes_flutter_evolution, fig_flutter_evolution, flutter_results, no_configurations, parameter_specs_dict, flutter_evolution_legend_lines, plot_parameter, flutter_point=idx)
    #         plot_flutter.save_plots(fig_flutter_evolution, fig_name=f"divergence_evolution_{idx}.pdf", output_dir="plots/parametric-study/evolution")

    

    for parameter in parameter_specs:
        for configuration in parameter.relevant_configurations_vgf:
            # VGF and Locus Plots
            fig_rlvgf, ax_rlvgf = plt.subplot_mosaic([["locus",  "damp"], ["locus",  "freq"]], figsize=(10, 6), layout="constrained")
            rlvgf_legend_lines = []
            rlvgf_legend_lines = plot_flutter.rlvgf_plot(ax_rlvgf, flutter_results, parameter, rlvgf_legend_lines)
            plot_flutter.rlvgf_axes(fig_rlvgf, ax_rlvgf, rlvgf_legend_lines)
            plot_flutter.save_plots(fig_rlvgf, fig_name="rlvgf_" + parameter.name + ".pdf", output_dir="plots/parametric-study")
            # plt.close(fig_rlvgf)

    # for parameter in parameter_specs:
    #     if parameter.relevant_configurations_static:
    #         # Static load case
    #         deformations_fig_hybrid, deformations_axes_hybrid = plt.subplot_mosaic([['disp'], ['rot']], figsize=(12, 8), layout='constrained')
    #         static_legend_lines = []
    #         static_legend_lines = plot_flutter.static_plot(deformations_axes_hybrid, results, parameter, static_legend_lines)
    #         plot_flutter.static_axes(deformations_fig_hybrid, deformations_axes_hybrid, static_legend_lines)
    #         # plot_flutter.save_plots(deformations_fig_hybrid, fig_name="static_case_" + parameter.name + ".pdf", output_dir="plots/parametric-study")
    #         # plt.close(deformations_fig_hybrid)



    print("Finished all parametric studies.")
    plt.show(block = True)



if __name__ == '__main__':
    main()