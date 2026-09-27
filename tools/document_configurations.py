"""Add complete human-readable metadata to every SHS JSON configuration."""
import json
from pathlib import Path

ROOTS=(Path("configs"),Path("tools/config"),Path("src/shs/config"))
UNITS={"_m":"metres","_s":"seconds","_K":"kelvin","_A":"amperes","_V":"volts","_W":"watts","_T":"tesla","_F":"farads","_Hz":"hertz","_ohm":"ohms","_ohm_m":"ohm-metres","_A_per_m2":"amperes per square metre"}
SPECIAL={
"name":"Human-readable identifier.","geometry":"Geometry JSON filename loaded from configs/geometry.","material":"Material JSON filename or material-property object.","simulation_config":"Base production simulation JSON built before this runner applies documented experiment controls.",
"temperature":"Initial uniform temperature in kelvin.","current":"Initial or commanded transport current.","duration":"Requested physical duration.","dt":"Physical integration timestep.","steps":"Number of physical steps executed by this runner.","sample_every":"Save diagnostics every this many accepted steps.",
"duration_multiplier":"Runs for this multiple of the referenced simulation duration; the simulation remains the sole owner of the physical duration and timestep.",
"directory":"Base output directory; runners reserve a numbered sibling if it already exists.","output_directory":"Base output directory; runners reserve a numbered sibling if it already exists.","enabled":"Enables this model or output feature.",
"nx":"Number of mesh nodes along x.","ny":"Number of mesh nodes along y.","width":"Physical width in metres.","height":"Physical height in metres.","thickness":"Film thickness in metres.","type":"Selects the geometry, boundary, contact, or model type.",
"critical_current_A":"Zero-field SIS critical current; null requests Ambegaokar-Baratoff calculation.","normal_resistance_ohm":"Normal/shunt resistance of the SIS junction.","capacitance_F":"Junction capacitance.","frequency_Hz":"Microwave drive frequency.","microwave_amplitudes_A":"Delivered AC-current amplitudes included in the sweep.","sweep_dc_values":"DC drive values included in the sweep.",
"source_contact":"Named contact where positive current enters.","sink_contact":"Named contact where positive current leaves.","drive_mode":"Selects current or voltage boundary control.","include_self_field":"Feeds the induced thin-film vector potential back into the coupled solve.","include_scalar_potential":"Includes electrostatic phase evolution in TDGL.",
"bath_temperature":"Thermal-reservoir temperature toward which the film relaxes.","thermal_relaxation_rate":"Linear heat-removal coefficient coupling the film to its bath.",
"normal_conductivity_model":"Selects how normal-carrier conductivity varies with superconducting state and temperature.","normalization":"Selects the dimensional scaling used by the TDGL equation.","temperature_model":"Selects how temperature modifies the TDGL coefficients.",
"u":"Dimensionless TDGL order-parameter relaxation constant.","gamma":"Generalized TDGL inelastic-relaxation parameter; zero selects the simpler gapless form.","kappa":"Ginzburg-Landau penetration-depth to coherence-length ratio.",
"gauge":"Vector-potential gauge used to represent the applied magnetic field.","source_stride":"Spatial subsampling for magnetic diagnostics; larger values are faster and less resolved.",
"amplitude_floor":"Minimum order-parameter amplitude at which phase winding is considered numerically reliable.","kymograph_centerline":"Chooses the horizontal or vertical line sampled through time.","vector_stride":"Grid subsampling for plotted current and field arrows.",
"fps":"Playback frames per second for generated GIFs.","dpi":"Raster resolution of generated figures.","charge":"Integer vortex winding: +1 is a vortex and -1 an antivortex.","core_radius_cells":"Initial vortex-core radius measured in mesh cells.",
"backend":"Selects the numerical implementation for this solver component.","omega":"Successive-over-relaxation weight for the iterative solve.","residual_check_interval":"Iterations between residual evaluations.","check_interval":"Iterations between convergence-controller decisions.",
"prediction_window":"Recent residual values used to estimate convergence behavior.","stall_ratio":"Residual-reduction ratio above which progress is stalled.","fast_ratio":"Residual-reduction ratio below which convergence is fast.","healthy_ratio":"Residual-reduction ratio below which convergence is healthy.","oscillation_window":"Residual-history length used to detect oscillation.",
"stability_safety_factor":"Safety margin used when choosing an explicit stable timestep.","screening_max_iterations":"Maximum self-consistent magnetic-screening iterations per coupled step.","screening_step_size":"Under-relaxation weight for magnetic-screening updates.","screening_source_stride":"Spatial subsampling for self-field feedback; one uses every source cell.","reasonable_iterations":"Iteration count treated as acceptable when scoring a benchmark.",
"region_type":"Behavior assigned to this geometric region.","contact_type":"Electrical boundary condition applied at this contact.","voltage_left":"Potential imposed at the left transport boundary.","voltage_right":"Potential imposed at the right transport boundary.","reference_voltage":"Reference that removes the arbitrary constant from scalar potential.",
"coherence_length":"Superconducting coherence length used for spatial TDGL scaling.","penetration_depth":"London penetration depth used for magnetic scaling and validity checks.","Tc":"Critical temperature above which equilibrium superconductivity vanishes.","gl_alpha":"Linear Ginzburg-Landau free-energy coefficient.","gl_beta":"Nonlinear Ginzburg-Landau free-energy coefficient.","tdgl_u":"Material default for the dimensionless TDGL relaxation constant.",
"thermal_conductivity":"Material heat-conduction coefficient.","heat_capacity":"Volumetric heat capacity used by thermal evolution.","normal_resistivity":"Normal-state electrical resistivity.","electrode":"Junction electrode in which the vortex is seeded.","microwave_steps_per_cycle":"Integration samples per microwave period; more samples improve resolution at greater cost.","initial_phase_rad":"Initial junction phase difference in radians.","rtol":"Relative numerical tolerance.","atol":"Absolute numerical tolerance in the associated quantity's units.",
"dynamic_plot_cycles":"Number of final microwave cycles expanded in transient and local-current diagnostic plots.",
"model":"Selects a single film temperature or separate electron and phonon temperatures.",
"electron_heat_capacity_fraction":"Fraction of the material heat capacity assigned to electrons; the remainder is assigned to phonons.",
"electron_thermal_conductivity_fraction":"Fraction of the material thermal conductivity assigned to electrons; the remainder is assigned to phonons.",
"electron_phonon_coupling_W_m3_K":"Linear electron-to-phonon energy-transfer coefficient in watts per cubic metre per kelvin.",
"phonon_escape_rate_per_s":"Rate at which film phonons relax toward the substrate bath; inverse seconds.",
"time_integrator":"Euler is first order; Heun uses a gauge-compatible predictor-corrector for improved small-step TDGL accuracy.",
}

def leaves(value,prefix=""):
    if isinstance(value,dict):
        for key,item in value.items():
            if key=="_documentation": continue
            yield from leaves(item,f"{prefix}.{key}" if prefix else key)
    elif isinstance(value,list):
        if value and isinstance(value[0],dict):
            for key in sorted(set().union(*(item.keys() for item in value if isinstance(item,dict)))):
                examples=[item[key] for item in value if isinstance(item,dict) and key in item]
                for example in examples:
                    yield from leaves(example,f"{prefix}[].{key}")
        else: yield prefix,value
    else: yield prefix,value

def explain(path,value):
    key=path.split(".")[-1].replace("[]","")
    if key in SPECIAL: return SPECIAL[key]
    unit=next((label for suffix,label in UNITS.items() if key.endswith(suffix)),None)
    words=key.replace("_"," ")
    if unit: return f"Sets {words}; units are {unit}."
    if key.startswith("maximum_") or key.startswith("max_"): return f"Upper limit for {words.removeprefix('maximum ').removeprefix('max ')}."
    if key.startswith("minimum_") or key.startswith("min_"): return f"Lower limit for {words.removeprefix('minimum ').removeprefix('min ')}."
    if "tolerance" in key: return f"Numerical or physical acceptance tolerance for {words.replace(' tolerance','')}."
    if isinstance(value,bool): return f"Enables or disables {words}."
    if "fraction" in key: return f"Dimensionless fraction controlling {words.replace(' fraction','')}."
    return f"Configures {words}."

def metadata(path,data):
    relative=path.as_posix()
    if path.parts[0]=="configs":
        kind=path.parts[1]
        purpose=f"Defines the reusable {kind[:-1]} named {data.get('name',path.stem)!r}."
        relies="Loaded by a simulation or experiment configuration; packaged defaults fill only omitted optional simulation controls."
        command="Not run directly; reference this filename from a simulation or tool configuration."
    elif path.parts[0]=="tools":
        purpose=f"Runs the {path.stem.replace('_',' ')} experiment or benchmark."
        relies=f"Uses {data.get('simulation_config','only the physics values in this file')}."
        runner="tdgl_diagnostics.py"
        if "sis_" in path.stem: runner="sis_vortex_microwave_diagnostics.py" if ("shapiro" in path.stem or "vortex" in path.stem) else "sis_junction_diagnostics.py"
        elif path.stem=="reduced_optical_tweezers": runner="reduced_vortex_diagnostics.py"
        elif path.stem=="coupled_convergence_benchmark": runner="benchmark_adaptive_coupled_updated.py"
        elif path.stem=="josephson_reversal_protocol": runner="josephson_reversal_protocol.py"
        elif path.stem=="josephson_current_sweep": runner="josephson_diagnostics.py"
        elif path.stem=="optical_vortex_two_temperature_search": runner="optical_vortex_search.py"
        command=f"python tools/{runner} --config {relative}"
    else:
        purpose="Defines packaged fallback values used only when a simulation omits an optional section or field."
        relies="Loaded by src/shs/utils/defaults.py; explicit values in a simulation configuration take precedence."
        command="Not run directly; edit only when changing package-wide fallback behavior."
    fields={p:explain(p,v) for p,v in leaves(data)}
    return {"purpose":purpose,"relies_on":relies,"run_command":command,
            "authority":"Values in this file are authoritative for its runner. A sweep value overrides the corresponding base-simulation value only for that named case; effective values are recorded in summary.json.","fields":fields}

def main():
    for root in ROOTS:
        for path in sorted(root.rglob("*.json")):
            data=json.loads(path.read_text(encoding="utf-8")); data["_documentation"]=metadata(path,data)
            ordered={"_documentation":data.pop("_documentation"),**data}
            path.write_text(json.dumps(ordered,indent=2)+"\n",encoding="utf-8")
            print(path)

if __name__=="__main__": main()
