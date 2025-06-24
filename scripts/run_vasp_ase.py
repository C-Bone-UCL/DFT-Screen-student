import os, sys, shutil, argparse, glob, subprocess, time
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun

class VaspError(Exception): pass

def check_finished(f="OUTCAR"):
    """
    Check if the VASP run has reached the required accuracy.
    """
    if not os.path.isfile(f):
        return False
    with open(f) as h:
        return any("reached required accuracy" in ln for ln in h)

def get_cycle():
    """ 
    Inspects the filesystem for previously completed cycle outputs
    By finding the highest-numbered cycle that has already run, resumes 
    from the next logical step, rather than starting over from cycle 0
    """
    outs = sorted(glob.glob("Cycle*ISIF3.OUTCAR"))
    return 0 if not outs else int(outs[-1].split('_')[0][5:]) + 1

def get_parallelization_settings():
    """
    Get proper parallelization settings based on available processors
    """
    nslots = int(os.environ.get("NSLOTS", "1"))
    
    # Conservative parallelization to avoid errors
    if nslots <= 1:
        return {"NPAR": 1, "KPAR": 1}
    elif nslots <= 4:
        return {"NPAR": 1, "KPAR": min(2, nslots)}
    elif nslots <= 8:
        return {"NPAR": 2, "KPAR": min(4, nslots//2)}
    else:
        return {"NPAR": 4, "KPAR": min(8, nslots//4)}

def validate_potcar_availability(structure, potcar_dir):
    """
    Check if POTCAR files exist for all elements in the structure
    """
    missing_elements = []
    for element in structure.symbol_set:
        potcar_path = os.path.join(potcar_dir, element, 'POTCAR')
        if not os.path.exists(potcar_path):
            missing_elements.append(element)
    
    if missing_elements:
        raise VaspError(f"Missing POTCAR files for elements: {missing_elements}")
    
    print(f"POTCAR validation passed for elements: {sorted(structure.symbol_set)}")

def run_step(structure, incar_settings, tag):
    """
    Run a single VASP step with the given structure and settings.
    
    1. Prepares the VASP input files (INCAR, KPOINTS, POSCAR, POTCAR).
    Using MPRelaxSet for standard settings.
    2. Executes the VASP command defined in the environment variable VASP_COMMAND.
    3. Checks for successful completion and copies output files.
    4. Returns the name of the CONTCAR file for the next step.
    """
    
    print(f"\n####Preparing VASP step: '{tag}'####\n")

    # Create dynamic POTCAR settings for all elements in the structure
    user_potcar_settings = {element: element for element in structure.symbol_set}

    calc_set = MPRelaxSet(
        structure,
        user_incar_settings=incar_settings,
        force_gamma=False,
        user_potcar_settings=user_potcar_settings  # Now handles all elements dynamically
    )
    # Default is the PBE fucntional
    
    # Make INCAR, KPOINTS, and POSCAR files normally
    calc_set.incar.write_file("INCAR")
    calc_set.kpoints.write_file("KPOINTS")
    structure.to(fmt="poscar", filename="POSCAR")
    
    potcar_dir = os.environ["PMG_VASP_PSP_DIR"]

    # This part collects all unique species in the structure
    # and writes a combined POTCAR file
    # Done manually because pymatgen's default POTCAR generation
    # does not support our potentials directory structure
    symbols_in_order = []
    for site in structure:
        symbol = site.specie.symbol
        if symbol not in symbols_in_order:
            symbols_in_order.append(symbol)
    
    # Check if POTCAR files exist before proceeding
    for symbol in symbols_in_order:
        potcar_path = os.path.join(potcar_dir, symbol, 'POTCAR')
        if not os.path.exists(potcar_path):
            raise VaspError(f"POTCAR file not found for element {symbol} at {potcar_path}")
            
    with open("POTCAR", 'wb') as potcar_file:
        for symbol in symbols_in_order:
            potcar_path = os.path.join(potcar_dir, symbol, 'POTCAR')
            with open(potcar_path, 'rb') as individual_potcar:
                shutil.copyfileobj(individual_potcar, potcar_file)

    vasp_command_str = os.environ["VASP_COMMAND"]
    command_list = vasp_command_str.split()

    print(f"####Executing command: {' '.join(command_list)}####")

    with open('vasp_out', 'w') as f_out:
        result = subprocess.run(command_list, stdout=f_out, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"####VASP CRASHED for tag {tag}####")
        print("####stderr from VASP:####")
        print(result.stderr)
        raise VaspError(f"VASP execution failed for tag {tag}. Check outputs.")

    # to record relaxation trajectory
    shutil.copy("OUTCAR",  f"{tag}.OUTCAR")
    shutil.copy("CONTCAR", f"{tag}.CONTCAR")
    shutil.copy("vasprun.xml", f"{tag}.vasprun.xml")
    
    return "CONTCAR"

def workflow(cif, potcar_dir):
    """
    Main workflow function to run VASP calculations on a given CIF file.

    1. Reads the CIF file to create a pymatgen Structure object.
    2. Sets up the VASP calculation parameters based on the environment.
    3. Runs the VASP calculations in a loop until convergence is achieved.
        a. Uses ISIF=2 for initial relaxation of ions while keeping cell fixed.
        b. Uses ISIF=3 for further relaxation of both ions and cell.
        c. Uses ISIF=3(s) with a smaller number of steps for final convergence.
    4. Outputs the final energy and band gap to a results file.
    """

    structure = Structure.from_file(cif)
    structure.comment = f"Structure {os.path.splitext(cif)[0]}"

    # Get dynamic parallelization settings
    parallel_settings = get_parallelization_settings()
    print(f"Using parallelization: {parallel_settings}")

    # for settings, https://github.com/materialsproject/pymatgen-io-validation/
    # as well as check for similar materials in Mpj
    # Materials Project compliant settings
    # Common settings for all VASP runs
    common_settings = dict(
        PREC="Accurate",    # Internal precision related params
        ENCUT=520,          # Energy cutoff for plane waves (size of basis set, scaling is expensive)
        EDIFF=1e-6,         # Energy convergence criterion
        EDIFFG=-0.05,       # Force convergence criterion (changed to negative for force-based)
        ISMEAR=0,          # smearing applied to electronic states (-5 = tetrahedron method)
        SIGMA=0.05,         # Width of smearing function
        IBRION=2,           # Ionic relaxation algorithm (2 = conjugate gradient)
        ISPIN=1,            # Spin polarization (1 = non-spin-polarized, 2 = if magnetic)
        LREAL=False,        # Whether to use real-space projection for faster calculations
        LWAVE=True,         # Write WAVECAR file (useful for restarting calculations)
        LCHARG=True,        # Write CHGCAR file (useful for charge density analysis)
        NELM=200,           # Max steps for electronic convergence
        ISYM=2,             # Symmetry handling (2 = automatic symmetry detection)
        SYMPREC=1e-5,       # Precision for symmetry detection
        **parallel_settings,  # Parallelization settings
        LASPH=True,         # Whether to use the ASR (Augmented Space Relaxation) method
        NELMIN=5,           # Minimum number of electronic steps before ionic relaxation (MP standard)
        ALGO="Normal",      # Algorithm for electronic minimization
        # IMAGES = 0, ICHARG = 1, NELMIN = 3 # potentially as these are on MPj?
    )

    isif2_settings  = {**common_settings, "ISIF": 2, "NSW": 200, "POTIM": 0.5}
    isif3_settings  = {**common_settings, "ISIF": 3, "NSW": 300, "POTIM": 0.75}
    isif3s_settings = {**common_settings, "ISIF": 3, "NSW": 10, "POTIM": 0.75}

    cyc = get_cycle()
    while cyc < 20:
        try:
            run_step(structure, isif2_settings, f"Cycle{cyc}_ISIF2")
            structure = Structure.from_file("CONTCAR")

            # Print intermediate results
            vr = Vasprun("vasprun.xml", parse_eigen=False)  # Skip eigenvalues for speed
            print(f"ISIF2 - Energy: {vr.final_energy:.6f} eV")

            run_step(structure, isif3_settings, f"Cycle{cyc}_ISIF3")
            structure = Structure.from_file("CONTCAR")

            # Print intermediate results
            vr = Vasprun("vasprun.xml", parse_eigen=False)
            print(f"ISIF3 - Energy: {vr.final_energy:.6f} eV")

            if check_finished(f"Cycle{cyc}_ISIF3.OUTCAR"):
                run_step(structure, isif3s_settings, f"Cycle{cyc}_ISIF3s")
                if check_finished(f"Cycle{cyc}_ISIF3s.OUTCAR"):
                    print("####Workflow converged successfully.####")
                    # Print final results
                    vr = Vasprun("vasprun.xml", parse_eigen=True)
                    print(f"ISIF3s Energy: {vr.final_energy:.6f} eV")
                    try:
                        band_gap = vr.get_band_structure().get_band_gap()['energy']
                        print(f"ISIF3s Band Gap: {band_gap:.4f} eV")
                    except:
                        print("Band gap calculation failed - likely metallic system")
                        band_gap = 0.0
                    break
            cyc += 1
        except Exception as e:
            print(f"Error in cycle {cyc}: {e}")
            if cyc == 0:  # If first cycle fails, exit
                raise e
            break
    
    # Final results
    vr = Vasprun("vasprun.xml", parse_eigen=True)
    try:
        gap = vr.get_band_structure().get_band_gap()["energy"]
    except:
        gap = 0.0  # Handle metallic systems or calculation errors
    return vr.final_energy, gap

def main():
    pr = argparse.ArgumentParser(description="Run VASP calculations with Materials Project compliance")
    pr.add_argument("cif_dir", help="Directory containing CIF files")
    args = pr.parse_args()

    # Validate environment
    potdir = os.environ.get("PMG_VASP_PSP_DIR")
    vasp_cmd = os.environ.get("VASP_COMMAND")
    
    if not potdir:
        sys.exit("Please set PMG_VASP_PSP_DIR environment variable")
    if not vasp_cmd:
        sys.exit("Please set VASP_COMMAND environment variable")
    if not os.path.exists(potdir):
        sys.exit(f"POTCAR directory not found: {potdir}")

    print(f"POTCAR directory: {potdir}")
    print(f"VASP command: {vasp_cmd}")
    print(f"Available processors (NSLOTS): {os.environ.get('NSLOTS', '1')}")

    output_dir = "vasp_outputs"
    os.makedirs(output_dir, exist_ok=True)

    # Process all CIF files
    cif_files = sorted([f for f in os.listdir(args.cif_dir) if f.lower().endswith(".cif")])
    print(f"Found {len(cif_files)} CIF files to process")

    for cif in cif_files:
        case = os.path.splitext(cif)[0]
        case_path = os.path.join(output_dir, case)
        
        print(f"\n{'='*60}")
        print(f"Processing: {case}")
        print(f"{'='*60}")
        
        if not os.path.isdir(case_path):
            os.makedirs(case_path)
            
        shutil.copy(os.path.join(args.cif_dir, cif), case_path)
        
        original_dir = os.getcwd()
        os.chdir(case_path)
        
        try:
            start_time = time.time()
            final_energy, band_gap = workflow(cif, potdir)
            end_time = time.time()
            duration = end_time - start_time
            
            # Write results
            with open("results.txt", "w") as f:
                f.write(f"Energy_eV {final_energy:.6f}\n")
                f.write(f"BandGap_eV {band_gap:.4f}\n")
                f.write(f"Runtime_s {duration:.2f}\n")
                f.write(f"Structure {case}\n")
            
            print(f"Completed {case} in {duration:.2f} seconds")
            
        except Exception as e:
            print(f"Failed to process {case}: {e}")
            with open("error.txt", "w") as f:
                f.write(f"Error: {e}\n")
        finally:
            os.chdir(original_dir)

if __name__ == "__main__":
    main()