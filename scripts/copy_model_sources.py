#!/usr/bin/env python3
"""
================================================================================
SCRIPT: copy_model_sources.py
PURPOSE:
    Safely populate the local 'modeles/' directory with clean source code for
    ORCHIDEE, XIOS, and IOIPSL from Toubkal's shared team installation.

DESCRIPTION:
    When cloning the project from GitHub, the large Fortran and C++ source trees
    for ORCHIDEE, XIOS, and IOIPSL are not included (they are excluded by
    .gitignore to keep the repository lightweight).
    
    This script safely copies the required upstream source files from the shared
    team repository into your personal workspace, ensuring that:
    1. Only genuine source code and build configs are copied.
    2. Local architecture files (e.g. GCC 12 configs in arch/) are PRESERVED and
       never overwritten.
    3. Legacy compilation artifacts (*.o, *.mod, *.a, build/ directories) are
       strictly EXCLUDED to avoid mixing incompatible compiler ABIs.
    4. Large NetCDF datasets (*.nc) and simulation run outputs are EXCLUDED.
    5. The shared team directory is accessed STRICTLY IN READ-ONLY mode.

DIRECTORIES:
    Source:      /srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles
    Destination: <your_workspace>/modipsl/modeles (e.g. modeles/IOIPSL, modeles/XIOS, modeles/ORCHIDEE)

COPIED COMPONENTS:
    - IOIPSL  : Fortran source files in src/, Makefile, and tool scripts.
    - XIOS    : C++ and Fortran sources in src/, extern/, inputs/, and make_xios.
    - ORCHIDEE: Fortran physics/driver sources (src_*), FCM build tools, and configs.

EXCLUDED ELEMENTS:
    - Object files (*.o), Fortran modules (*.mod), static libraries (*.a), executables (*.exe).
    - Previous build directories: build/, obj/, ppsrc/, lib/, inc/, bin/, flags/, done/.
    - Lock files: fcm.bld.lock.
    - Version control metadata: .svn/, .cache/.
    - NetCDF files: *.nc, *.nc4.

USAGE EXAMPLES:
    1. Preview mode (Dry-run, recommended first step):
       python3 scripts/copy_model_sources.py --dry-run

    2. Actual copy execution:
       python3 scripts/copy_model_sources.py

    3. Verbose preview (detailed file-by-file listing):
       python3 scripts/copy_model_sources.py --dry-run --verbose

AUTHOR:
    Mohammad El Aabaribaoune — IWRI / UM6P (Mohammed VI Polytechnic University)
    Date: November 2026
================================================================================
"""

import os
import sys
import shutil
import argparse
from pathlib import Path

# ==============================================================================
# CONFIGURATION & FILTERING RULES
# ==============================================================================

# Default path to the certified shared modeles directory on Toubkal HPC
DEFAULT_SHARED_DIR = "/srv/lustre01/project/climat-um6p-st-iwri-7ksifkvwkuy/shared/TEAM/models/modipsl/modeles"

# The three essential modeling components required for ORCHIDEE Offline
REQUIRED_COMPONENTS = ["IOIPSL", "XIOS", "ORCHIDEE"]

# Directory names to completely skip during traversal (prevents copying old builds)
EXCLUDED_DIR_NAMES = {
    "build", "build_gcc7", "bin", "bin_gcc7", "lib", "lib_gcc7_backup",
    "inc", "inc_gcc7_backup", "obj", "ppsrc", "flags", "done", ".cache",
    ".svn", ".void_dir", ".fcm", "__pycache__", "simulations"
}

# File extensions to skip (compiled binaries, libraries, data, and logs)
EXCLUDED_EXTENSIONS = {
    ".o", ".mod", ".a", ".so", ".exe", ".nc", ".nc4", ".log",
    ".out", ".err", ".swp", ".done"
}

# Exact file names to skip
EXCLUDED_FILENAMES = {
    "fcm.bld.lock", "dim2_driver.exe", "xios_server.exe", "orchidee_ol",
    "core", ".DS_Store"
}

# ==============================================================================
# COMMAND LINE ARGUMENTS
# ==============================================================================

def parse_args():
    """Parse user command-line arguments with sensible default paths."""
    # Automatically locate local 'modeles' directory relative to this script
    script_directory = Path(__file__).resolve().parent
    default_target_dir = (script_directory.parent / "modeles").resolve()

    parser = argparse.ArgumentParser(
        description="Safely populate ORCHIDEE/XIOS/IOIPSL source trees from shared Toubkal repository."
    )
    parser.add_argument(
        "--source-dir",
        default=DEFAULT_SHARED_DIR,
        help=f"Path to shared modeles source directory (default: {DEFAULT_SHARED_DIR})"
    )
    parser.add_argument(
        "--target-dir",
        default=str(default_target_dir),
        help=f"Path to local target modeles directory (default: {default_target_dir})"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview actions and statistics without copying or creating any files"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed file-by-file actions to terminal"
    )
    return parser.parse_args()

# ==============================================================================
# PRE-FLIGHT VALIDATION
# ==============================================================================

def validate_paths(source_dir: Path, target_dir: Path):
    """Ensure shared source directory exists, is readable, and contains components."""
    if not source_dir.is_dir():
        print(f"ERROR: Shared source directory not found: {source_dir}", file=sys.stderr)
        sys.exit(1)

    # Verify read permission on shared space
    if not os.access(source_dir, os.R_OK):
        print(f"ERROR: No read permission on source directory: {source_dir}", file=sys.stderr)
        sys.exit(1)

    # Check presence of all required component folders
    missing = [c for c in REQUIRED_COMPONENTS if not (source_dir / c).is_dir()]
    if missing:
        print("ERROR: Missing required components in shared source directory:", file=sys.stderr)
        for m in missing:
            print(f"  - {source_dir / m}", file=sys.stderr)
        sys.exit(1)

    # Verify target directory exists
    if not target_dir.is_dir():
        print(f"ERROR: Target modeles directory not found: {target_dir}", file=sys.stderr)
        sys.exit(1)

# ==============================================================================
# COPY COMPONENT LOGIC
# ==============================================================================

def copy_component(comp_name: str, src_comp: Path, dst_comp: Path, dry_run: bool, verbose: bool):
    """
    Safely copy a single component tree from source to destination.
    Preserves existing files (e.g. arch/) and skips unwanted build artifacts.
    """
    copied_files = 0
    copied_bytes = 0
    preserved_files = 0
    excluded_files = 0
    excluded_dirs = 0

    print(f"\n---> Processing component: {comp_name}")
    print(f"     Source : {src_comp}")
    print(f"     Target : {dst_comp}")

    # Traverse directory tree
    for current_root, dirs, files in os.walk(src_comp, followlinks=False):
        current_rel = Path(current_root).relative_to(src_comp)

        # ----------------------------------------------------------------------
        # 1. Filter subdirectories in-place:
        # Modifying 'dirs' in-place tells os.walk NOT to descend into excluded dirs.
        # This completely avoids scanning large build/ and obj/ folders.
        # ----------------------------------------------------------------------
        excluded_in_current = []
        for d in dirs:
            if d in EXCLUDED_DIR_NAMES or d.startswith(("build_", "bin_", "lib_", "inc_")):
                excluded_in_current.append(d)

        for d in excluded_in_current:
            dirs.remove(d)
            excluded_dirs += 1
            if verbose:
                print(f"  [EXCLUDE DIR] {current_rel / d}")

        # ----------------------------------------------------------------------
        # 2. Process each file in the current directory
        # ----------------------------------------------------------------------
        for f in files:
            src_file = Path(current_root) / f
            rel_file_path = current_rel / f
            dst_file = dst_comp / rel_file_path

            # Check if file should be excluded
            suffix = Path(f).suffix.lower()
            if (f in EXCLUDED_FILENAMES or 
                suffix in EXCLUDED_EXTENSIONS or 
                f.endswith("~") or 
                f.startswith(".#")):
                excluded_files += 1
                if verbose:
                    print(f"  [EXCLUDE FILE] {rel_file_path}")
                continue

            # Check if destination file already exists (e.g. GCC 12 architecture files)
            # Crucial: never overwrite existing files in your local workspace!
            if dst_file.exists():
                preserved_files += 1
                if verbose:
                    print(f"  [PRESERVED]    {rel_file_path}")
                continue

            # Record file to copy
            try:
                file_size = src_file.stat().st_size
            except OSError:
                file_size = 0

            copied_files += 1
            copied_bytes += file_size

            if verbose:
                print(f"  [COPY]         {rel_file_path} ({file_size / 1024:.1f} KB)")

            # Execute actual file copy only if NOT in dry-run mode
            if not dry_run:
                dst_file.parent.mkdir(parents=True, exist_ok=True)
                if src_file.is_symlink():
                    link_target = os.readlink(src_file)
                    try:
                        if dst_file.is_symlink() or dst_file.exists():
                            dst_file.unlink()
                        os.symlink(link_target, dst_file)
                    except OSError:
                        shutil.copy2(src_file, dst_file)
                else:
                    shutil.copy2(src_file, dst_file)

    return {
        "copied_files": copied_files,
        "copied_bytes": copied_bytes,
        "preserved_files": preserved_files,
        "excluded_files": excluded_files,
        "excluded_dirs": excluded_dirs
    }

# ==============================================================================
# MAIN EXECUTION FLOW
# ==============================================================================

def main():
    args = parse_args()
    source_dir = Path(args.source_dir).resolve()
    target_dir = Path(args.target_dir).resolve()

    mode_label = "DRY-RUN (PREVIEW ONLY)" if args.dry_run else "EXECUTE (ACTUAL COPY)"

    print("==============================================================================")
    print(" ORCHIDEE Offline — Model Sources Safe Staging Utility")
    print(f" Mode   : {mode_label}")
    print(f" Source : {source_dir}")
    print(f" Target : {target_dir}")
    print("==============================================================================")

    # 1. Validate paths and component directories
    validate_paths(source_dir, target_dir)

    # 2. Iterate through required components and stage sources
    total_copied_files = 0
    total_copied_bytes = 0
    total_preserved_files = 0
    total_excluded_files = 0
    total_excluded_dirs = 0

    results = {}
    for comp in REQUIRED_COMPONENTS:
        src_comp = source_dir / comp
        dst_comp = target_dir / comp
        res = copy_component(comp, src_comp, dst_comp, args.dry_run, args.verbose)
        results[comp] = res
        total_copied_files += res["copied_files"]
        total_copied_bytes += res["copied_bytes"]
        total_preserved_files += res["preserved_files"]
        total_excluded_files += res["excluded_files"]
        total_excluded_dirs += res["excluded_dirs"]

    # 3. Print clear, formatted summary table
    print("\n==============================================================================")
    print(f" SUMMARY OF ACTIONS ({mode_label})")
    print("==============================================================================")
    print(f"{'Component':12s} | {'Files to Copy':14s} | {'Volume (MB)':12s} | {'Preserved':10s} | {'Excluded Files':14s}")
    print("-" * 78)
    for comp, res in results.items():
        mb = res["copied_bytes"] / (1024 * 1024)
        print(f"{comp:12s} | {res['copied_files']:14d} | {mb:12.2f} | {res['preserved_files']:10d} | {res['excluded_files']:14d}")
    print("-" * 78)
    total_mb = total_copied_bytes / (1024 * 1024)
    print(f"{'TOTAL':12s} | {total_copied_files:14d} | {total_mb:12.2f} | {total_preserved_files:10d} | {total_excluded_files:14d}")
    print(f"\nExcluded Directories Count : {total_excluded_dirs}")
    print(f"Shared Directory State     : UNTOUCHED (Strictly Read-Only)")

    if args.dry_run:
        print("\n--> DRY-RUN COMPLETED SUCCESSFULLY.")
        print("    No files were created or modified. Run without --dry-run to perform actual copy:")
        print("      python3 scripts/copy_model_sources.py")
    else:
        print("\n--> COPY COMPLETED SUCCESSFULLY.")
        print("    Model source trees are now staged and ready for GCC 12 compilation.")
    print("==============================================================================")

if __name__ == "__main__":
    main()
