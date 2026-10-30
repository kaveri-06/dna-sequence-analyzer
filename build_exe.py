"""
DNA Sequence Analyzer - Executable Build Script
==============================================
Handles PyInstaller bundling with automated support for Python 3.14+
Tcl 9 virtual filesystem extraction and dependency packaging.
"""

import os
import shutil
import subprocess
import sys
import zipfile


def prepare_tcl_tk_assets() -> tuple[str, str]:
    """
    Extracts Tcl/Tk data files from Python 3.14+ zipfs archives
    so PyInstaller can bundle them into _tcl_data and _tk_data.
    """
    tcl_dir = os.path.join(sys.prefix, "tcl")
    bundle_dir = os.path.abspath("tcl_bundle")
    tcl_out = os.path.join(bundle_dir, "_tcl_data")
    tk_out = os.path.join(bundle_dir, "_tk_data")
    
    os.makedirs(tcl_out, exist_ok=True)
    os.makedirs(tk_out, exist_ok=True)
    
    # 1. Check for Python 3.14 zip-based Tcl
    z_tcl_path = os.path.join(tcl_dir, "libtcl9.0.4.zip")
    if not os.path.exists(z_tcl_path):
        # Search for any libtcl*.zip
        for f in os.listdir(tcl_dir) if os.path.exists(tcl_dir) else []:
            if f.startswith("libtcl") and f.endswith(".zip"):
                z_tcl_path = os.path.join(tcl_dir, f)
                break
                
    if os.path.exists(z_tcl_path):
        print(f"Extracting Tcl data from {os.path.basename(z_tcl_path)}...")
        with zipfile.ZipFile(z_tcl_path) as z:
            for member in z.namelist():
                if member.startswith("tcl_library/"):
                    rel_path = member[len("tcl_library/"):]
                    if rel_path:
                        target = os.path.join(tcl_out, rel_path)
                        if member.endswith("/"):
                            os.makedirs(target, exist_ok=True)
                        else:
                            os.makedirs(os.path.dirname(target), exist_ok=True)
                            with open(target, "wb") as f:
                                f.write(z.read(member))
        shutil.copy(z_tcl_path, tcl_out)
        
    # 2. Check for Python 3.14 zip-based Tk
    z_tk_path = os.path.join(tcl_dir, "libtk9.0.4.zip")
    if not os.path.exists(z_tk_path):
        for f in os.listdir(tcl_dir) if os.path.exists(tcl_dir) else []:
            if f.startswith("libtk") and f.endswith(".zip"):
                z_tk_path = os.path.join(tcl_dir, f)
                break
                
    if os.path.exists(z_tk_path):
        print(f"Extracting Tk data from {os.path.basename(z_tk_path)}...")
        with zipfile.ZipFile(z_tk_path) as z:
            for member in z.namelist():
                if member.startswith("tk_library/"):
                    rel_path = member[len("tk_library/"):]
                    if rel_path:
                        target = os.path.join(tk_out, rel_path)
                        if member.endswith("/"):
                            os.makedirs(target, exist_ok=True)
                        else:
                            os.makedirs(os.path.dirname(target), exist_ok=True)
                            with open(target, "wb") as f:
                                f.write(z.read(member))
        shutil.copy(z_tk_path, tk_out)
        
    # Also copy tk9.0 directory if exists
    tk9_dir = os.path.join(tcl_dir, "tk9.0")
    if os.path.exists(tk9_dir):
        shutil.copytree(tk9_dir, os.path.join(tk_out, "tk9.0"), dirs_exist_ok=True)
        
    return tcl_out, tk_out


def build():
    print("=" * 60)
    print("       Building Standalone DNA Sequence Analyzer.exe")
    print("=" * 60)
    
    tcl_out, tk_out = prepare_tcl_tk_assets()
    
    cmd = [
        sys.executable,
        "-m", "PyInstaller",
        "--onefile",
        "--windowed",
        "--name", "DNA Sequence Analyzer",
        "--clean",
        "--add-data", f"{tcl_out};_tcl_data",
        "--add-data", f"{tk_out};_tk_data",
        "app.py"
    ]
    
    print("\nRunning PyInstaller packaging...")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"[ERROR] PyInstaller failed with code {result.returncode}")
        sys.exit(result.returncode)
        
    dist_exe = os.path.join("dist", "DNA Sequence Analyzer.exe")
    target_exe = "DNA Sequence Analyzer.exe"
    if os.path.exists(dist_exe):
        shutil.copy2(dist_exe, target_exe)
        print(f"\n[SUCCESS] Executable copied to root: {os.path.abspath(target_exe)}")
    else:
        print(f"[ERROR] Could not find built executable at: {dist_exe}")
        sys.exit(1)


if __name__ == "__main__":
    build()
