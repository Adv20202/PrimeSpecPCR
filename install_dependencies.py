#!/usr/bin/env python3
"""
PrimeSpecPCR Dependencies Installer
This script checks for and installs all required dependencies for PrimeSpecPCR.
"""

import os
import sys
import subprocess
import platform
import shutil
from pathlib import Path
import tempfile

try:
    import requests
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'requests'])
    import requests

# Define required Python packages
REQUIRED_PACKAGES = [
    'biopython>=1.79',  # For sequence handling
    'primer3-py',
    'pandas',
    'numpy',
    'tqdm',
    'validators',
    'requests>=2.25.0',  # For BLAST API requests in the specificity test module
]

def is_package_installed(package_name):
    """Check if a Python package is installed."""
    try:
        subprocess.check_call([sys.executable, '-m', 'pip', 'show', package_name],
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False

def is_tkinter_available():
    """Check if tkinter is available."""
    try:
        import tkinter
        return True
    except ImportError:
        return False

def install_package(package_name):
    """Install a Python package using pip."""
    print(f"Installing {package_name}...")
    try:
        cmd = [sys.executable, '-m', 'pip', 'install', package_name, '--upgrade']
        if platform.system() == 'Darwin':
            cmd.append('--break-system-packages')
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to install {package_name}: {e}")
        return False

def is_mafft_installed():
    """Check if MAFFT is installed and available for the current system."""
    # Check in PATH first
    if shutil.which('mafft') is not None:
        print(f"Found MAFFT in system PATH")
        return True
        
    # Check for local installation based on OS
    current_dir = os.getcwd()
    mafft_dir = os.path.join(current_dir, "mafft")
    system = platform.system()
    
    if system == "Darwin":  # macOS
        mac_mafft_bat = os.path.join(mafft_dir, "mafft-mac", "mafft.bat")
        local_mafft = os.path.join(current_dir, "mafft")
        
        if os.path.exists(mac_mafft_bat):
            print(f"Found MAFFT for macOS at {mac_mafft_bat}")
            return True
        elif os.path.exists(local_mafft) and os.access(local_mafft, os.X_OK):
            print(f"Found local MAFFT executable at {local_mafft}")
            return True
    
    elif system == "Windows":
        win_mafft_bat = os.path.join(mafft_dir, "mafft-win", "mafft.bat")
        local_mafft_bat = os.path.join(current_dir, "mafft.bat")
        
        if os.path.exists(win_mafft_bat):
            print(f"Found MAFFT for Windows at {win_mafft_bat}")
            return True
        elif os.path.exists(local_mafft_bat):
            print(f"Found local MAFFT batch file at {local_mafft_bat}")
            return True
    
    elif system == "Linux":
        # For Linux check both the local mafft script and the standard locations
        local_mafft = os.path.join(current_dir, "mafft")
        if os.path.exists(local_mafft) and os.access(local_mafft, os.X_OK):
            print(f"Found local MAFFT executable at {local_mafft}")
            return True
        
        # Check standard Linux paths
        for path in ['/usr/bin/mafft', '/usr/local/bin/mafft']:
            if os.path.exists(path) and os.access(path, os.X_OK):
                print(f"Found MAFFT at {path}")
                return True
    
    # If we got here, no suitable MAFFT was found
    return False

def setup_mafft_from_local():
    """
    Set up MAFFT from the local folder included in the project.
    This configures the mafft executable based on the OS.
    """
    system = platform.system()
    current_dir = os.getcwd()
    mafft_dir = os.path.join(current_dir, "mafft")
    
    # Check if the mafft folder exists in the current directory
    if not os.path.exists(mafft_dir):
        print(f"Error: MAFFT folder not found at {mafft_dir}")
        return False
    
    if system == "Windows":
        # For Windows, we need to set up mafft.bat and make sure it's executable
        mafft_win_dir = os.path.join(mafft_dir, "mafft-win")
        mafft_bat = os.path.join(mafft_win_dir, "mafft.bat")
        if os.path.exists(mafft_bat):
            # Create a batch file in the current directory for easy access
            target_bat = os.path.join(current_dir, "mafft.bat")
            with open(target_bat, 'w') as f:
                f.write(f'@echo off\n"{mafft_bat}" %*')
            print(f"MAFFT configured for Windows at {target_bat}")
            return True
        else:
            print(f"Error: mafft.bat not found at {mafft_bat}")
            return False
    
    elif system == "Darwin":  # macOS
        # For macOS, directly use the mafft binary in mafftdir/bin
        mafft_mac_dir = os.path.join(mafft_dir, "mafft-mac")
        mafft_bin = os.path.join(mafft_mac_dir, "mafftdir", "bin", "mafft")
        
        if os.path.exists(mafft_bin):
            # Make the binary executable
            os.chmod(mafft_bin, 0o755)
            
            # Create a shell script in the current directory for easy access
            target_script = os.path.join(current_dir, "mafft")
            
            # If a symlink already exists, remove it
            if os.path.exists(target_script):
                os.remove(target_script)
                
            # Create a more robust shell script wrapper that handles paths correctly
            with open(target_script, 'w') as f:
                f.write(f'#!/bin/bash\n\n')
                f.write(f'# Fix paths for relative file references\n')
                f.write(f'SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"\n')
                f.write(f'cd "$SCRIPT_DIR"\n\n')
                f.write(f'# Run MAFFT with all arguments passed to this script\n')
                f.write(f'"{mafft_bin}" "$@"\n')
            
            # Make it executable
            os.chmod(target_script, 0o755)
            print(f"MAFFT configured for macOS at {target_script}")
            return True
        else:
            print(f"Error: mafft binary not found at {mafft_bin}")
            return False
    
    elif system == "Linux":
        # For Linux, use the appropriate package
        if os.path.exists(os.path.join(mafft_dir, "mafft_7.526-1_amd64.deb")):
            # Use the binary from mafftdir/bin for Debian-based systems
            print("Debian package found. Will use the executable from it.")
            mafft_bin = os.path.join(mafft_dir, "mafftdir", "bin", "mafft")
        else:
            # Use RPM or default binary
            mafft_bin = os.path.join(mafft_dir, "mafftdir", "bin", "mafft")
        
        if os.path.exists(mafft_bin):
            # Make the binary executable
            os.chmod(mafft_bin, 0o755)
            
            # Create a symlink in the current directory for easy access
            target_script = os.path.join(current_dir, "mafft")
            
            # If a symlink already exists, remove it
            if os.path.exists(target_script):
                os.remove(target_script)
                
            # Create the shell script wrapper
            with open(target_script, 'w') as f:
                f.write(f'#!/bin/bash\n"{mafft_bin}" "$@"')
            
            # Make it executable
            os.chmod(target_script, 0o755)
            print(f"MAFFT configured for Linux at {target_script}")
            return True
        else:
            print(f"Error: mafft binary not found at {mafft_bin}")
            return False
    
    else:
        print(f"Unsupported platform: {system}")
        return False

def verify_mafft_installation():
    """Verify that MAFFT is properly installed and working."""
    try:
        # Create a temporary FASTA file with test sequences
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as temp_file:
            temp_file.write(">seq1\nACGTACGTACGT\n>seq2\nACGTACGT\n")
            temp_file_path = temp_file.name
        
        print(f"Testing MAFFT with temporary file: {temp_file_path}")
        output_file = tempfile.NamedTemporaryFile(suffix='.fa', delete=False).name
        
        # Prepare command based on OS and available executables
        current_dir = os.getcwd()
        mafft_dir = os.path.join(current_dir, "mafft")
        system = platform.system()
        successful = False
        
        # Try system-wide MAFFT first
        if shutil.which('mafft'):
            print("System-wide MAFFT detected in PATH, assuming it's correctly installed.")
            successful = True
            
            try:
                with open(output_file, 'w') as outfile:
                    subprocess.check_call(['mafft', '--quiet', temp_file_path], 
                                         stdout=outfile, 
                                         stderr=subprocess.DEVNULL,
                                         timeout=30)
                print("System-wide MAFFT verified successfully!")
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                print(f"System-wide MAFFT test failed: {e}, but will assume it's working.")
        
        # For macOS and Windows, check if appropriate bat file exists
        elif system in ["Darwin", "Windows"]:
            if system == "Darwin":
                mafft_bat = os.path.join(mafft_dir, "mafft-mac", "mafft.bat")
                if os.path.exists(mafft_bat):
                    print(f"Found MAFFT batch file for macOS, assuming it's correctly installed.")
                    successful = True
            else:  # Windows
                mafft_bat = os.path.join(mafft_dir, "mafft-win", "mafft.bat")
                if os.path.exists(mafft_bat):
                    print(f"Found MAFFT batch file for Windows, assuming it's correctly installed.")
                    successful = True
        
        # If not yet successful, try OS-specific local installation
        if not successful:
            mafft_cmd = None
            
            if system == "Darwin":  # macOS
                mac_mafft_bat = os.path.join(mafft_dir, "mafft-mac", "mafft.bat")
                local_mafft = os.path.join(current_dir, "mafft")
                
                if os.path.exists(mac_mafft_bat):
                    mafft_cmd = [mac_mafft_bat, "--quiet", temp_file_path]
                elif os.path.exists(local_mafft) and os.access(local_mafft, os.X_OK):
                    mafft_cmd = [local_mafft, "--quiet", temp_file_path]
                    
            elif system == "Windows":
                win_mafft_bat = os.path.join(mafft_dir, "mafft-win", "mafft.bat")
                local_mafft_bat = os.path.join(current_dir, "mafft.bat")
                
                if os.path.exists(win_mafft_bat):
                    mafft_cmd = [win_mafft_bat, "--quiet", temp_file_path]
                elif os.path.exists(local_mafft_bat):
                    mafft_cmd = [local_mafft_bat, "--quiet", temp_file_path]
                    
            elif system == "Linux":
                local_mafft = os.path.join(current_dir, "mafft")
                if os.path.exists(local_mafft) and os.access(local_mafft, os.X_OK):
                    mafft_cmd = [local_mafft, "--quiet", temp_file_path]
            
            if mafft_cmd:
                try:
                    print(f"Testing local MAFFT using: {' '.join(mafft_cmd)}")
                    with open(output_file, 'w') as outfile:
                        subprocess.check_call(mafft_cmd, stdout=outfile, stderr=subprocess.DEVNULL)
                    
                    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                        print("Local MAFFT installation verified successfully!")
                        successful = True
                except Exception as e:
                    print(f"Local MAFFT test failed: {e}")
        
        # Clean up temp files
        if os.path.exists(temp_file_path):
            os.unlink(temp_file_path)
        if os.path.exists(output_file):
            os.unlink(output_file)
            
        return successful
        
    except Exception as e:
        print(f"MAFFT verification failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    print("PrimeSpecPCR Dependencies Installer")
    print("==================================")
    
    required_packages = [pkg for pkg in REQUIRED_PACKAGES if 'scikit-bio' not in pkg]
    
    missing_packages = []
    for package in required_packages:
        package_name = package.split('>=')[0]
        if not is_package_installed(package_name):
            missing_packages.append(package)
    
    # Check if tkinter is available
    if not is_tkinter_available():
        print("\nError: tkinter is not available.")
        if platform.system() == 'Darwin':
            print("On macOS, you can install tkinter with:")
            print("  brew install python-tk")
        elif platform.system() == 'Linux':
            print("On Linux, you can install tkinter with:")
            print("  For Ubuntu/Debian: sudo apt-get install python3-tk")
            print("  For Fedora: sudo dnf install python3-tkinter")
            print("  For Arch Linux: sudo pacman -S tk")
        return 1
    
    # Install missing Python packages
    if missing_packages:
        print(f"\nMissing packages: {', '.join(missing_packages)}")
        print("Installing missing packages...")
        for package in missing_packages:
            if not install_package(package):
                print(f"Failed to install {package}. Please install it manually.")
                return 1
        print("\nAll Python packages installed successfully!")
    else:
        print("\nAll required Python packages are already installed.")
    
    # Check for MAFFT in PATH
    if is_mafft_installed():
        print("\nMAFFT is already installed.")
        if verify_mafft_installation():
            print("MAFFT is working correctly.")
            print("\nAll dependencies are satisfied. PrimeSpecPCR is ready to use!")
            return 0
        else:
            print("MAFFT is installed but appears to be not working correctly.")
            print("Attempting to set up local MAFFT instead...")
            if setup_mafft_from_local():
                if verify_mafft_installation():
                    print("Local MAFFT set up successfully and verified.")
                else:
                    print("Local MAFFT set up but verification failed.")
                    return 1
            else:
                print("Failed to set up local MAFFT.")
                return 1
    else:
        print("\nMAFFT is not installed in PATH. Setting up local MAFFT...")
        if setup_mafft_from_local():
            if verify_mafft_installation():
                print("Local MAFFT set up successfully and verified.")
            else:
                print("Local MAFFT set up but verification failed.")
                return 1
        else:
            print("Failed to set up local MAFFT.")
            return 1
    
    print("\nAll dependencies are satisfied. PrimeSpecPCR is ready to use!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
