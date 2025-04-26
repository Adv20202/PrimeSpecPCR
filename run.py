#!/usr/bin/env python3
"""
PrimeSpecPCR Setup and Run Script
This script checks for required dependencies and installs them if necessary
before launching the PrimeSpecPCR application.
"""

import os
import sys
import subprocess
import platform
import shutil
import time
import venv
import tempfile
from pathlib import Path

# Define required Python packages
REQUIRED_PACKAGES = [
    'tkinter',  # For GUI
    'biopython',  # For sequence handling
    'primer3-py',  # For primer design
    'pandas',  # For data processing
    'numpy',  # For numerical operations
    'tqdm',  # For progress bars
    'validators',  # For input validation
]

# Define required system tools
REQUIRED_TOOLS = [
    'mafft',  # For multiple sequence alignment
]

# Function to check if running in a virtual environment
def is_virtual_env():
    """Check if the script is running in a virtual environment."""
    return hasattr(sys, 'real_prefix') or (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix)

# Function to check if a Python package is installed
def is_package_installed(package_name):
    """Check if a Python package is installed."""
    if package_name == 'tkinter':
        try:
            import tkinter
            return True
        except ImportError:
            return False
    else:
        try:
            subprocess.check_call([sys.executable, '-m', 'pip', 'show', package_name],
                                stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)
            return True
        except subprocess.CalledProcessError:
            return False

# Function to check if a system tool is installed
# Function to check if a system tool is installed
def is_tool_installed(tool_name):
    """Check if a system tool is available in PATH or in local directories."""
    # For tools other than MAFFT, use standard PATH check
    if tool_name != 'mafft':
        return shutil.which(tool_name) is not None
    
    # Special handling for MAFFT
    # First check if MAFFT is in PATH
    if shutil.which('mafft') is not None:
        return True
    
    # Then check local directories based on OS
    current_dir = os.getcwd()
    mafft_dir = os.path.join(current_dir, "mafft")  # Main mafft directory
    system = platform.system()
    
    if system == "Darwin":  # macOS
        # Check for macOS specific MAFFT installation
        mac_mafft_bat = os.path.join(mafft_dir, "mafft-mac", "mafft.bat")
        local_mafft = os.path.join(current_dir, "mafft")  # Script in main directory
        
        if os.path.exists(mac_mafft_bat):
            return True
        elif os.path.exists(local_mafft) and os.access(local_mafft, os.X_OK):
            return True
    
    elif system == "Windows":
        # Check for Windows specific MAFFT installation
        win_mafft_bat = os.path.join(mafft_dir, "mafft-win", "mafft.bat")
        local_mafft_bat = os.path.join(current_dir, "mafft.bat")
        
        if os.path.exists(win_mafft_bat):
            return True
        elif os.path.exists(local_mafft_bat):
            return True
    
    elif system == "Linux":
        # Check for Linux specific MAFFT installation
        local_mafft = os.path.join(current_dir, "mafft")
        if os.path.exists(local_mafft) and os.access(local_mafft, os.X_OK):
            return True
        
        # Check standard Linux paths
        for path in ['/usr/bin/mafft', '/usr/local/bin/mafft']:
            if os.path.exists(path) and os.access(path, os.X_OK):
                return True
    
    # If we get here, MAFFT is not installed
    return False

# Function to install Python packages
def install_package(package_name):
    """Install a Python package using pip."""
    print(f"Installing {package_name}...")
    try:
        if package_name == 'tkinter':
            if platform.system() == 'Darwin':  # macOS
                # Try to install python-tk using homebrew
                try:
                    subprocess.check_call(['brew', 'install', 'python-tk'],
                                         stdout=subprocess.DEVNULL,
                                         stderr=subprocess.DEVNULL)
                    print(f"Installed {package_name} using Homebrew")
                    return True
                except (subprocess.CalledProcessError, FileNotFoundError):
                    print("Could not install tkinter using Homebrew")
                    # Fallback to instructions
                    print("\nTkinter is not installed. Please install it using one of these methods:")
                    print("  - Install via Homebrew: brew install python-tk")
                    print("  - Download Python from python.org (includes tkinter)")
                    return False
            elif platform.system() == 'Linux':
                # Try to detect the distro and install tkinter
                try:
                    if os.path.exists('/etc/debian_version'):  # Debian/Ubuntu
                        subprocess.check_call(['sudo', 'apt-get', 'update'],
                                            stdout=subprocess.DEVNULL,
                                            stderr=subprocess.DEVNULL)
                        subprocess.check_call(['sudo', 'apt-get', 'install', '-y', 'python3-tk'],
                                            stdout=subprocess.DEVNULL,
                                            stderr=subprocess.DEVNULL)
                        return True
                    elif os.path.exists('/etc/fedora-release'):  # Fedora
                        subprocess.check_call(['sudo', 'dnf', 'install', '-y', 'python3-tkinter'],
                                            stdout=subprocess.DEVNULL,
                                            stderr=subprocess.DEVNULL)
                        return True
                    elif os.path.exists('/etc/arch-release'):  # Arch
                        subprocess.check_call(['sudo', 'pacman', '-S', '--noconfirm', 'tk'],
                                            stdout=subprocess.DEVNULL,
                                            stderr=subprocess.DEVNULL)
                        return True
                    else:
                        print("\nCould not automatically install tkinter. Please install it manually.")
                        print("For example, on Ubuntu/Debian: sudo apt-get install python3-tk")
                        return False
                except (subprocess.CalledProcessError, FileNotFoundError):
                    print("\nCould not automatically install tkinter. Please install it manually.")
                    print("For example, on Ubuntu/Debian: sudo apt-get install python3-tk")
                    return False
            else:  # Windows or other
                print("\nTkinter should be included with Python on this platform.")
                print("If it's missing, consider reinstalling Python with the tcl/tk option selected.")
                return False
        else:
            # For other packages, use pip
            cmd = [sys.executable, '-m', 'pip', 'install', package_name]
            
            # Check if we need --user flag (when not in a venv and not using sudo)
            if not is_virtual_env() and os.geteuid() != 0 if hasattr(os, 'geteuid') else False:
                cmd.append('--user')
                
            # On macOS with system Python, might need --break-system-packages
            if platform.system() == 'Darwin' and not is_virtual_env():
                cmd.append('--break-system-packages')
                
            subprocess.check_call(cmd)
            return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to install {package_name}: {e}")
        return False

# Function to install system tools
def install_tool(tool_name):
    """Install a system tool using the appropriate package manager."""
    print(f"Installing {tool_name}...")
    try:
        if platform.system() == 'Darwin':  # macOS
            try:
                subprocess.check_call(['brew', 'install', tool_name],
                                     stdout=subprocess.DEVNULL,
                                     stderr=subprocess.DEVNULL)
                return True
            except (subprocess.CalledProcessError, FileNotFoundError):
                print(f"Could not install {tool_name} using Homebrew")
                print(f"Please install {tool_name} manually")
                return False
        elif platform.system() == 'Linux':
            try:
                if os.path.exists('/etc/debian_version'):  # Debian/Ubuntu
                    subprocess.check_call(['sudo', 'apt-get', 'update'],
                                        stdout=subprocess.DEVNULL,
                                        stderr=subprocess.DEVNULL)
                    subprocess.check_call(['sudo', 'apt-get', 'install', '-y', tool_name],
                                        stdout=subprocess.DEVNULL,
                                        stderr=subprocess.DEVNULL)
                    return True
                elif os.path.exists('/etc/fedora-release'):  # Fedora
                    subprocess.check_call(['sudo', 'dnf', 'install', '-y', tool_name],
                                        stdout=subprocess.DEVNULL,
                                        stderr=subprocess.DEVNULL)
                    return True
                elif os.path.exists('/etc/arch-release'):  # Arch
                    subprocess.check_call(['sudo', 'pacman', '-S', '--noconfirm', tool_name],
                                        stdout=subprocess.DEVNULL,
                                        stderr=subprocess.DEVNULL)
                    return True
                else:
                    print(f"\nCould not automatically install {tool_name}. Please install it manually.")
                    return False
            except (subprocess.CalledProcessError, FileNotFoundError):
                print(f"\nCould not automatically install {tool_name}. Please install it manually.")
                return False
        else:  # Windows or other
            print(f"\nAutomatic installation of {tool_name} is not supported on this platform.")
            print(f"Please install {tool_name} manually.")
            return False
    except Exception as e:
        print(f"Failed to install {tool_name}: {e}")
        return False

# Function to create a virtual environment
def create_virtual_environment():
    """Create a virtual environment for running the application."""
    print("Creating a virtual environment...")
    venv_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'venv')
    
    try:
        venv.create(venv_dir, with_pip=True)
        
        # Get path to Python in the virtual environment
        if platform.system() == "Windows":
            python_path = os.path.join(venv_dir, "Scripts", "python.exe")
        else:
            python_path = os.path.join(venv_dir, "bin", "python")
        
        # Upgrade pip
        subprocess.check_call([python_path, '-m', 'pip', 'install', '--upgrade', 'pip'],
                             stdout=subprocess.DEVNULL)
        
        # Install required packages in the virtual environment
        for package in REQUIRED_PACKAGES:
            if package != 'tkinter':  # Skip tkinter as it can't be installed via pip
                subprocess.check_call([python_path, '-m', 'pip', 'install', package],
                                    stdout=subprocess.DEVNULL)
        
        return venv_dir
    except Exception as e:
        print(f"Error creating virtual environment: {e}")
        return None

# Function to run PrimeSpecPCR
def run_application(use_venv=False, venv_dir=None):
    """Run the PrimeSpecPCR application."""
    # Determine the path to the main script
    primespecpcr_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PrimeSpecPCR.py')
    gui_application_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gui_application.py')
    
    # If PrimeSpecPCR.py doesn't exist but gui_application.py does, copy it
    if not os.path.exists(primespecpcr_script) and os.path.exists(gui_application_script):
        shutil.copy(gui_application_script, primespecpcr_script)
        print(f"Created {primespecpcr_script} from {gui_application_script}")
    
    if not os.path.exists(primespecpcr_script):
        print(f"Error: Could not find {primespecpcr_script}")
        return False
    
    try:
        if use_venv and venv_dir:
            # Run using the virtual environment Python
            if platform.system() == "Windows":
                python_path = os.path.join(venv_dir, "Scripts", "python.exe")
            else:
                python_path = os.path.join(venv_dir, "bin", "python")
            subprocess.check_call([python_path, primespecpcr_script])
        else:
            # Run using the current Python interpreter
            subprocess.check_call([sys.executable, primespecpcr_script])
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running application: {e}")
        return False

def main():
    """Main function to set up and run PrimeSpecPCR."""
    print("PrimeSpecPCR Setup and Run")
    print("==========================")
    
    # Check for required packages
    print("\nChecking required Python packages...")
    missing_packages = []
    for package in REQUIRED_PACKAGES:
        if not is_package_installed(package):
            print(f"  - {package}: Not installed")
            missing_packages.append(package)
        else:
            print(f"  - {package}: Installed")
    
    # Check for required tools
    print("\nChecking required system tools...")
    missing_tools = []
    for tool in REQUIRED_TOOLS:
        if not is_tool_installed(tool):
            print(f"  - {tool}: Not installed")
            missing_tools.append(tool)
        else:
            print(f"  - {tool}: Installed")
    
    # If there are missing dependencies, handle them
    if missing_packages or missing_tools:
        print("\nSome dependencies are missing.")
        
        # Try automatic installation if user agrees
        response = input("Do you want to automatically install missing dependencies? (y/n): ")
        if response.lower() == 'y':
            # Install missing packages
            package_install_failed = False
            for package in missing_packages:
                if not install_package(package):
                    package_install_failed = True
            
            # Install missing tools
            tool_install_failed = False
            for tool in missing_tools:
                if not install_tool(tool):
                    tool_install_failed = True
            
            # If any installation failed, offer virtual environment option
            if package_install_failed:
                print("\nSome Python packages could not be installed.")
                response = input("Do you want to create a virtual environment with the required packages? (y/n): ")
                if response.lower() == 'y':
                    venv_dir = create_virtual_environment()
                    if venv_dir:
                        print("\nVirtual environment created successfully.")
                        print("Running PrimeSpecPCR in the virtual environment...")
                        run_application(use_venv=True, venv_dir=venv_dir)
                        return
                    else:
                        print("Failed to create virtual environment.")
            
            # If some tool installations failed, notify user
            if tool_install_failed:
                print("\nSome system tools could not be installed automatically.")
                print("Please install them manually before continuing.")
                return
            
            # If all dependencies are now installed, run the application
            print("\nAll dependencies installed successfully.")
            run_application()
        else:
            print("Setup canceled by user.")
    else:
        print("\nAll dependencies are installed.")
        print("Running PrimeSpecPCR...")
        run_application()

if __name__ == "__main__":
    main()
