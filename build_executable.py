import os
import sys
import platform
import subprocess
import shutil
from pathlib import Path
import venv
import tempfile

def check_pyinstaller_in_env(env_dir):
    """Check if PyInstaller is installed in the virtual environment."""
    python_path = get_venv_python_path(env_dir)
    try:
        subprocess.check_call([python_path, '-m', 'pip', 'show', 'pyinstaller'], 
                             stdout=subprocess.DEVNULL, 
                             stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False

def get_venv_python_path(env_dir):
    """Get the Python executable path in the virtual environment."""
    if platform.system() == "Windows":
        return os.path.join(env_dir, "Scripts", "python.exe")
    else:
        return os.path.join(env_dir, "bin", "python")

def create_virtual_environment():
    """Create a virtual environment for building the executable."""
    print("Creating a virtual environment for building...")
    temp_dir = tempfile.mkdtemp(prefix="primespecpcr_build_")
    
    try:
        venv.create(temp_dir, with_pip=True)
        python_path = get_venv_python_path(temp_dir)
        
        # Upgrade pip
        subprocess.check_call([python_path, '-m', 'pip', 'install', '--upgrade', 'pip'])
        
        # Install PyInstaller
        print("Installing PyInstaller in virtual environment...")
        subprocess.check_call([python_path, '-m', 'pip', 'install', 'pyinstaller'])
        
        return temp_dir
    except Exception as e:
        print(f"Error creating virtual environment: {e}")
        shutil.rmtree(temp_dir, ignore_errors=True)
        return None

def build_executable(env_dir):
    """Build executable file using PyInstaller."""
    system = platform.system()
    
    print(f"Building application for {system}...")
    
    # Executable file name
    app_name = "PrimeSpecPCR"
    if system == "Windows":
        app_name += ".exe"
    
    # Path to main script
    main_script = "PrimeSpecPCR.py"
    
    # Rename GUI script if needed
    if not os.path.exists(main_script) and os.path.exists("gui_application.py"):
        shutil.copy("gui_application.py", main_script)
        print(f"Renamed gui_application.py to {main_script}")
    
    # Check if script exists
    if not os.path.exists(main_script):
        print(f"File {main_script} not found!")
        return False
    
    # List of files to include
    files_to_include = [
        "1_Genetic_Background.py",
        "2_MSA_Alignment.py",
        "3_PCR_Primers_Design.py",
        "4_Primers_Specificity.py"
    ]
    
    # Check if all files exist
    missing_files = [f for f in files_to_include if not os.path.exists(f)]
    if missing_files:
        print(f"Missing files: {', '.join(missing_files)}")
        return False
    
    # Create dist directory if it doesn't exist
    dist_dir = Path("dist")
    dist_dir.mkdir(exist_ok=True)
    
    # Prepare PyInstaller arguments
    python_path = get_venv_python_path(env_dir)
    
    # Data files separator and format based on platform
    separator = ";" if system == "Windows" else ":"
    
    # Build individual --add-data arguments for each file
    data_args = []
    for file in files_to_include:
        if system == "Windows":
            data_args.extend(["--add-data", f"{file}{separator}."])
        else:
            data_args.extend(["--add-data", f"{file}{separator}."])
    
    pyinstaller_args = [
        python_path, 
        '-m', 
        'PyInstaller',
        '--onefile',
        '--windowed',
        '--name', app_name,
    ]
    
    # Add data arguments
    pyinstaller_args.extend(data_args)
    
    # Add icon for Windows
    if system == "Windows":
        pyinstaller_args.extend(['--icon', 'NONE'])
    
    # Add main script
    pyinstaller_args.append(main_script)
    
    try:
        # Run PyInstaller
        subprocess.check_call(pyinstaller_args)
        
        # Move executable to main directory
        executable_path = os.path.join("dist", app_name)
        if os.path.exists(executable_path):
            shutil.copy(executable_path, ".")
            print(f"Executable created: {app_name}")
            return True
        else:
            print("Executable file not found after building.")
            return False
    except subprocess.CalledProcessError as e:
        print(f"Error building executable: {e}")
        return False

def main():
    print("PrimeSpecPCR Build Tool")
    print("======================")
    
    if platform.system() == "Darwin":  # macOS
        print("Detected macOS - using virtual environment approach to avoid system limitations.")
        env_dir = create_virtual_environment()
        if not env_dir:
            print("Cannot proceed without virtual environment.")
            return
            
        success = build_executable(env_dir)
        
        # Clean up
        print("Cleaning up temporary files...")
        shutil.rmtree(env_dir, ignore_errors=True)
    else:
        # For Windows and Linux, try a simpler approach first
        try:
            # Check if PyInstaller is already installed
            subprocess.check_call([sys.executable, '-m', 'pip', 'show', 'pyinstaller'], 
                                stdout=subprocess.DEVNULL, 
                                stderr=subprocess.DEVNULL)
            print("PyInstaller is already installed.")
            success = build_executable(".")
        except subprocess.CalledProcessError:
            # PyInstaller not installed, use virtual environment approach
            print("PyInstaller not found. Using virtual environment approach.")
            env_dir = create_virtual_environment()
            if not env_dir:
                print("Cannot proceed without virtual environment.")
                return
                
            success = build_executable(env_dir)
            
            # Clean up
            print("Cleaning up temporary files...")
            shutil.rmtree(env_dir, ignore_errors=True)
    
    if success:
        print("\nApplication built successfully!")
        print(f"You can now run PrimeSpecPCR by clicking on the {'PrimeSpecPCR.exe' if platform.system() == 'Windows' else 'PrimeSpecPCR'} file.")
    else:
        print("\nError occurred during build process.")
        print("Alternative options:")
        print("1. Install PyInstaller manually: pip install pyinstaller")
        print("2. Run the application directly: python PrimeSpecPCR.py")

if __name__ == "__main__":
    main()
