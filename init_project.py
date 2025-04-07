import subprocess
import sys
from pathlib import Path

def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}")
        print(f"Error: {e}")
        sys.exit(1)

def main():
    # Create necessary directories
    Path("data").mkdir(exist_ok=True)
    Path("data/models").mkdir(exist_ok=True)
    
    # Create virtual environment
    print("Creating virtual environment...")
    run_command("python -m venv venv")
    
    # Activate virtual environment and install dependencies
    if sys.platform == "win32":
        run_command("venv\\Scripts\\activate && pip install -r requirements.txt")
    else:
        run_command("source venv/bin/activate && pip install -r requirements.txt")
    
    # Install frontend dependencies
    print("\nInstalling frontend dependencies...")
    run_command("cd frontend && npm install")
    
    # Generate dataset and train model
    print("\nGenerating dataset and training model...")
    if sys.platform == "win32":
        run_command("venv\\Scripts\\python generate_data_and_train.py")
    else:
        run_command("./venv/bin/python generate_data_and_train.py")
    
    print("\nInitialization complete! You can now run the project using:")
    if sys.platform == "win32":
        print("run.bat")
    else:
        print("./run.sh")

if __name__ == "__main__":
    main() 