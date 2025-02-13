import os
import subprocess
import yaml

# Automatically find the repo root (cross-platform safe)
def get_repo_root():
    """Find the root directory of the git repository."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True, text=True, check=True
        )
        return os.path.abspath(result.stdout.strip())
    except subprocess.CalledProcessError:
        print("⚠️ Not inside a git repository. Using the current working directory.")
        return os.path.abspath(os.getcwd())

# Set the repo root as a global variable
REPO_ROOT = get_repo_root()

# Path to the environment file (relative to repo root)
env_file_path = os.path.join(REPO_ROOT, "configs", "auto-prot-env.yml")

def get_active_env():
    """Returns the name of the currently active conda environment."""
    try:
        return os.environ['CONDA_DEFAULT_ENV']
    except KeyError:
        print("No active conda environment detected.")
        return None

def load_yaml_env(file_path):
    """Loads the dependencies from the specified YAML file."""
    with open(file_path, 'r') as file:
        env_data = yaml.safe_load(file)
        return set(env_data.get('dependencies', []))

def get_installed_packages():
    """Returns a set of installed packages with versions in the current conda environment."""
    result = subprocess.run(
        ["conda", "list", "--export"], 
        capture_output=True, text=True, check=True
    )
    packages = set()
    for line in result.stdout.splitlines():
        if not line.startswith("#") and "=" in line:
            # Capture both package name and version
            package, version = line.split('=')[:2]
            packages.add(f"{package}={version}")
    return packages

def compare_envs():
    """Compares the active conda environment with the reference YAML."""
    active_env = get_active_env()
    if not active_env:
        return
    print(f"📁 Using repo root: {REPO_ROOT}")
    print(f"Active conda environment: {active_env}")
    yaml_packages = load_yaml_env(env_file_path)
    installed_packages = get_installed_packages()
    # Compare packages
    if yaml_packages == installed_packages:
        print("✅ The active environment matches the YAML file.")
    else:
        print("⚠️ WARNING: The active environment differs from the YAML file! Consider updating the YAML")

if __name__ == "__main__":
    compare_envs()
