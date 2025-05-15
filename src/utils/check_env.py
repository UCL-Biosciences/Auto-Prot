import os
import platform
import subprocess

import yaml


# Automatically find the repo root (cross-platform safe)
def get_repo_root():
    """
    Find the root directory of the current git repository.

    Returns:
        str: Absolute path to the repository root, or current working directory if not in a git repo.
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True,
            text=True,
            check=True,
        )
        return os.path.abspath(result.stdout.strip())
    except subprocess.CalledProcessError:
        print("⚠️ Not inside a git repository. Using the current working directory.")
        return os.path.abspath(os.getcwd())


# Set the repo root as a global variable
REPO_ROOT = get_repo_root()

# which environment to use depends on the operating system
Active_OS = platform.system()

if Active_OS == "Windows":
    env_file_path = os.path.join(REPO_ROOT, "configs", "auto-prot-env-windowsOS.yml")
elif Active_OS == "Darwin":  # MAC is detected as Darwin
    env_file_path = os.path.join(REPO_ROOT, "configs", "auto-prot-env-macOS.yml")


# Path to the environment file (relative to repo root)


def get_active_env():
    """
    Get the name of the currently active conda environment.

    Returns:
        str or None: Name of the active environment, or None if not inside one.
    """
    try:
        return os.environ["CONDA_DEFAULT_ENV"]
    except KeyError:
        print("No active conda environment detected.")
        return None


def load_yaml_env(file_path):
    """
    Load dependencies from a Conda YAML environment file.

    Args:
        file_path (str): Path to the environment YAML file.

    Returns:
        set[str]: A set of 'package=version' strings for all dependencies, including pip entries.
    """
    with open(file_path) as file:
        env_data = yaml.safe_load(file)
        deps = env_data.get("dependencies", [])
        all_packages = set()
        # return set(env_data.get('dependencies', []))
        for dep in deps:
            if isinstance(dep, str):
                all_packages.add(dep)
            elif isinstance(dep, dict) and "pip" in dep:
                for pip_dep in dep["pip"]:
                    # Convert 'package==version' to 'package=version' for comparison
                    normalized = pip_dep.replace("==", "=")
                    all_packages.add(normalized)
    return all_packages


def get_installed_packages():
    """
    List installed packages in the current conda environment.

    Returns:
        set[str]: A set of 'package=version' strings for all currently installed packages.
    """
    result = subprocess.run(
        ["conda", "list", "--export"], capture_output=True, text=True, check=True
    )
    packages = set()
    for line in result.stdout.splitlines():
        if not line.startswith("#") and "=" in line:
            # Capture both package name and version
            package, version = line.split("=")[:2]
            packages.add(f"{package}={version}")
    return packages


def compare_envs():
    """
    Compare the active conda environment against a reference YAML file.

    Side effects:
        Prints messages indicating whether the environments match or differ.
    """
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
        print(
            "⚠️ WARNING: The active environment differs from the YAML file! Consider updating the YAML"
        )


if __name__ == "__main__":
    compare_envs()
