import hashlib
import json
import subprocess
from datetime import datetime

def file_hash(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()

def get_git_info(diff_file=None):
    try:
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip()
        dirty = subprocess.call(['git', 'diff-index', '--quiet', 'HEAD', '--']) != 0
        if dirty and diff_file:
            diff = subprocess.check_output(['git', 'diff']).decode()
            with open(diff_file, 'w') as f:
                f.write(diff)
        return {'commit': commit, 'dirty': dirty, 'diff_file': diff_file if dirty else None}
    except Exception:
        return {'commit': None, 'dirty': None, 'diff_file': None}

def get_conda_info():
    """Capture Conda environment details and package list."""
    try:
        env_info = subprocess.check_output(['conda', 'info', '--show']).decode()
        packages = subprocess.check_output(['conda', 'list']).decode()
        return {'conda_info': env_info, 'conda_list': packages}
    except Exception:
        return {'conda_info': None, 'conda_list': None}

def log_run_metadata(input_files, args, run_metadata_file='run_metadata.json'):
    start_time = datetime.utcnow().isoformat()
    run_metadata = {
        'start_time': start_time,
        'input_files': {f: file_hash(f) for f in input_files},
        'args': args,
        'git': get_git_info(diff_file='run_diff.patch'),
        'conda': get_conda_info()
    }

    with open(run_metadata_file, 'w') as f:
        json.dump(run_metadata, f, indent=2)

    return run_metadata_file, run_metadata

def finalise_run_metadata(run_metadata_file, run_metadata):
    run_metadata['end_time'] = datetime.utcnow().isoformat()
    with open(run_metadata_file, 'w') as f:
        json.dump(run_metadata, f, indent=2)
