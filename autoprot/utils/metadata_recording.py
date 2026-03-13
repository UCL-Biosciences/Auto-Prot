import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from typing import Optional

import yaml

try:
    import psutil
    _PSUTIL_AVAILABLE = True
except ImportError:
    _PSUTIL_AVAILABLE = False


def file_hash(path: str) -> str:
    """Compute the SHA-256 hex digest of a file.

    Args:
        path: Path to the file to hash.

    Returns:
        Lowercase hex string of the SHA-256 digest.
    """
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()


def get_git_info(diff_file: Optional[str] = None) -> dict:
    """Capture git repository state at the time of the run.

    If the working tree is dirty and `diff_file` is provided, the full
    unstaged diff is written to that path so the exact code state can be
    reproduced later.

    Args:
        diff_file: Path to write the diff patch to when the repo is dirty.
                   If None, no patch is written.

    Returns:
        Dict with keys:
            commit    – full SHA of HEAD (str or None on failure)
            tag       – nearest tag from `git describe --tags --always` (str or None)
            dirty     – True if the working tree has uncommitted changes (bool or None)
            diff_file – path to the written patch file, or None if not applicable
    """
    try:
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip()
        dirty = subprocess.call(['git', 'diff-index', '--quiet', 'HEAD', '--']) != 0
        if dirty and diff_file:
            diff = subprocess.check_output(['git', 'diff']).decode()
            with open(diff_file, 'w') as f:
                f.write(diff)
        try:
            tag = subprocess.check_output(
                ['git', 'describe', '--tags' ],
                stderr=subprocess.DEVNULL
            ).decode().strip()
        except Exception:
            tag = None
        return {
            'commit': commit,
            'tag': tag,
            'dirty': dirty,
            'diff_file': diff_file if dirty else None,
        }
    except Exception:
        return {'commit': None, 'tag': None, 'dirty': None, 'diff_file': None}


def get_conda_env(env_name: Optional[str] = None) -> Optional[str]:
    """Return the raw output of `conda list` for a given environment.

    Args:
        env_name: Name of the conda environment to query. If None, the
                  currently active environment is used. For main.py to be running, auto-prot will be active so conda env details will be recorded.

    Returns:
        Raw text output of `conda list`, or None if the command fails.
    """
    try:
        cmd = ['conda', 'list']
        if env_name:
            cmd += ['-n', env_name]
        return subprocess.check_output(cmd).decode()
    except Exception:
        return None


def _utcnow() -> str:
    """Return the current UTC time as an ISO-8601 string."""
    return datetime.now(timezone.utc).isoformat()


def _write_yaml(run_meta: dict, run_metadata_file: str) -> None:
    """Serialise run_meta to YAML and write to disk, preserving key order."""
    with open(run_metadata_file, 'w') as f:
        yaml.dump(run_meta, f, default_flow_style=False, allow_unicode=True, sort_keys=False)


def log_run_metadata(
    input_files: list,
    args: list,
    config: dict,
    run_metadata_file: str = 'run_metadata.yaml',
    repo_root: Optional[str] = None,
) -> tuple:
    """Initialise and write the run metadata file at pipeline start.

    Collects input file hashes, command-line arguments, config contents,
    git state, conda package lists, and platform info, then writes them
    to a YAML file under <meta_dir>/run_metadata.yaml. Conda environment
    dumps are written to <meta_dir>/conda_envs/<env_name>.txt (one file
    per environment) and referenced by name only in the main metadata.
    If the git repo is dirty, a diff patch is written to
    <meta_dir>/run_diff.patch and referenced in the metadata.

    Args:
        input_files:       List of input file paths to hash.
        args:              Command-line arguments (typically sys.argv).
        config:            Full parsed config dictionary for this run.
        run_metadata_file: Path to write the YAML metadata file to.
        repo_root:         Repository root used to compute relative paths for
                           input_files. If None, absolute paths are stored.

    Returns:
        Tuple of (run_metadata_file, run_meta) where run_meta is the
        dict that was serialised — pass both to `record_step_complete`
        and `finalise_run_metadata`.
    """
    meta_dir = os.path.dirname(os.path.abspath(run_metadata_file))
    os.makedirs(meta_dir, exist_ok=True)

    # Write git diff inside the run_metadata dir
    diff_file = os.path.join(meta_dir, 'run_diff.patch')
    git_info = get_git_info(diff_file=diff_file)
    if git_info.get('diff_file'):
        git_info['diff_file'] = os.path.relpath(git_info['diff_file'], meta_dir)

    # Write conda env dumps to separate files; record names only in main metadata
    conda_env_names = []
    conda_envs = {
        'auto-proteomics': get_conda_env(),
        'r-limma-env': get_conda_env('r-limma-env'),
    }
    if any(v is not None for v in conda_envs.values()):
        conda_envs_dir = os.path.join(meta_dir, 'conda_envs')
        os.makedirs(conda_envs_dir, exist_ok=True)
        for env_name, env_text in conda_envs.items():
            if env_text is not None:
                conda_env_names.append(env_name)
                with open(os.path.join(conda_envs_dir, f'{env_name}.txt'), 'w') as f:
                    f.write(env_text)

    # Use relative input file paths (relative to repo_root when provided)
    if repo_root:
        input_files_rel = {
            os.path.relpath(f, repo_root).replace('\\', '/'): file_hash(f)
            for f in input_files
        }
    else:
        input_files_rel = {f: file_hash(f) for f in input_files}

    run_meta = {
        'run': {
            'exit_status': None,
            'start_time': _utcnow(),
            'end_time': None,
            'duration_s': None,
            'args': list(args),
        },
        'steps': [],
        'config': config,
        'provenance': {
            'platform': platform.platform(),
            'git': git_info,
            'conda_env_names': conda_env_names,
            'input_files': input_files_rel,
        },
    }

    _write_yaml(run_meta, run_metadata_file)
    return run_metadata_file, run_meta


def load_run_metadata(run_metadata_file: str) -> tuple:
    """Load an existing run metadata file for resume.

    Supports both YAML (.yaml/.yml) and legacy JSON (.json) formats,
    selected automatically by file extension.

    Args:
        run_metadata_file: Path to the metadata file.

    Returns:
        Tuple of (run_meta dict, completed_steps set of step name strings).
    """
    with open(run_metadata_file) as f:
        ext = os.path.splitext(run_metadata_file)[1].lower()
        if ext in ('.yaml', '.yml'):
            run_meta = yaml.safe_load(f)
        else:
            run_meta = json.load(f)
    completed_steps = {s['step'] for s in run_meta.get('steps', []) if s['status'] == 'success'}
    return run_meta, completed_steps


def record_step_complete(
    step_name: str,
    status: str,
    run_meta: dict,
    run_metadata_file: str,
    details: Optional[dict] = None,
) -> None:
    """Append a timestamped step record to the run metadata and flush to disk.

    Call this immediately after each major pipeline step completes so that
    the YAML file always reflects the latest progress, even if the run is
    interrupted later.

    Args:
        step_name:         Human-readable name for the completed step
                           (e.g. 'data_processing', 'subset_Control').
        status:            Outcome of the step, typically 'success' or 'error'.
        run_meta:          The in-memory metadata dict returned by
                           `log_run_metadata`.
        run_metadata_file: Path to the YAML file to update.
        details:           Optional dict of extra metadata to include in the step entry.
    """
    entry = {'step': step_name, 'status': status, 'timestamp': _utcnow()}
    if details:
        entry['details'] = details
    run_meta['steps'].append(entry)
    _write_yaml(run_meta, run_metadata_file)


def finalise_run_metadata(
    run_metadata_file: str,
    run_meta: dict,
    exit_status: str = 'success',
    error: Optional[str] = None,
) -> None:
    """Write final fields to the run metadata file at pipeline end.

    Adds end timestamp, duration, exit status, optional error message, and
    resource usage (peak RSS memory and total CPU time) if psutil is available.

    Args:
        run_metadata_file: Path to the YAML metadata file to update.
        run_meta:          The in-memory metadata dict returned by
                           `log_run_metadata`.
        exit_status:       'success' or 'error'.
        error:             Exception message string, included only when
                           exit_status is 'error'.
    """
    end_time = _utcnow()
    run_block = run_meta.setdefault('run', {})
    run_block['end_time'] = end_time
    run_block['exit_status'] = exit_status

    try:
        start = datetime.fromisoformat(run_block['start_time'])
        end = datetime.fromisoformat(end_time)
        run_block['duration_s'] = round((end - start).total_seconds(), 3)
    except Exception:
        run_block['duration_s'] = None

    if error is not None:
        run_block['error'] = error

    if _PSUTIL_AVAILABLE:
        try:
            proc = psutil.Process()
            run_block['resources'] = {
                'peak_memory_mb': proc.memory_info().rss / 1e6,
                'cpu_time_s': sum(proc.cpu_times()[:2]),
            }
        except Exception:
            pass

    _write_yaml(run_meta, run_metadata_file)
