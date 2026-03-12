"""Unit tests for run_metadata.py.

Structure mirrors the module: one section per public function, in the order
they are typically called during a pipeline run.

Subprocess calls (git, conda) are always mocked so tests run without a git
repo or conda installation. File I/O uses pytest's `tmp_path` fixture, which
provides a fresh temporary directory per test that is cleaned up automatically.

Mocking approach
----------------
@patch replaces the named object for the duration of a single test. When
stacking multiple @patch decorators, the innermost decorator maps to the
first extra argument, i.e.:

    @patch("module.outer")   ->  second extra arg (mock_outer)
    @patch("module.inner")   ->  first extra arg  (mock_inner)
    def test_foo(mock_inner, mock_outer): ...

side_effect on a mock returns different values on successive calls, used to
simulate multiple subprocess.check_output calls (e.g. git rev-parse, then
git diff, then git describe).
"""

import json
import os
from unittest.mock import MagicMock, patch

import pytest

import autoprot.utils.metadata_recording as rm


# ---------------------------------------------------------------------------
# file_hash
# ---------------------------------------------------------------------------

def test_file_hash_known_value(tmp_path):
    """file_hash produces the correct SHA-256 digest for a known input.

    We independently compute the expected hash with hashlib and compare,
    confirming the function is not just returning a constant or truncating.
    """
    import hashlib
    content = b"hello world"
    f = tmp_path / "test.txt"
    f.write_bytes(content)
    expected = hashlib.sha256(content).hexdigest()
    assert rm.file_hash(str(f)) == expected


def test_file_hash_empty_file(tmp_path):
    """file_hash handles an empty file correctly.

    The SHA-256 of an empty byte string is a well-known constant, so this
    also serves as a cross-check that the function reads the file rather than
    hashing the path string or some other artefact.
    """
    f = tmp_path / "empty.txt"
    f.write_bytes(b"")
    assert rm.file_hash(str(f)) == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"


def test_file_hash_different_files_differ(tmp_path):
    """file_hash returns distinct digests for files with different content.

    Guards against a broken implementation that returns the same value
    regardless of content.
    """
    f1 = tmp_path / "a.txt"
    f2 = tmp_path / "b.txt"
    f1.write_bytes(b"aaa")
    f2.write_bytes(b"bbb")
    assert rm.file_hash(str(f1)) != rm.file_hash(str(f2))


# ---------------------------------------------------------------------------
# _utcnow
# ---------------------------------------------------------------------------

def test_utcnow_is_iso_string():
    """_utcnow returns a timezone-aware ISO-8601 timestamp string.

    We parse the result with datetime.fromisoformat and check that tzinfo is
    set, confirming the string is both valid ISO-8601 and explicitly UTC
    rather than a naive local timestamp.
    """
    from datetime import datetime, timezone
    result = rm._utcnow()
    parsed = datetime.fromisoformat(result)
    assert parsed.tzinfo is not None


# ---------------------------------------------------------------------------
# get_git_info
# ---------------------------------------------------------------------------

# A fake 36-character commit SHA used throughout git tests.
GOOD_COMMIT = "abc123def456" * 3


@patch("autoprot.utils.metadata_recording.subprocess.check_output")
@patch("autoprot.utils.metadata_recording.subprocess.call")
def test_get_git_info_clean_repo(mock_call, mock_output):
    """get_git_info returns commit, tag, dirty=False when working tree is clean.

    subprocess.call (used for git diff-index) returns 0 to signal a clean
    tree. check_output is called twice: once for rev-parse HEAD and once for
    git describe. We verify all four returned keys are populated correctly and
    diff_file is None (no patch should be written for a clean repo).
    """
    mock_output.side_effect = [
        (GOOD_COMMIT + "\n").encode(),   # first call:  git rev-parse HEAD
        b"v1.2.3\n",                     # second call: git describe --tags
    ]
    mock_call.return_value = 0           # 0 = clean working tree

    result = rm.get_git_info()

    assert result["commit"] == GOOD_COMMIT
    assert result["dirty"] is False
    assert result["diff_file"] is None
    assert result["tag"] == "v1.2.3"


@patch("autoprot.utils.metadata_recording.subprocess.check_output")
@patch("autoprot.utils.metadata_recording.subprocess.call")
def test_get_git_info_dirty_writes_patch(mock_call, mock_output, tmp_path):
    """get_git_info writes a .patch file and records its path when repo is dirty.

    subprocess.call returns 1 (non-zero = dirty). check_output is called three
    times: rev-parse HEAD, git diff (output written to the patch file), and
    git describe. We assert the patch file exists and contains exactly the
    mocked diff content, confirming the write path works end-to-end.
    """
    diff_content = "diff --git a/foo.py b/foo.py\n+new line\n"
    mock_output.side_effect = [
        (GOOD_COMMIT + "\n").encode(),   # git rev-parse HEAD
        diff_content.encode(),           # git diff -> written to patch file
        b"v1.2.3-1-gabcdef\n",           # git describe
    ]
    mock_call.return_value = 1           # 1 = dirty working tree

    patch_file = str(tmp_path / "run.patch")
    result = rm.get_git_info(diff_file=patch_file)

    assert result["dirty"] is True
    assert result["diff_file"] == patch_file
    assert os.path.exists(patch_file)
    assert open(patch_file).read() == diff_content


@patch("autoprot.utils.metadata_recording.subprocess.check_output", side_effect=Exception("no git"))
@patch("autoprot.utils.metadata_recording.subprocess.call", side_effect=Exception("no git"))
def test_get_git_info_no_git(mock_call, mock_output):
    """get_git_info returns all-None dict gracefully when git is unavailable.

    Covers the case where a user has downloaded the code as a zip (no .git
    directory) or git is simply not installed. The function must not raise;
    it should return a safe sentinel dict so log_run_metadata can still write
    valid JSON.
    """
    result = rm.get_git_info()
    assert result == {"commit": None, "tag": None, "dirty": None, "diff_file": None}


# ---------------------------------------------------------------------------
# get_conda_env
# ---------------------------------------------------------------------------

@patch("autoprot.utils.metadata_recording.subprocess.check_output")
def test_get_conda_env_active(mock_output):
    """get_conda_env with no argument queries the currently active environment.

    Checks that 'conda list' is called without a -n flag, and that the raw
    text output is returned as a decoded string (not bytes).
    """
    mock_output.return_value = b"# packages in environment\nnumpy 1.26.0\n"
    result = rm.get_conda_env()
    assert "numpy" in result
    mock_output.assert_called_once_with(["conda", "list"])


@patch("autoprot.utils.metadata_recording.subprocess.check_output")
def test_get_conda_env_named(mock_output):
    """get_conda_env passes -n <env_name> when an environment name is given.

    The pipeline records both the active auto-proteomics env and the separate
    r-limma-env. This test confirms the -n flag is forwarded correctly so the
    right environment is queried.
    """
    mock_output.return_value = b"pandas 2.0.0\n"
    result = rm.get_conda_env("r-limma-env")
    mock_output.assert_called_once_with(["conda", "list", "-n", "r-limma-env"])
    assert result == "pandas 2.0.0\n"


@patch("autoprot.utils.metadata_recording.subprocess.check_output", side_effect=Exception("conda not found"))
def test_get_conda_env_failure(mock_output):
    """get_conda_env returns None rather than raising when conda is unavailable.

    Mirrors the git fallback behaviour: the caller can check for None and
    still write valid metadata rather than crashing the whole pipeline.
    """
    assert rm.get_conda_env() is None


# ---------------------------------------------------------------------------
# log_run_metadata
# ---------------------------------------------------------------------------

def _mock_git_info(**kwargs):
    """Return a minimal git info dict, overriding specific keys via kwargs.

    Used to keep log_run_metadata tests readable — we swap in a predictable
    git state without repeating the full dict in every test.
    """
    defaults = {"commit": GOOD_COMMIT, "tag": "v1.0", "dirty": False, "diff_file": None}
    defaults.update(kwargs)
    return defaults


@patch("autoprot.utils.metadata_recording.get_conda_env", return_value="numpy 1.26\n")
@patch("autoprot.utils.metadata_recording.get_git_info", return_value=_mock_git_info())
def test_log_run_metadata_creates_file(mock_git, mock_conda, tmp_path):
    """log_run_metadata writes a valid JSON file containing all expected keys.

    git and conda calls are mocked so the test does not require a real repo
    or conda installation. We verify that:
      - the file is created on disk and is valid JSON
      - input file hashes, config, git info, and platform are all present
      - steps is initialised as an empty list (no steps have run yet)
      - start_time is recorded
    """
    input_file = tmp_path / "input.txt"
    input_file.write_bytes(b"data")
    out_path = str(tmp_path)
    meta_file = str(tmp_path / "run_metadata.json")

    result_file, result_meta = rm.log_run_metadata(
        input_files=[str(input_file)],
        args=["main.py", "--config", "cfg.yaml"],
        config={"threshold": 0.05},
        out_path=out_path,
        run_metadata_file=meta_file,
    )

    assert os.path.exists(result_file)
    with open(result_file) as f:
        on_disk = json.load(f)

    assert on_disk["git"]["commit"] == GOOD_COMMIT
    assert on_disk["config"] == {"threshold": 0.05}
    assert str(input_file) in on_disk["input_files"]
    assert on_disk["steps"] == []
    assert "start_time" in on_disk
    assert "platform" in on_disk


# ---------------------------------------------------------------------------
# load_run_metadata
# ---------------------------------------------------------------------------

def _write_meta(path, steps):
    """Write a minimal run_metadata JSON fixture to disk.

    Helper used by load_run_metadata tests to avoid repeating boilerplate.
    Only includes the fields that load_run_metadata actually reads.
    """
    meta = {"start_time": "2024-01-01T00:00:00+00:00", "steps": steps}
    with open(path, "w") as f:
        json.dump(meta, f)
    return meta


def test_load_run_metadata_completed_steps(tmp_path):
    """load_run_metadata returns only 'success' steps in the completed set.

    This set is used during pipeline resume to decide which steps to skip.
    We include two successful steps and one failed step, and assert that the
    failed step is excluded — re-running a failed step on resume is
    intentional behaviour.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    _write_meta(meta_file, [
        {"step": "data_processing", "status": "success", "timestamp": "t1"},
        {"step": "normalisation",   "status": "success", "timestamp": "t2"},
        {"step": "subset_ctrl",     "status": "error",   "timestamp": "t3"},
    ])
    meta, completed = rm.load_run_metadata(meta_file)
    assert completed == {"data_processing", "normalisation"}
    assert "subset_ctrl" not in completed


def test_load_run_metadata_no_steps(tmp_path):
    """load_run_metadata returns an empty set when no steps have been recorded.

    Covers a fresh run interrupted before any step completed. The function
    should not crash and should signal that all steps need to be (re-)run.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    _write_meta(meta_file, [])
    _, completed = rm.load_run_metadata(meta_file)
    assert completed == set()


# ---------------------------------------------------------------------------
# record_step_complete
# ---------------------------------------------------------------------------

def test_record_step_complete_appends_and_flushes(tmp_path):
    """record_step_complete appends each step and keeps the JSON file in sync.

    We call it twice and verify that both the in-memory dict and the on-disk
    JSON reflect both entries. The flush-to-disk behaviour is critical: if
    the run is interrupted after this call, the completed step is not lost.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    run_meta = {"steps": []}

    rm.record_step_complete("data_processing", "success", run_meta, meta_file)
    rm.record_step_complete("normalisation", "success", run_meta, meta_file)

    assert len(run_meta["steps"]) == 2
    with open(meta_file) as f:
        on_disk = json.load(f)
    assert len(on_disk["steps"]) == 2
    assert on_disk["steps"][0]["step"] == "data_processing"


def test_record_step_complete_with_details(tmp_path):
    """record_step_complete stores an optional details dict inside the step entry.

    The details argument lets callers attach step-specific metadata (e.g.
    number of proteins quantified). We confirm the dict is preserved verbatim.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    run_meta = {"steps": []}

    rm.record_step_complete(
        "subset_ctrl", "success", run_meta, meta_file,
        details={"n_proteins": 1200}
    )

    assert run_meta["steps"][0]["details"] == {"n_proteins": 1200}


def test_record_step_complete_error_status(tmp_path):
    """record_step_complete stores 'error' status faithfully.

    Failed steps must be recorded (not silently dropped) so that resume
    logic can identify which steps need to be re-run.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    run_meta = {"steps": []}

    rm.record_step_complete("normalisation", "error", run_meta, meta_file)

    assert run_meta["steps"][0]["status"] == "error"


# ---------------------------------------------------------------------------
# finalise_run_metadata
# ---------------------------------------------------------------------------

def test_finalise_adds_end_time_and_status(tmp_path):
    """finalise_run_metadata adds end_time and exit_status and flushes to disk.

    Both the in-memory dict and the JSON file on disk are checked, confirming
    the final state is durable even if the process exits immediately after.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    run_meta = {"steps": []}

    rm.finalise_run_metadata(meta_file, run_meta, exit_status="success")

    assert "end_time" in run_meta
    assert run_meta["exit_status"] == "success"
    with open(meta_file) as f:
        on_disk = json.load(f)
    assert on_disk["exit_status"] == "success"


def test_finalise_records_error_message(tmp_path):
    """finalise_run_metadata stores the exception message when exit_status is 'error'.

    The error string lets someone reading the metadata understand what went
    wrong without having to dig through log files.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    run_meta = {"steps": []}

    rm.finalise_run_metadata(meta_file, run_meta, exit_status="error", error="ValueError: bad input")

    assert run_meta["error"] == "ValueError: bad input"


def test_finalise_with_psutil(tmp_path):
    """finalise_run_metadata records memory and CPU usage when psutil is available.

    psutil.Process is mocked to return fixed RSS memory (500 MB) and CPU times
    (10 s user + 2 s system = 12 s total). We use pytest.approx because the
    MB conversion involves floating-point division. _PSUTIL_AVAILABLE is patched
    to True so the resource block is entered regardless of whether psutil is
    actually installed in the test environment.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    run_meta = {"steps": []}

    mock_proc = MagicMock()
    mock_proc.memory_info.return_value.rss = 500_000_000   # bytes -> 500 MB
    mock_proc.cpu_times.return_value = (10.0, 2.0, 0.0, 0.0)  # user, system, ...

    with patch("autoprot.utils.metadata_recording._PSUTIL_AVAILABLE", True), \
         patch("autoprot.utils.metadata_recording.psutil.Process", return_value=mock_proc):
        rm.finalise_run_metadata(meta_file, run_meta)

    assert run_meta["resources"]["peak_memory_mb"] == pytest.approx(500.0)
    assert run_meta["resources"]["cpu_time_s"] == pytest.approx(12.0)


def test_finalise_without_psutil(tmp_path):
    """finalise_run_metadata omits the resources key when psutil is not installed.

    psutil is an optional dependency. If absent, finalise should complete
    normally and simply not record resource usage — downstream code reading
    the metadata must handle the missing key gracefully.
    """
    meta_file = str(tmp_path / "run_metadata.json")
    run_meta = {"steps": []}

    with patch("autoprot.utils.metadata_recording._PSUTIL_AVAILABLE", False):
        rm.finalise_run_metadata(meta_file, run_meta)

    assert "resources" not in run_meta