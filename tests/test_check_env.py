# Allow pytest to import your project’s code
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import os
import subprocess

# Import the function we want to test
from src.utils.check_env import get_repo_root


def test_get_repo_root_when_git_succeeds(tmp_path, monkeypatch):
    """
    If `git rev-parse --show-toplevel` works, get_repo_root()
    should return exactly that path.
    """
    # 1. Prepare a fake repository root path
    fake_root = tmp_path / "myrepo"
    # 2. Create a CompletedProcess object as subprocess.run would
    completed = subprocess.CompletedProcess(
        args=["git", "rev-parse"],  # the command that was “run”
        returncode=0,  # indicates success
        stdout=str(fake_root) + "\n",  # stdout is the fake path + newline
        stderr="",
    )
    # 3. Monkeypatch subprocess.run so it always returns our fake result
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: completed)

    # 4. Call the function under test
    result = get_repo_root()

    # 5. Assert it returns exactly the path we faked
    assert result == str(fake_root)


def test_get_repo_root_when_not_git_repo(tmp_path, monkeypatch):
    """
    If `git rev-parse` fails (raises CalledProcessError),
    get_repo_root() should fall back to os.getcwd().
    """

    # 1. Make subprocess.run raise the git-error
    def fake_run(*args, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=args[0])

    monkeypatch.setattr(subprocess, "run", fake_run)

    # 2. Fake the current working directory
    fake_cwd = tmp_path / "current_dir"
    monkeypatch.setattr(os, "getcwd", lambda: str(fake_cwd))

    # 3. Run the function
    result = get_repo_root()

    # 4. It should now be the absolute path to our fake cwd
    assert result == os.path.abspath(str(fake_cwd))
