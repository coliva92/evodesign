import os
import subprocess
from typing import Generator, List


def stream_command_line(command: List[str]) -> Generator[str, None, None]:
    """Runs a command and yields combined stdout and stderr lines in real time."""
    # Force child Python/Docker processes to flush output immediately
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # Interleave stderr into stdout to prevent pipe deadlock
        text=True,
        bufsize=1,
        env=env,
    )

    if process.stdout is not None:
        for line in process.stdout:
            yield line

    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)


def run_subprocess(command: List[str]) -> None:
    """
    Executes a command and streams its standard output to the console in real-time.

    Args:
        command: The command to execute as a list of strings.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
    """
    for line in stream_command_line(command):
        print(line, end="")
    return
