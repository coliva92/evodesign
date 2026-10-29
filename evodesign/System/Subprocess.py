import subprocess
from typing import Generator, List


def stream_command_line(command: List[str]) -> Generator[str, None, None]:
    """Runs a command and yields stdout lines one-by-one in real time."""
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    if process.stdout is not None:
        for line in process.stdout:
            yield line
    return_code = process.wait()
    if return_code != 0:
        err = process.stderr.read() if process.stderr else ""
        raise subprocess.CalledProcessError(return_code, command, stderr=err)


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
