from typing import List
import subprocess


def _open_command_line(command: List[str]):
    popen = subprocess.Popen(command, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, command)


def run_subprocess(command: List[str]):
    for line in _open_command_line(command):
        print(line, end="")
