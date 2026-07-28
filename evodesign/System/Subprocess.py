import subprocess
from typing import List




def _open_command_line(command: List[str]):
    popen = subprocess.Popen(command, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, command)



def run_subprocess(command: List[str]) -> None:
    """
    Executes a command and streams its standard output to the console in real-time.
    
    Args:
        command: The command to execute as a list of strings.
    
    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
    """
    for line in _open_command_line(command):
        print(line, end="")
    return
