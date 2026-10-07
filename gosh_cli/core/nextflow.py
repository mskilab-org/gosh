import os
import select
import subprocess
import sys
import time
from datetime import datetime


def _process_is_finished(pid):
    try:
        with open(f"/proc/{pid}/stat", "r", encoding="utf-8") as handle:
            return handle.read().split()[2] == "Z"
    except FileNotFoundError:
        return True


def _wait_for_process(pid):
    if hasattr(os, "pidfd_open"):
        try:
            pid_fd = os.pidfd_open(pid)
        except ProcessLookupError:
            return
        try:
            poller = select.poll()
            poller.register(pid_fd, select.POLLIN)
            while not poller.poll(1000):
                if _process_is_finished(pid):
                    return
        finally:
            os.close(pid_fd)
        return

    while not _process_is_finished(pid):
        time.sleep(5)


def _finalize_background_run(pid_file, gosh_command):
    with open(pid_file, "r", encoding="utf-8") as handle:
        pid = int(handle.read().strip())

    try:
        _wait_for_process(pid)
        subprocess.run([gosh_command, "run", "post"], check=True)
        subprocess.run([gosh_command, "protect"], check=True)
    finally:
        try:
            os.remove(pid_file)
        except FileNotFoundError:
            pass


def _background_finalizer_main():
    _finalize_background_run(sys.argv[1], sys.argv[2])
class NextflowRunner:
    def __init__(self):
        self.cmd = 'nextflow'

    def get_timestamp(self):
        """Get current timestamp in YYYYMMDD_HHMMSS format"""
        return datetime.now().strftime('%Y%m%d_%H%M%S')

    def run(self, command):
        """Run nextflow command with given command string"""
        try:
            subprocess.run(command, shell=True, check=True, executable = "/bin/bash")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Pipeline execution failed: {e}")

    def finalize_after_background_run(self, pid_file, execution_log):
        """Start a detached process that finalizes a background Nextflow run."""
        log_handle = open(execution_log, "a", encoding="utf-8")
        try:
            subprocess.Popen(
                [sys.executable, __file__, pid_file, "gosh"],
                stdin=subprocess.DEVNULL,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        finally:
            log_handle.close()


if __name__ == "__main__":
    _background_finalizer_main()
