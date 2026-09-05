"""Warnings this run raises, and what was done about each one.

Two warnings appear on every MD stage. Both were chased down before being
silenced, because a suppressed warning that mattered is the worst outcome
available here, and both are recorded with the check that cleared them so
nobody has to chase them twice.

**PresetChargesAndVirtualSitesWarning**, from openff-interchange:

    Preset charges were provided (via `charge_from_molecules`) alongside a force
    field that includes virtual site parameters.

This one deserved the look. If it were true, the ligand's electrostatics in the
production run would not be the charges the toolkit computed, which would be a
real and invisible error. It is not true here. The warning fires on the mere
presence of a `VirtualSites` handler in the force field, and openmmforcefields
registers that handler **empty**: openff-2.1.0 defines no virtual site
parameters at all, so none can match and none are created. Measured on the EGFR
run, with erlotinib:

* the built system has 52 particles for 52 atoms, so no virtual site exists
* every ligand charge is bit-identical to a clean AM1BCC assignment made with no
  preset charges anywhere (max difference 0.00e+00 e), and the total is 0.000000 e

So the charges are exactly what they would be without preset charges, and the
warning is about a thing that did not happen. Re-run those two checks before
changing `ligand_forcefield` to a force field that does define virtual sites,
because then it would stop being a false alarm.

**FutureWarning: `torch.distributed.reduce_op` is deprecated**:

Raised from openff-interchange's cache-clearing walk, which touches every
attribute of every loaded module looking for `functools.lru_cache` wrappers and
trips torch's deprecation shim on the way past. Nothing in this bundle uses
`torch.distributed`; ESMFold runs single-process and OpenMM does not use torch
at all. Cosmetic, and not ours to fix.

Everything else stays visible. Anything not on the list below is reformatted
rather than hidden: one line, in the run's own colours, so a warning that
matters reads like the rest of the output instead of three lines of Python
traceback furniture in the middle of a progress bar.
"""

from __future__ import annotations

import contextlib
import os
import sys
import warnings
from typing import Optional

# (message fragment, category name). Matched on the fragment because the
# categories live in modules that are expensive to import and may not be
# installed in every environment this file is read from.
SILENCED = [
    ("Preset charges were provided", "PresetChargesAndVirtualSitesWarning"),
    ("torch.distributed.reduce_op", "FutureWarning"),
]


def install(console=None) -> None:
    """Silence the two cleared warnings; route the rest through the console."""
    for fragment, _category in SILENCED:
        # re.escape is not used: these fragments are literal enough, and the
        # filter matches the START of the message, which is what `message=`
        # means here.
        warnings.filterwarnings("ignore", message=r".*" + fragment + r".*")

    if console is None:
        return

    def show(message, category, filename, lineno, file=None, line=None):
        text = str(message).split("\n")[0].strip()
        where = filename.split("site-packages/")[-1]
        console.warn(f"{category.__name__}: {text}")
        console.detail(f"raised at {where}:{lineno}")

    warnings.showwarning = show


@contextlib.contextmanager
def hush_c_stdout():
    """Swallow what C libraries print straight to file descriptor 1.

    MDTraj's DCD plugin writes two lines every time a trajectory is opened:

        dcdplugin) detected standard 32-bit DCD file of native endianness
        dcdplugin) CHARMM format DCD file (also NAMD 2.1 and later)

    It comes from C, so no warnings filter and no `sys.stdout` swap can reach
    it, and it lands in the middle of whatever progress bar is drawing. The file
    descriptor is the only lever there is. Anything else written while this is
    held is not thrown away: it is returned to the caller through the list it
    yields, so a library that says something worth hearing is not lost.
    """
    kept: list[str] = []
    sys.stdout.flush()
    saved = os.dup(1)
    read_fd, write_fd = os.pipe()
    os.dup2(write_fd, 1)
    os.close(write_fd)
    try:
        yield kept
    finally:
        # sys.stdout.flush() empties Python's buffer, not the C library's, and
        # the dcdplugin lines are written by C. Left unflushed they sit in libc
        # until the process exits, by which time the real descriptor is back and
        # they print after the run has finished: which is exactly what happened,
        # the two lines arriving under the summary table. fflush(NULL) flushes
        # every open C stream.
        sys.stdout.flush()
        try:
            import ctypes

            ctypes.CDLL(None).fflush(None)
        except (OSError, AttributeError):        # no libc handle: nothing to do
            pass
        os.dup2(saved, 1)
        os.close(saved)
        os.set_blocking(read_fd, False)
        try:
            captured = os.read(read_fd, 1 << 20).decode("utf-8", "replace")
        except BlockingIOError:
            captured = ""
        os.close(read_fd)
        kept.extend(line for line in captured.splitlines()
                    if line.strip() and not line.startswith("dcdplugin)"))
