"""What the run looks like while it is running.

A GOBSMACKED run is five stages, most of a minute to most of an hour each, on
someone else's laptop. The only thing standing between a working run and a
killed one is whether the person watching can tell that it is still working, so
this file exists to answer three questions at a glance: which stage, how far
through, and how much longer.

Hand-rolled rather than rich or tqdm, because the bundle ships a locked pixi
environment: a progress bar is not worth a dependency in a lock file that has to
resolve on someone else's machine. It is a few hundred lines and it does exactly
what is needed.

Everything degrades. No colour when the output is a pipe, a file, a dumb
terminal or NO_COLOR is set; ASCII when the terminal cannot encode the box
characters; one plain line per update instead of a redrawn bar when there is no
terminal to redraw. The log file always gets the plain text, because a log full
of escape sequences is worse than no colour at all.
"""

from __future__ import annotations

import os
import re
import shutil
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Optional

ANSI = re.compile(r"\x1b\[[0-9;]*m")

# The site's palette, as close as 256-colour terminals get. Cyan is the run
# itself, amber anything optional or skipped, red a failure, green a finish:
# the same four meanings they carry on the results page.
COLOURS = {
    "phos":   "\x1b[38;5;80m",
    "amber":  "\x1b[38;5;215m",
    "red":    "\x1b[38;5;203m",
    "green":  "\x1b[38;5;114m",
    "grey":   "\x1b[38;5;246m",
    "dim":    "\x1b[38;5;240m",
    "bold":   "\x1b[1m",
    "off":    "\x1b[0m",
}

ICONS = {"done": "✓", "run": "▸", "skip": "·", "fail": "✗", "warn": "⚑",
         "clock": "⏱", "bar_full": "█", "bar_empty": "░",
         "tl": "╭", "tr": "╮", "bl": "╰", "br": "╯", "h": "─", "v": "│",
         "dot": "•"}
ASCII_ICONS = {"done": "+", "run": ">", "skip": ".", "fail": "x", "warn": "!",
               "clock": "~", "bar_full": "#", "bar_empty": "-",
               "tl": "+", "tr": "+", "bl": "+", "br": "+", "h": "-", "v": "|",
               "dot": "*"}

SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
ASCII_SPINNER = "|/-\\"


def human(seconds: float) -> str:
    """A duration someone can read at a glance: 45s, 6m 32s, 1h 14m."""
    seconds = max(0.0, float(seconds))
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        return f"{int(seconds // 60)}m {int(seconds % 60):02d}s"
    return f"{int(seconds // 3600)}h {int((seconds % 3600) // 60):02d}m"


class Console:
    """The run's voice. Callable, so `log("...")` from a stage still works."""

    def __init__(self, log_path: Optional[Path] = None, stream=None):
        self.stream = stream or sys.stdout
        self.log_path = log_path
        self.stage: Optional[str] = None
        self.colour = self._wants_colour()
        self.icons = ICONS if self._wants_unicode() else ASCII_ICONS
        self.spinner = SPINNER if self._wants_unicode() else ASCII_SPINNER
        self.live = bool(getattr(self.stream, "isatty", lambda: False)())
        self._open_bar: Optional["Bar"] = None

    # -- capability detection ------------------------------------------------
    def _wants_colour(self) -> bool:
        if os.environ.get("NO_COLOR") is not None:
            return False
        if os.environ.get("FORCE_COLOR"):
            return True
        if os.environ.get("TERM", "") in ("dumb", ""):
            return False
        return bool(getattr(self.stream, "isatty", lambda: False)())

    def _wants_unicode(self) -> bool:
        encoding = getattr(self.stream, "encoding", None) or ""
        try:
            "".join(ICONS.values()).encode(encoding or "ascii")
            return True
        except (UnicodeEncodeError, LookupError):
            return False

    @property
    def width(self) -> int:
        return max(52, min(100, shutil.get_terminal_size((88, 24)).columns))

    # -- the primitives ------------------------------------------------------
    def paint(self, text: str, *names: str) -> str:
        if not self.colour or not names:
            return text
        return "".join(COLOURS[n] for n in names) + text + COLOURS["off"]

    def _to_log(self, text: str) -> None:
        """Append one plain line to run.log, and never raise while doing it.

        Two ways this killed a run before it was written down. `open()` with no
        encoding uses the locale's, which inside a pixi task is ASCII, so the
        first bullet character in a bar's closing line raised UnicodeEncodeError
        seventeen minutes into an MD stage. The handler that reported the
        failure then wrote a cross through the same path and raised again, so
        the traceback the person saw was about the console rather than about
        their run.

        The encoding is now named rather than inherited, and any failure here is
        swallowed: this is a progress display, and there is no state it can be
        in that justifies discarding an hour of molecular dynamics.
        """
        if self.log_path is None:
            return
        try:
            stamp = datetime.now().strftime("%H:%M:%S")
            with open(self.log_path, "a", encoding="utf-8", errors="replace") as fh:
                fh.write(f"[{stamp}] {ANSI.sub('', text).rstrip()}\n")
        except Exception:                          # noqa: BLE001 - never fatal
            pass

    def _speakable(self, text: str) -> str:
        """Text the terminal can actually encode.

        The icons are chosen from the stream's encoding at startup, but a
        stage's own message is not: a residue name with an Angstrom sign in it
        would go to an ASCII terminal and raise from inside print().
        """
        encoding = getattr(self.stream, "encoding", None) or "utf-8"
        try:
            text.encode(encoding)
            return text
        except (UnicodeEncodeError, LookupError):
            return text.encode(encoding or "ascii", "replace").decode(encoding or "ascii")

    def write(self, line: str = "", to_log: bool = True) -> None:
        """One line to the terminal, the same line without escapes to the log."""
        self._clear_bar()
        print(self._speakable(line), file=self.stream, flush=True)
        if to_log:
            self._to_log(line)

    def __call__(self, message: str) -> None:
        """A stage's own log line.

        Stages prefix their messages with their own name, which was the only
        structure the output had before this file existed. Now that every line
        already sits under a stage heading, the prefix is redundant: it is
        stripped here rather than edited out of forty call sites, and kept in
        the log where there are no headings to sit under.
        """
        text = message
        for name in ("fold", "prep", "dock", "md", "summarise"):
            # Any stage's prefix, not only the running one: pack() logs after
            # the stage that called it has already ended.
            if text.lower().startswith(f"{name}: "):
                text = text[len(name) + 2:]
                break
        self.detail(text, logged=message)

    def detail(self, message: str, logged: Optional[str] = None) -> None:
        self._clear_bar()
        print(self._speakable("    " + self.paint(message, "grey")),
              file=self.stream, flush=True)
        self._to_log(logged or message)

    # -- the furniture -------------------------------------------------------
    def banner(self, job_id: str, title: str) -> None:
        i, w = self.icons, self.width
        inner = w - 4
        top = i["tl"] + i["h"] * inner + i["tr"]
        bottom = i["bl"] + i["h"] * inner + i["br"]
        name = "GOBSMACKED"
        label = f"  {name}   {i['dot']}   {job_id}"
        self.write("")
        self.write("  " + self.paint(top, "phos"))
        self.write("  " + self.paint(i["v"], "phos")
                   + self.paint(label.ljust(inner)[:inner], "phos", "bold")
                   + self.paint(i["v"], "phos"))
        self.write("  " + self.paint(bottom, "phos"))
        if title:
            self.write("  " + self.paint(title, "grey"))

    def plan(self, rows: list[tuple[str, str, str]], minutes: float) -> None:
        """The stage table, printed before anything expensive begins."""
        i = self.icons
        self.write("")
        for index, (name, state, note) in enumerate(rows, start=1):
            colour = {"run": "phos", "done": "green", "skip": "dim"}[state]
            icon = {"run": i["run"], "done": i["done"], "skip": i["skip"]}[state]
            self.write("   " + self.paint(icon, colour) + " "
                       + self.paint(f"{index}", "dim") + "  "
                       + self.paint(name.ljust(11), colour)
                       + self.paint(note, "grey"))
        if minutes > 0:
            finish = datetime.now() + timedelta(minutes=minutes)
            self.write("")
            self.write("     " + self.paint(
                f"{i['clock']} about {human(minutes * 60)} in total, "
                f"finishing around {finish.strftime('%H:%M')}", "amber"))
        self.write("")

    def stage_start(self, index: int, total: int, name: str, note: str,
                    estimate_s: Optional[float]) -> None:
        self.stage = name
        i, w = self.icons, self.width
        left = ("  " + self.paint(i["run"], "phos") + " "
                + self.paint(f"{index}/{total}", "dim") + "  "
                + self.paint(name.upper(), "phos", "bold"))
        right = self.paint(f"{i['clock']} ~{human(estimate_s)}", "dim") if estimate_s else ""
        pad = max(1, w - len(ANSI.sub("", left)) - len(ANSI.sub("", right)))
        self.write("")
        self.write(left + " " * pad + right)
        if note:
            self.write("     " + self.paint(note, "dim"))

    def stage_end(self, name: str, seconds: float, note: str = "") -> None:
        self.stage = None
        text = ("  " + self.paint(self.icons["done"], "green") + "  "
                + self.paint(f"{name} finished in {human(seconds)}", "green"))
        if note:
            text += self.paint(f"  {self.icons['dot']}  {note}", "grey")
        self.write(text)

    def warn(self, message: str) -> None:
        self.write("  " + self.paint(self.icons["warn"] + "  " + message, "amber"))

    def fail(self, message: str) -> None:
        self.write("  " + self.paint(self.icons["fail"] + "  " + message, "red", "bold"))

    def summary(self, timings: dict[str, float], warnings: list[str],
                archive: Optional[Path], job_id: str) -> None:
        i, w = self.icons, self.width
        self.write("")
        self.write("  " + self.paint(i["h"] * (w - 4), "dim"))
        total = sum(timings.values())
        for name, seconds in timings.items():
            share = seconds / total if total else 0.0
            cells = int(round(share * 24))
            bar = i["bar_full"] * cells + i["bar_empty"] * (24 - cells)
            self.write("   " + self.paint(name.ljust(11), "grey")
                       + self.paint(bar, "phos") + "  "
                       + self.paint(human(seconds).rjust(8), "grey"))
        self.write("   " + self.paint("total".ljust(11), "bold")
                   + " " * 24 + "  " + self.paint(human(total).rjust(8), "bold"))
        if warnings:
            self.write("")
            for warning in warnings:
                self.warn(warning)
        if archive is not None:
            self.write("")
            self.write("  " + self.paint(i["done"], "green") + "  "
                       + self.paint("Done.", "green", "bold") + " "
                       + self.paint("Upload", "grey") + " "
                       + self.paint(str(archive), "phos") + " "
                       + self.paint("on the Analyze tab.", "grey"))
        self.write("")

    # -- progress ------------------------------------------------------------
    def bar(self, label: str, total: Optional[float] = None,
            estimate_s: Optional[float] = None, unit: str = "") -> "Bar":
        return Bar(self, label, total, estimate_s, unit)

    def _clear_bar(self) -> None:
        if self._open_bar is not None:
            self._open_bar._erase()


class Bar:
    """A progress bar, a spinner, or a plain line, depending on what is known.

    Three cases, and the difference between them is honesty rather than looks:

    * `total` given: real progress, and the time left is measured from the rate
      the run is actually going at.
    * only `estimate_s` given: no real progress exists, so the bar tracks
      elapsed time against a prior guess and says so. Past the estimate it stops
      pretending and reports elapsed alone rather than a bar stuck at 99 %.
    * neither: a spinner with the elapsed time, which claims nothing at all.
    """

    MIN_REDRAW_S = 0.1

    def __init__(self, console: Console, label: str, total: Optional[float],
                 estimate_s: Optional[float], unit: str):
        self.console = console
        self.label = label
        self.total = total
        self.estimate_s = estimate_s
        self.unit = unit
        self.value = 0.0
        self.note = ""
        self.started = time.time()
        self.drawn = False
        self._last_draw = 0.0
        self._tick = 0
        self._last_plain = 0.0

    def __enter__(self) -> "Bar":
        self.console._open_bar = self
        self.draw(force=True)
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._erase()
        self.console._open_bar = None
        if exc_type is None:
            elapsed = time.time() - self.started
            self.console.detail(f"{self.label}: {human(elapsed)}"
                                + (f"  {self.console.icons['dot']}  {self.note}"
                                   if self.note else ""))

    def update(self, value: Optional[float] = None, note: str = "") -> None:
        if value is not None:
            self.value = value
        if note:
            self.note = note
        self.draw()

    def advance(self, delta: float = 1.0, note: str = "") -> None:
        self.update(self.value + delta, note)

    # -- drawing -------------------------------------------------------------
    def _erase(self) -> None:
        if self.drawn and self.console.live:
            print("\r\x1b[2K", end="", file=self.console.stream, flush=True)
            self.drawn = False

    def draw(self, force: bool = False) -> None:
        now = time.time()
        if not self.console.live:
            self._plain(now, force)
            return
        if not force and now - self._last_draw < self.MIN_REDRAW_S:
            return
        self._last_draw = now
        print("\r\x1b[2K" + self._render(now), end="", file=self.console.stream, flush=True)
        self.drawn = True

    def _plain(self, now: float, force: bool) -> None:
        """No terminal to redraw: one line a minute, so a log stays readable."""
        if not force and now - self._last_plain < 60:
            return
        self._last_plain = now
        self.console.detail(ANSI.sub("", self._render(now)).strip())

    def _render(self, now: float) -> str:
        """One line, never wider than the terminal.

        Width is the whole game here. A line one character too long wraps, and
        then the `\r` that starts the next redraw returns to the beginning of
        the SECOND row: `\x1b[2K` clears that row and the first one stays on
        screen. Every redraw leaves another copy, which is how a progress bar
        turns into the wall of repeated text it was added to prevent. So the
        pieces are measured as plain text, given budgets in priority order, and
        only coloured once they fit.
        """
        c, i = self.console, self.console.icons
        elapsed = now - self.started
        width = c.width

        label = self.label
        if len(label) > width // 2:
            label = label[:max(8, width // 2 - 1)] + ("…" if c.icons is ICONS else "~")
        head = "     " + label + " "

        if self.total:
            fraction = min(1.0, self.value / self.total) if self.total else 0.0
            rate = self.value / elapsed if elapsed > 0.5 and self.value else 0.0
            left = (self.total - self.value) / rate if rate > 0 else None
            tail_parts = [(f"{fraction * 100:4.0f}%", "phos"), ("  ", None),
                          (human(elapsed), "dim")]
            if left:
                tail_parts += [(f" {i['dot']} {human(left)} left", "amber")]
        elif self.estimate_s:
            fraction = min(1.0, elapsed / self.estimate_s)
            if elapsed > self.estimate_s:
                tail_parts = [(human(elapsed), "dim"),
                              (f" {i['dot']} longer than the "
                               f"{human(self.estimate_s)} estimate", "amber")]
            else:
                tail_parts = [(human(elapsed), "dim"),
                              (f" {i['dot']} ~{human(self.estimate_s - elapsed)} left "
                               f"(estimate)", "dim")]
        else:
            fraction = None
            tail_parts = [(c.spinner[self._tick % len(c.spinner)], "phos"),
                          (" ", None), (human(elapsed), "dim")]
            self._tick += 1

        tail = "".join(text for text, _ in tail_parts)

        # Budgets, in the order things may be dropped: the bar shrinks first,
        # then the note is cut, and the label and the numbers always survive.
        spare = width - len(head) - len(tail) - 2
        cells = 0 if fraction is None else max(0, min(30, spare - 2))
        if cells and cells < 8:                       # a stub bar reads as noise
            cells = 0
        note_room = width - len(head) - len(tail) - cells - (4 if cells else 2)
        note = ""
        if self.note and note_room > 6:
            note = self.note if len(self.note) <= note_room else self.note[:note_room - 1] + "…"
            note = "  " + note

        line = c.paint(head, "grey")
        if cells:
            full = int(round(fraction * cells))
            line += (c.paint(i["bar_full"] * full, "phos")
                     + c.paint(i["bar_empty"] * (cells - full), "dim") + "  ")
        line += "".join(c.paint(text, colour) if colour else text
                        for text, colour in tail_parts)
        line += c.paint(note, "dim")
        return line


class NullBar:
    """What a stage gets when it is driven by something that is not a Console:
    the unit tests, or a caller that passed a plain function for `log`. Every
    method exists and none of them draw anything."""

    def __enter__(self) -> "NullBar":
        return self

    def __exit__(self, *exc) -> None:
        return None

    def update(self, value: Optional[float] = None, note: str = "") -> None:
        return None

    def advance(self, delta: float = 1.0, note: str = "") -> None:
        return None


def bar_for(log: Any, label: str, total: Optional[float] = None,
            estimate_s: Optional[float] = None, unit: str = "") -> Any:
    """A progress bar if `log` can draw one, and a silent stand-in if it cannot."""
    maker = getattr(log, "bar", None)
    if maker is None:
        return NullBar()
    return maker(label, total=total, estimate_s=estimate_s, unit=unit)
