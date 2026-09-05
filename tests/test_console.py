"""The bundle's terminal output.

Worth testing because it degrades in four directions at once (no colour, no
unicode, no terminal to redraw, no total to count towards) and every one of
those paths runs on somebody else's machine where nobody will see it fail.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bundle_template"))

from gobsmacked_run import console as con  # noqa: E402


def plain(stream: io.StringIO) -> str:
    return con.ANSI.sub("", stream.getvalue())


@pytest.fixture()
def piped(monkeypatch):
    """A console writing to a pipe: no colour, no redraw, ASCII only if needed."""
    monkeypatch.delenv("FORCE_COLOR", raising=False)
    monkeypatch.setenv("NO_COLOR", "1")
    stream = io.StringIO()
    return con.Console(stream=stream), stream


def test_a_pipe_gets_no_escape_sequences(piped):
    console, stream = piped
    console.banner("gs_test", "a title")
    console.stage_start(1, 5, "dock", "docking", 300)
    console.stage_end("dock", 12.0, "10 poses")
    assert "\x1b[" not in stream.getvalue()


def test_the_log_file_never_gets_escape_sequences(tmp_path, monkeypatch):
    monkeypatch.setenv("FORCE_COLOR", "1")
    stream = io.StringIO()
    console = con.Console(log_path=tmp_path / "run.log", stream=stream)
    console.warn("something to say")
    assert "\x1b[" in stream.getvalue()          # the terminal did get colour
    assert "\x1b[" not in (tmp_path / "run.log").read_text()


@pytest.mark.parametrize("prefix", ["fold", "prep", "dock", "md", "summarise"])
def test_a_stage_prefix_is_stripped_from_its_own_lines(piped, prefix):
    console, stream = piped
    console.stage = prefix
    console(f"{prefix}: something happened")
    assert plain(stream).strip() == "something happened"


def test_a_prefix_is_stripped_even_after_the_stage_ended(piped):
    """pack() logs after the stage that called it has already finished."""
    console, stream = piped
    console.stage = None
    console("summarise: wrote results.tar.gz")
    assert plain(stream).strip() == "wrote results.tar.gz"


def test_a_bar_never_runs_past_the_terminal_width(monkeypatch):
    """A bar wider than the terminal wraps, and then \\r clears only the second
    row: the first stays on screen and every redraw leaves another copy."""
    monkeypatch.setenv("FORCE_COLOR", "1")
    stream = io.StringIO()
    console = con.Console(stream=stream)
    console.live = True
    for width in (52, 60, 88, 100):
        monkeypatch.setattr(type(console), "width", property(lambda self, w=width: w))
        with console.bar("a label of some length", total=1000) as bar:
            bar.update(371, note="24.7 ns/day, 400 of 1000 ps")
            rendered = con.ANSI.sub("", bar._render(bar.started + 30))
            assert len(rendered) <= width, f"{len(rendered)} > {width}"


def test_a_finished_bar_leaves_one_line_and_no_bar(monkeypatch):
    monkeypatch.setenv("NO_COLOR", "1")
    stream = io.StringIO()
    console = con.Console(stream=stream)
    console.live = True
    with console.bar("measuring", total=3) as bar:
        bar.update(3, note="done")
    lines = [line for line in plain(stream).split("\n") if line.strip()]
    assert len(lines) == 1
    assert "measuring" in lines[0]


def test_an_estimate_bar_admits_when_it_is_wrong(monkeypatch):
    monkeypatch.setenv("NO_COLOR", "1")
    console = con.Console(stream=io.StringIO())
    console.live = True
    bar = con.Bar(console, "docking", None, estimate_s=100.0, unit="")
    assert "left (estimate)" in con.ANSI.sub("", bar._render(bar.started + 30))
    assert "longer than" in con.ANSI.sub("", bar._render(bar.started + 300))


def test_a_bar_with_no_total_and_no_estimate_claims_nothing(monkeypatch):
    monkeypatch.setenv("NO_COLOR", "1")
    console = con.Console(stream=io.StringIO())
    console.live = True
    bar = con.Bar(console, "importing", None, None, "")
    rendered = con.ANSI.sub("", bar._render(bar.started + 5))
    assert "%" not in rendered and "left" not in rendered


class AsciiStream(io.StringIO):
    """StringIO's `encoding` is read-only, so the ASCII case needs its own."""
    encoding = "ascii"

    def isatty(self):
        return False


def test_ascii_only_terminals_get_ascii():
    stream = AsciiStream()
    console = con.Console(stream=stream)
    assert console.icons is con.ASCII_ICONS
    console.stage_end("dock", 3.0)
    stream.getvalue().encode("ascii")            # raises if a box character got through


def test_stages_without_a_console_still_run():
    """Unit tests pass a plain function for `log`; the bars must not care."""
    with con.bar_for(lambda message: None, "anything", total=10) as bar:
        bar.update(5, note="fine")
        bar.advance()


@pytest.mark.parametrize("seconds,expected", [
    (0, "0s"), (45, "45s"), (60, "1m 00s"), (392, "6m 32s"), (4211, "1h 10m"),
])
def test_durations_read_the_way_a_person_would_say_them(seconds, expected):
    assert con.human(seconds) == expected
