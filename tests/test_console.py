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


# --- the console must never be able to take a run down -----------------------

def test_the_log_names_its_encoding_rather_than_inheriting_it(tmp_path, monkeypatch):
    """An ASCII locale killed a real MD run seventeen minutes in.

    `open()` with no encoding uses the locale's, and inside a pixi task that is
    ASCII, so the bullet in a bar's closing line raised UnicodeEncodeError from
    the logging call and the handler reporting it raised again on the cross.

    Setting LC_ALL in a test proves nothing: Python reads the locale once at
    startup. So this asserts the thing that actually matters, that the call
    never leaves the encoding to be inherited, by making an inherited encoding
    fail the way the user's machine made it fail.
    """
    import builtins

    real_open = builtins.open

    def strict_open(file, mode="r", *args, **kwargs):
        if "b" not in mode and kwargs.get("encoding") is None:
            raise UnicodeEncodeError("ascii", "", 0, 1, "locale encoding inherited")
        return real_open(file, mode, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", strict_open)
    console = con.Console(log_path=tmp_path / "run.log", stream=io.StringIO())
    console.write("equilibrating 100 ps: 1m 34s  \u2022  restraint 200 kJ/mol/nm^2")
    console.fail("md failed: \u2717 something")
    monkeypatch.undo()
    written = (tmp_path / "run.log").read_text(encoding="utf-8")
    assert "\u2022" in written and "\u2717" in written


def test_no_bundle_file_leaves_its_text_encoding_to_the_locale():
    """The same trap, everywhere else it could be sprung.

    A campaign titled with a Greek beta would have failed on the read of
    campaign.yaml, before a single stage ran. Parsed rather than grepped: the
    first version of this test matched the sentence describing the bug in a
    docstring and failed on its own prose.
    """
    import ast

    bundle = Path(__file__).resolve().parents[1] / "bundle_template"
    offenders = []
    for path in sorted(list(bundle.glob("*.py")) + list(bundle.glob("*/*.py"))):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if isinstance(func, ast.Name) and func.id == "open":
                name = "open"
            elif isinstance(func, ast.Attribute) and func.attr in ("read_text", "write_text"):
                name = func.attr
            else:
                continue                    # tarfile.open and everything else
            mode = next((a.value for a in node.args
                         if isinstance(a, ast.Constant) and isinstance(a.value, str)), "")
            if name == "open" and "b" in mode:
                continue
            if any(k.arg == "encoding" for k in node.keywords):
                continue
            offenders.append(f"{path.name}:{node.lineno} {name}")
    assert not offenders, "text I/O with no encoding: " + ", ".join(offenders)


def test_no_subprocess_decodes_output_with_the_locales_encoding():
    """The same trap one level out, and it has now bitten twice.

    `text=True` decodes the child's output with the locale's encoding. Inside a
    pixi task that is ASCII, and boltz prints UTF-8, so a progress bar's box
    characters raised UnicodeDecodeError and killed the affinity stage after MD
    had already run. Scanned across the app as well as the bundle: PLIP is a
    subprocess there and prints whatever it likes.
    """
    import ast

    root = Path(__file__).resolve().parents[1]
    offenders = []
    for folder in ("bundle_template", "app"):
        for path in sorted((root / folder).rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
                if not isinstance(node, ast.Call):
                    continue
                func = node.func
                name = getattr(func, "attr", None) or getattr(func, "id", None)
                if name not in ("run", "Popen", "check_output"):
                    continue
                kwargs = {k.arg: k.value for k in node.keywords}
                text_mode = kwargs.get("text")
                if not (isinstance(text_mode, ast.Constant) and text_mode.value is True):
                    continue
                if "encoding" not in kwargs:
                    offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, ("subprocess in text mode with no encoding: "
                           + ", ".join(offenders))


def test_a_broken_log_never_reaches_the_stage(tmp_path):
    """There is no state a progress display can be in that justifies
    discarding an hour of molecular dynamics."""
    console = con.Console(log_path=tmp_path / "no-such-directory" / "run.log",
                          stream=io.StringIO())
    console.write("this must not raise")
    console.detail("nor this")
    console.fail("nor this")


def test_a_terminal_that_cannot_encode_the_message_still_gets_a_line():
    """The icons follow the stream's encoding, but a stage's own text does not:
    an Angstrom sign in a message would raise from inside print()."""
    stream = AsciiStream()
    console = con.Console(stream=stream)
    console.write("pocket volume 412 Å³")
    assert "pocket volume 412" in stream.getvalue()
