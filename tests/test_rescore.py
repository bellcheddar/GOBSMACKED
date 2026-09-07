"""The second opinion on pose ranking.

It exists because the search was reliably better than the ranking: on five
starting structures a near-native pose was in the list every time, and docking a
ligand back into its own crystal ranked the 1.51 A pose eighth of ten. What it
must never do is act on that opinion, because Vinardo was better than the
shipped ranking on four of those five receptors and worse on the one that
actually worked.
"""

from __future__ import annotations

import subprocess
import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bundle_template"))

from gobsmacked_run import dock, rescore  # noqa: E402

SMINA_OUT = """Affinity: -6.311 (kcal/mol)
Affinity: -6.694 (kcal/mol)
Affinity: -7.102 (kcal/mol)
"""


def _fake_run(stdout="", returncode=0):
    def run(cmd, **kwargs):
        return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr="")
    return run


def test_it_reports_a_disagreement_without_acting_on_it(tmp_path, monkeypatch):
    receptor, poses = tmp_path / "r.pdb", tmp_path / "p.sdf"
    receptor.write_text("ATOM\n"); poses.write_text("x\n")
    monkeypatch.setattr(rescore.subprocess, "run", _fake_run(SMINA_OUT))
    monkeypatch.setattr(rescore.shutil, "which", lambda name: "/usr/bin/pixi")

    block = rescore.run(receptor, poses, tmp_path / "rescore.csv", lambda m: None)
    assert block["ran"] is True
    assert block["top_pose"] == 3          # -7.102 is the best of the three
    assert block["agrees"] is False
    assert block["ranking"] == [3, 2, 1]
    # The written file keeps the engine's own pose numbering, so the two score
    # columns can be read side by side.
    rows = (tmp_path / "rescore.csv").read_text(encoding="utf-8").strip().splitlines()
    assert rows[0].startswith("pose_id,vinardo")
    assert rows[1].startswith("pose1,-6.311,3")


def test_agreement_is_reported_too(tmp_path, monkeypatch):
    receptor, poses = tmp_path / "r.pdb", tmp_path / "p.sdf"
    receptor.write_text("ATOM\n"); poses.write_text("x\n")
    monkeypatch.setattr(rescore.subprocess, "run",
                        _fake_run("Affinity: -9.0 (kcal/mol)\nAffinity: -1.0 (kcal/mol)\n"))
    monkeypatch.setattr(rescore.shutil, "which", lambda name: "/usr/bin/pixi")
    block = rescore.run(receptor, poses, tmp_path / "r.csv", lambda m: None)
    assert block["agrees"] is True and block["top_pose"] == 1


def test_every_failure_is_a_missing_column_not_a_dead_run(tmp_path, monkeypatch):
    """Docking has already cost minutes by this point. A decoration that could
    end the run would be a bad trade, so each way of failing is checked."""
    receptor, poses = tmp_path / "r.pdb", tmp_path / "p.sdf"
    receptor.write_text("ATOM\n"); poses.write_text("x\n")
    log = lambda m: None

    # No inputs.
    assert rescore.run(tmp_path / "nope.pdb", poses, tmp_path / "a.csv", log)["ran"] is False

    monkeypatch.setattr(rescore.shutil, "which", lambda name: "/usr/bin/pixi")
    # Non-zero exit.
    monkeypatch.setattr(rescore.subprocess, "run", _fake_run("boom", returncode=1))
    assert rescore.run(receptor, poses, tmp_path / "b.csv", log)["ran"] is False
    # Zero exit but nothing parseable.
    monkeypatch.setattr(rescore.subprocess, "run", _fake_run("hello"))
    assert rescore.run(receptor, poses, tmp_path / "c.csv", log)["ran"] is False
    # The subprocess itself blowing up.
    def explode(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd, 1)
    monkeypatch.setattr(rescore.subprocess, "run", explode)
    assert rescore.run(receptor, poses, tmp_path / "d.csv", log)["ran"] is False
    # And with no pixi at all.
    monkeypatch.setattr(rescore.shutil, "which", lambda name: None)
    assert rescore.run(receptor, poses, tmp_path / "e.csv", log)["ran"] is False


def test_dock_imports_rescore_at_module_level():
    """It was briefly inserted inside a function, which dedented the line after
    it and left dock.py unparseable. The stage that calls it must be importable."""
    import ast
    tree = ast.parse(Path(dock.__file__).read_text(encoding="utf-8"))
    top = {ast.unparse(n) for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))}
    assert "from . import rescore" in top
