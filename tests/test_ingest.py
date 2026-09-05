"""A results archive is untrusted input, and is checked like one."""

from __future__ import annotations

import io
import json
import tarfile
from pathlib import Path

import pytest

from app.services import ingest


def repack(source: Path, dest: Path, drop: set[str] | None = None,
           replace: dict[str, bytes] | None = None) -> Path:
    """Copy an archive, optionally omitting or rewriting members."""
    drop = drop or set()
    replace = replace or {}
    with tarfile.open(source) as src, tarfile.open(dest, "w:gz") as out:
        for member in src.getmembers():
            rel = member.name.split("/", 1)[-1]
            if rel in drop:
                continue
            if rel in replace:
                data = replace[rel]
                info = tarfile.TarInfo(member.name)
                info.size = len(data)
                out.addfile(info, io.BytesIO(data))
                continue
            fh = src.extractfile(member)
            if fh is None:
                continue
            out.addfile(member, fh)
    return dest


def test_valid_archive_passes(egfr_archive, tmp_path):
    results = ingest.extract(egfr_archive, tmp_path / "out")
    assert results.job_id
    assert results.campaign["protein"]["uniprot"] == "P00533"
    assert results.summary["frames"] == 20
    assert not ingest.validate(egfr_archive)


def test_missing_summary_names_the_file(egfr_archive, tmp_path):
    broken = repack(egfr_archive, tmp_path / "broken.tar.gz", drop={"traj/summary.json"})
    with pytest.raises(ingest.IngestError) as caught:
        ingest.extract(broken, tmp_path / "out")
    assert "traj/summary.json" in str(caught.value)


def test_missing_several_files_lists_them(egfr_archive, tmp_path):
    broken = repack(egfr_archive, tmp_path / "broken.tar.gz",
                    drop={"complex_min.pdb", "poses/scores.csv"})
    with pytest.raises(ingest.IngestError) as caught:
        ingest.extract(broken, tmp_path / "out")
    message = str(caught.value)
    assert "complex_min.pdb" in message and "poses/scores.csv" in message


def test_future_schema_is_refused(egfr_archive, tmp_path):
    with tarfile.open(egfr_archive) as tar:
        manifest = json.loads(tar.extractfile("results/manifest.json").read())
    manifest["schema"] = "9.0"
    broken = repack(egfr_archive, tmp_path / "future.tar.gz",
                    replace={"manifest.json": json.dumps(manifest).encode()})
    with pytest.raises(ingest.IngestError) as caught:
        ingest.extract(broken, tmp_path / "out")
    assert "schema" in str(caught.value)


def test_manifest_that_is_not_an_object_is_refused(egfr_archive, tmp_path):
    # json.loads("3") is a valid int, and .get() on it raises AttributeError
    # rather than saying anything useful, so the type is checked explicitly.
    broken = repack(egfr_archive, tmp_path / "int.tar.gz",
                    replace={"manifest.json": b"3"})
    with pytest.raises(ingest.IngestError) as caught:
        ingest.extract(broken, tmp_path / "out")
    assert "manifest.json" in str(caught.value)


def test_extra_files_are_ignored(egfr_archive, tmp_path):
    with tarfile.open(egfr_archive) as src, tarfile.open(tmp_path / "extra.tar.gz", "w:gz") as out:
        for member in src.getmembers():
            fh = src.extractfile(member)
            if fh is not None:
                out.addfile(member, fh)
        data = b"nothing to see"
        info = tarfile.TarInfo("results/notes.txt")
        info.size = len(data)
        out.addfile(info, io.BytesIO(data))
    results = ingest.extract(tmp_path / "extra.tar.gz", tmp_path / "out")
    assert (results.root / "notes.txt").exists()
    assert results.summary["frames"] == 20


def test_a_repacked_archive_says_so(egfr_archive, tmp_path):
    """The manifest and the campaign both carry a job id. Disagreement means the
    archive was rebuilt after the run, which the results page should say."""
    with tarfile.open(egfr_archive) as tar:
        manifest = json.loads(tar.extractfile("results/manifest.json").read())
    manifest["job_id"] = "gs_20260905_somethingelse"
    repacked = repack(egfr_archive, tmp_path / "repacked.tar.gz",
                      replace={"manifest.json": json.dumps(manifest).encode()})
    results = ingest.extract(repacked, tmp_path / "out")
    assert results.job_id == "gs_20260905_somethingelse"
    assert any("repacked" in w for w in results.warnings)
