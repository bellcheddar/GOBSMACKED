"""The binding-mode classifiers, against structures whose labels are known.

These are the tests that would catch a subpocket definition or a distance
threshold drifting: 1M17 is a textbook type I, 1IEP a textbook type II, 2RH1 the
canonical inactive GPCR and 3SN6 the canonical active one.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from app.services import annotate, fetch, modes
from tests.conftest import needs_network

CACHE = Path(__file__).resolve().parent / "fixtures" / "_structures"


def structure(pdb_id: str) -> Path:
    """A cached mmCIF, downloaded once. Kept out of git: it is regenerable."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{pdb_id}.cif"
    if not path.exists():
        import requests

        path.write_bytes(requests.get(
            f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=60).content)
    return path


@pytest.fixture(scope="module")
def egfr():
    sequence = fetch.fetch_uniprot("P00533")["sequence"]
    return sequence, annotate.klifs_annotation("P00533", sequence, "EGFR")


@pytest.fixture(scope="module")
def abl():
    sequence = fetch.fetch_uniprot("P00519")["sequence"]
    return sequence, annotate.klifs_annotation("P00519", sequence, "ABL1")


@pytest.fixture(scope="module")
def adrb2():
    sequence = fetch.fetch_uniprot("P07550")["sequence"]
    return sequence, annotate.gpcrdb_annotation("P07550")


@needs_network
def test_klifs_pocket_maps_onto_the_canonical_sequence(egfr):
    sequence, klifs = egfr
    named = klifs["named_positions"]
    # Every one of these is checkable against the literature, which is the point.
    assert sequence[named["beta3_lysine"] - 1] == "K" and named["beta3_lysine"] == 745
    assert sequence[named["alphaC_glutamate"] - 1] == "E" and named["alphaC_glutamate"] == 762
    assert named["gatekeeper"] == 790                      # Thr790, the T790M site
    assert named["hinge"][2] == 793                        # Met793, the hinge H-bond
    assert named["dfg_asp"] == 855 and named["dfg_phe"] == 856
    assert len(klifs["pocket_map"]) >= 80


@needs_network
def test_1m17_is_type_one_dfg_in(egfr):
    sequence, klifs = egfr
    result = modes.classify_kinase(structure("1M17"), sequence, klifs["pocket_map"],
                                   chain_name="A", ligand_ccd="AQ4")
    assert result["label"] == "I"
    assert result["dfg"] == "in"
    assert result["alphac"] == "in"
    # 1M17 numbers EGFR from the mature protein, 24 lower than UniProt, and the
    # gatekeeper must follow the structure's own numbering.
    assert result["gatekeeper_residue"] == "Thr766"
    assert result["occupancy"]["hinge"] is True
    assert result["occupancy"]["back pocket II"] is False


@needs_network
def test_1iep_is_type_two_dfg_out(abl):
    sequence, klifs = abl
    result = modes.classify_kinase(structure("1IEP"), sequence, klifs["pocket_map"],
                                   chain_name="A", ligand_ccd="STI")
    assert result["label"] == "II"
    assert result["dfg"] == "out"
    assert result["gatekeeper_residue"] == "Thr315"        # the imatinib resistance site
    assert result["hinge_hbonds"] >= 1


@needs_network
def test_2rh1_is_orthosteric_and_inactive(adrb2):
    sequence, gpcrdb = adrb2
    result = modes.classify_gpcr(structure("2RH1"), sequence, gpcrdb["generic"],
                                 chain_name="A", ligand_ccd="CAU")
    assert result["label"] == "orthosteric"
    assert result["state"] == "inactive-like"
    assert result["toggle_residue"] == "Trp286"            # the 6.48 toggle switch
    assert "6.48" in result["contact_generic"]


@needs_network
def test_3sn6_is_active(adrb2):
    """The other side of the TM3-TM6 threshold, on the same receptor."""
    sequence, gpcrdb = adrb2
    result = modes.classify_gpcr(structure("3SN6"), sequence, gpcrdb["generic"],
                                 chain_name="R")
    assert result["state"] == "active-like"
    assert result["tm3_tm6"] > modes.TM3_TM6_CUT


@needs_network
def test_the_threshold_separates_the_two_populations(adrb2):
    """Six structures, measured by this code: the cut must sit in the gap."""
    sequence, gpcrdb = adrb2
    inactive = modes.classify_gpcr(structure("2RH1"), sequence, gpcrdb["generic"],
                                   chain_name="A", ligand_ccd="CAU")["tm3_tm6"]
    active = modes.classify_gpcr(structure("3SN6"), sequence, gpcrdb["generic"],
                                 chain_name="R")["tm3_tm6"]
    assert inactive < modes.TM3_TM6_CUT < active


def test_verdict_reports_a_match_and_a_difference():
    predicted = {"family": "kinase", "label": "I", "dfg": "in", "alphac": "in",
                 "occupancy": {"hinge": True}}
    same = dict(predicted)
    assert modes.compare_modes(predicted, same)["match"] is True
    other = {**predicted, "label": "II", "dfg": "out"}
    verdict = modes.compare_modes(predicted, other)
    assert verdict["match"] is False
    assert "DFG" in verdict["detail"]


def test_no_reference_is_unverified_not_a_mismatch():
    verdict = modes.compare_modes({"family": "kinase", "label": "I"}, None)
    assert verdict["match"] is None
    assert "unverified" in verdict["verdict"].lower()


@needs_network
def test_hinge_hbonds_counts_the_hinge_not_the_complex(egfr):
    """PLIP's H-bond count for the whole complex was being reported under a
    hinge label, so a pose that never touches the hinge showed three hinge
    hydrogen bonds."""
    sequence, klifs = egfr
    rows = [
        {"type": "hbond", "resnr": 769},        # Met769: the hinge, in 1M17 numbering
        {"type": "hbond", "resnr": 831},        # elsewhere in the pocket
        {"type": "hydrophobic", "resnr": 769},  # not a hydrogen bond
    ]
    result = modes.classify_kinase(structure("1M17"), sequence, klifs["pocket_map"],
                                   chain_name="A", ligand_ccd="AQ4", plip_rows=rows)
    assert 769 in result["hinge_residues"]
    assert result["hinge_hbonds"] == 1
    assert result["hinge_hbond_source"].startswith("PLIP")

    geometric = modes.classify_kinase(structure("1M17"), sequence, klifs["pocket_map"],
                                      chain_name="A", ligand_ccd="AQ4")
    assert geometric["hinge_hbond_source"].startswith("geometry")
