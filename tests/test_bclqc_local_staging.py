"""Unit tests for local NVMe staging, local CRAM/FASTA resolution, and bounded QC in bclqc.py."""

import os
import sys
import types
from unittest import mock

# Stub plotting/interop imports that are only needed for InterOp PNG generation
for mod_name in (
    "matplotlib",
    "matplotlib.pyplot",
    "seaborn",
    "interop",
):
    if mod_name not in sys.modules:
        sys.modules[mod_name] = types.ModuleType(mod_name)
sys.modules["interop"].py_interop_run_metrics = mock.MagicMock()
sys.modules["interop"].py_interop_run = mock.MagicMock()
sys.modules["interop"].py_interop_table = mock.MagicMock()

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import bclqc


def test_resolve_local_bam_prefers_local_and_staged_normals(tmp_path):
    bams_dir = tmp_path / "run1"
    s1_dir = bams_dir / "S1"
    s1_dir.mkdir(parents=True)
    local_cram = s1_dir / "S1.cram"
    local_cram.write_bytes(b"cram")

    staged_dir = bams_dir / "_staged_normals" / "N1"
    staged_dir.mkdir(parents=True)
    staged_cram = staged_dir / "N1.cram"
    staged_cram.write_bytes(b"norm")

    assert bclqc.resolve_local_bam(str(bams_dir), "S1", "/mnt/pns/bams/run1/S1/S1.cram") == str(local_cram)
    assert bclqc.resolve_local_bam(str(bams_dir), "N1", "/mnt/pns/bams/old/N1/N1.cram") == str(staged_cram)
    assert bclqc.resolve_local_bam(str(bams_dir), "N2", "/mnt/pns/bams/old/N2/N2.cram") == "/mnt/pns/bams/old/N2/N2.cram"


def test_resolve_picard_ref(tmp_path):
    local_fa = tmp_path / "hg38.fa"
    local_fa.write_text(">chr1\nACGT\n")

    assert bclqc.resolve_picard_ref("/explicit/ref.fa") == "/explicit/ref.fa"
    with mock.patch.object(bclqc, "LOCAL_PICARD_REF", str(local_fa)):
        assert bclqc.resolve_picard_ref() == str(local_fa)
    with mock.patch.object(bclqc, "LOCAL_PICARD_REF", str(tmp_path / "missing.fa")):
        assert bclqc.resolve_picard_ref() == bclqc.PICARD_REF


def test_qcsum_command_passes_xmx_and_writes_hsm_atomically(tmp_path):
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    cmds = []

    def fake_exec(cmd, **kwargs):
        cmds.append(cmd)
        if cmd[0] == "picard":
            o_arg = next(a for a in cmd if a.startswith("O="))
            with open(o_arg[2:], "w") as fh:
                fh.write("# metrics\n")

    fake_cfg = {
        "bait_intervals": "/baits.ilist",
        "target_intervals": "/targets.ilist",
        "pipeline_version": "HemeDx-v2.2.0",
        "platform": "NovaSeq6000",
        "pass_min_align_pct": "99.0",
        "fail_min_align_pct": "99.0",
        "covered": "500",
        "pass_min_roi_pct": "85",
        "fail_min_roi_pct": "85",
        "pass_min_avgcov": "1070",
        "fail_min_avgcov": "1070",
        "pass_min_reads": "19000000",
        "fail_min_reads": "19000000",
        "capture": "goal",
        "capture_version": "v1.0",
    }

    with mock.patch.object(bclqc.QCSumInfo, "get_config", return_value=fake_cfg), \
         mock.patch.object(bclqc, "exec_command", side_effect=fake_exec):
        bclqc.qcsum_command("/local/S1.cram", "S1", str(sample_dir), "Comprehensive Heme Panel", picard_xmx="4g", picard_ref="/local/hg38.fa")

    assert cmds[0][:3] == ["picard", "-Xmx4g", "CollectHsMetrics"]
    assert "R=/local/hg38.fa" in cmds[0]
    assert (sample_dir / "S1.hsm.txt").exists()
    assert not (sample_dir / "S1.hsm.txt.tmp").exists()


def test_multiqc_ignores_staged_normals(tmp_path):
    cmds = []
    with mock.patch.object(bclqc, "exec_command", side_effect=lambda cmd: cmds.append(cmd)):
        bclqc.multiqc("/staging/hot/reads/run1", "/staging/hot/bams/run1")
    mqc_cmd = cmds[-1]
    assert mqc_cmd[0] == "multiqc"
    assert "--ignore" in mqc_cmd
    assert "_staged_normals" in mqc_cmd
