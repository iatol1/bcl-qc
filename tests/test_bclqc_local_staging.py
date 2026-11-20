"""Unit tests for local NVMe staging, local CRAM/FASTA resolution, and bounded QC in bclqc.py."""

import os
import sys
import types
from unittest import mock

import pytest

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


def test_resolve_local_bam_prefers_nonempty_local_and_handles_missing_fallback(tmp_path):
    bams_dir = tmp_path / "run1"
    s1_dir = bams_dir / "S1"
    s1_dir.mkdir(parents=True)
    local_cram = s1_dir / "S1.cram"
    local_cram.write_bytes(b"cram")

    # 0-byte local CRAM must be ignored in favor of fallback
    empty_dir = bams_dir / "S_EMPTY"
    empty_dir.mkdir(parents=True)
    (empty_dir / "S_EMPTY.cram").write_bytes(b"")

    staged_dir = bams_dir / "_staged_normals" / "N1"
    staged_dir.mkdir(parents=True)
    staged_cram = staged_dir / "N1.cram"
    staged_cram.write_bytes(b"norm")

    assert bclqc.resolve_local_bam(str(bams_dir), "S1", "/mnt/pns/bams/run1/S1/S1.cram") == str(local_cram)
    assert bclqc.resolve_local_bam(str(bams_dir), "S_EMPTY", "/mnt/pns/bams/run1/S_EMPTY/S_EMPTY.cram") == "/mnt/pns/bams/run1/S_EMPTY/S_EMPTY.cram"
    assert bclqc.resolve_local_bam(str(bams_dir), "N1", "/mnt/pns/bams/old/N1/N1.cram") == str(staged_cram)
    assert bclqc.resolve_local_bam(str(bams_dir), "N2", "/mnt/pns/bams/old/N2/N2.cram") == "/mnt/pns/bams/old/N2/N2.cram"
    assert bclqc.resolve_local_bam(str(bams_dir), "N3", None) is None
    assert bclqc.resolve_local_bam(str(bams_dir), "N3", float("nan")) is None
    assert bclqc.resolve_local_bam(str(bams_dir), "N3", "   ") is None


def test_resolve_picard_ref_requires_fai_and_dict(tmp_path):
    local_fa = tmp_path / "hg38.fa"
    local_fa.write_text(">chr1\nACGT\n")

    assert bclqc.resolve_picard_ref("/explicit/ref.fa") == "/explicit/ref.fa"

    # Missing .fai and .dict -> must fall back to NAS PICARD_REF
    with mock.patch.object(bclqc, "LOCAL_PICARD_REF", str(local_fa)):
        assert bclqc.resolve_picard_ref() == bclqc.PICARD_REF

    # Only .fai present -> still must fall back
    (tmp_path / "hg38.fa.fai").write_text("chr1\t4\t6\t4\t5\n")
    with mock.patch.object(bclqc, "LOCAL_PICARD_REF", str(local_fa)):
        assert bclqc.resolve_picard_ref() == bclqc.PICARD_REF

    # Both .fai and <stem>.dict present -> local FASTA selected
    stem_dict = tmp_path / "hg38.dict"
    stem_dict.write_text("@HD\tVN:1.6\n")
    with mock.patch.object(bclqc, "LOCAL_PICARD_REF", str(local_fa)):
        assert bclqc.resolve_picard_ref() == str(local_fa)

    # Alternative <ref>.dict naming also accepted
    stem_dict.unlink()
    (tmp_path / "hg38.fa.dict").write_text("@HD\tVN:1.6\n")
    with mock.patch.object(bclqc, "LOCAL_PICARD_REF", str(local_fa)):
        assert bclqc.resolve_picard_ref() == str(local_fa)

    with mock.patch.object(bclqc, "LOCAL_PICARD_REF", str(tmp_path / "missing.fa")):
        assert bclqc.resolve_picard_ref() == bclqc.PICARD_REF


def test_has_enough_staging_space_handles_relative_and_missing_mounts(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Non-existent relative child under existing cwd
    assert bclqc.has_enough_staging_space("sub/local_bams", 0.0) is True
    assert bclqc.has_enough_staging_space("sub/local_bams", 1e12) is False
    # Unmounted top-level path walking all the way up to '/' must return False
    assert bclqc.has_enough_staging_space("/definitely_unmounted_nvme_xyz/hot/bams", 0.0) is False


def test_qcsum_command_passes_xmx_normalized_and_writes_hsm_atomically(tmp_path):
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
        bclqc.qcsum_command("/local/S1.cram", "S1", str(sample_dir), "Comprehensive Heme Panel", picard_xmx="-Xmx4g", picard_ref="/local/hg38.fa")

    assert cmds[0][:3] == ["picard", "-Xmx4g", "CollectHsMetrics"]
    assert "R=/local/hg38.fa" in cmds[0]
    assert (sample_dir / "S1.hsm.txt").exists()
    assert not (sample_dir / "S1.hsm.txt.tmp").exists()


def test_qcsum_bounded_pool_atomic_csv_and_missing_bam_error(tmp_path):
    bams_dir = tmp_path / "bams"
    s1_dir = bams_dir / "S1"
    s2_dir = bams_dir / "S2"
    s1_dir.mkdir(parents=True)
    s2_dir.mkdir(parents=True)
    (s1_dir / "S1.cram").write_bytes(b"c1")
    (s2_dir / "S2.cram").write_bytes(b"c2")

    si = tmp_path / "sampleinfo.tsv"
    si.write_text(
        "Samples\tPanel\tBAM Path\n"
        "S1\tComprehensive Heme Panel\t/nas/S1.cram\n"
        "S2\tComprehensive Heme Panel\t/nas/S2.cram\n"
    )

    pool_sizes = []

    class FakePool:
        def __init__(self, processes):
            pool_sizes.append(processes)

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def starmap(self, fn, args_iter):
            for bam_path, sample_id, output_dir, *_ in args_iter:
                with open(os.path.join(output_dir, f"{sample_id}.qcsum.txt"), "w") as fh:
                    fh.write("header\n")
                    fh.write(f"{sample_id},,NovaSeq6000,v2,PASS,PASS\n")

    with mock.patch.object(bclqc, "Pool", FakePool):
        bclqc.qcsum(str(bams_dir), str(si), qc_workers=4, picard_xmx="4g")

    assert pool_sizes == [2]  # min(qc_workers=4, num_samples=2)
    csv_path = bams_dir / "qcsum_mqc.csv"
    assert csv_path.exists()
    assert not (bams_dir / "qcsum_mqc.csv.tmp").exists()
    content = csv_path.read_text()
    assert "S1,,NovaSeq6000,v2,PASS,PASS" in content
    assert "S2,,NovaSeq6000,v2,PASS,PASS" in content

    # Missing BAM Path without local CRAM raises ValueError instead of TypeError
    bad_si = tmp_path / "bad_si.tsv"
    bad_si.write_text("Samples\tPanel\tBAM Path\nS_MISSING\tComprehensive Heme Panel\t\n")
    with pytest.raises(ValueError, match="S_MISSING"):
        bclqc.qcsum(str(bams_dir), str(bad_si), qc_workers=2)


def test_align_samples_requires_cram_and_crai_and_uses_force_on_partial(tmp_path):
    fastqs_dir = tmp_path / "reads" / "run1"
    i10_dir = fastqs_dir / "I10"
    i10_dir.mkdir(parents=True)
    (i10_dir / "fastq_list.csv").write_text(
        "RGID,RGSM,RGLB,Lane,Read1File,Read2File\n"
        "1,DONE,L1,1,r1,r2\n"
        "2,PARTIAL,L1,1,r1,r2\n"
        "3,FRESH,L1,1,r1,r2\n"
    )
    bams_dir = tmp_path / "bams" / "run1"
    # DONE has non-empty CRAM + CRAI -> skipped
    done_dir = bams_dir / "DONE"
    done_dir.mkdir(parents=True)
    (done_dir / "DONE.cram").write_bytes(b"cram")
    (done_dir / "DONE.cram.crai").write_bytes(b"crai")
    # PARTIAL has 0-byte CRAM without CRAI -> must re-run with force=True
    partial_dir = bams_dir / "PARTIAL"
    partial_dir.mkdir(parents=True)
    (partial_dir / "PARTIAL.cram").write_bytes(b"")

    calls = []
    with mock.patch.object(bclqc, "align", side_effect=lambda fl, out, sid, panel, force=False: calls.append((sid, panel, force))):
        bclqc.align_samples(str(fastqs_dir), str(bams_dir))

    assert calls == [
        ("FRESH", "PCP", False),
        ("PARTIAL", "PCP", True),
    ]


def test_qc_samples_creates_bams_dir_before_plot(tmp_path):
    bams_dir = tmp_path / "new_bams_dir"
    seen_exists = []

    def fake_plot(run_dir, out_dir):
        seen_exists.append(os.path.isdir(out_dir))

    with mock.patch.object(bclqc, "save_occ_pf_plot", side_effect=fake_plot), \
         mock.patch.object(bclqc, "qcsum"), \
         mock.patch.object(bclqc, "multiqc"):
        bclqc.qc_samples(str(tmp_path / "run"), str(tmp_path / "reads"), str(bams_dir), str(tmp_path / "si.tsv"))

    assert seen_exists == [True]


def test_sync_local_bams_to_nas_guards_missing_and_uses_cifs_flags(tmp_path):
    cmds = []
    with mock.patch.object(bclqc, "exec_command", side_effect=lambda cmd: cmds.append(cmd)):
        # Missing local dir (e.g. exome-only run) must no-op without invoking rsync
        bclqc.sync_local_bams_to_nas(str(tmp_path / "nonexistent"), str(tmp_path / "nas"))
        assert cmds == []

        local_dir = tmp_path / "local"
        local_dir.mkdir()
        bclqc.sync_local_bams_to_nas(str(local_dir), str(tmp_path / "nas"))

    assert len(cmds) == 1
    cmd = cmds[0]
    assert cmd[:2] == ["rsync", "-rlt"]
    assert "-a" not in cmd
    for flag in ("--no-perms", "--no-owner", "--no-group", "--omit-dir-times"):
        assert flag in cmd
    assert ".*.done" in cmd
    assert "_staged_normals" in cmd
    assert "*.tmp" in cmd


def test_multiqc_ignores_staged_normals(tmp_path):
    cmds = []
    with mock.patch.object(bclqc, "exec_command", side_effect=lambda cmd: cmds.append(cmd)):
        bclqc.multiqc("/staging/hot/reads/run1", "/staging/hot/bams/run1")
    mqc_cmd = cmds[-1]
    assert mqc_cmd[0] == "multiqc"
    assert "--ignore" in mqc_cmd
    assert "_staged_normals" in mqc_cmd


def test_bclqc_run_stage_local_and_no_sync_behavior(tmp_path):
    run_dir = tmp_path / "250101_A01718_0001_AH00000001"
    run_dir.mkdir()
    local_root = tmp_path / "local_bams"
    nas_root = tmp_path / "nas_bams"

    argv = [
        "bclqc.py",
        "--run-dir", str(run_dir) + "//",
        "--fastqs-dir", str(tmp_path / "reads"),
        "--bams-dir", str(nas_root),
        "--local-bams-dir", str(local_root),
        "--stage-local",
        "--staging-min-gb", "0",
        "--steps", "align",
    ]
    with mock.patch.object(sys, "argv", argv), \
         mock.patch.object(bclqc, "align_samples") as m_align, \
         mock.patch.object(bclqc, "sync_local_bams_to_nas") as m_sync:
        bclqc.bclqc_run()

    expected_local = str(local_root / "250101_A01718_0001_AH00000001")
    expected_nas = str(nas_root / "250101_A01718_0001_AH00000001")
    m_align.assert_called_once_with(str(tmp_path / "reads" / "250101_A01718_0001_AH00000001"), expected_local)
    m_sync.assert_called_once_with(expected_local, expected_nas)

    # With --no-sync, sync_local_bams_to_nas must not be called
    with mock.patch.object(sys, "argv", argv + ["--no-sync"]), \
         mock.patch.object(bclqc, "align_samples"), \
         mock.patch.object(bclqc, "sync_local_bams_to_nas") as m_sync2:
        bclqc.bclqc_run()
    m_sync2.assert_not_called()
