"""Behavioral tests using only Python's standard library and the built CLI."""

import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BINARY = Path(sys.argv.pop(1)).resolve()
FIXTURES = Path(__file__).resolve().parent / "fixtures"


def row(name, length, start, end, chromosome, ref_start, mapq, strand="+", cs=None, ref_length=None):
    span = end - start
    if ref_length is None:
        ref_length = span
    if cs is None:
        cs = f":{span}"
    # Tests supplying a longer reference span use an all-match alignment plus a deletion.
    fields = [name, length, start, end, strand, chromosome, 300000000,
              ref_start, ref_start + ref_length, span, max(span, ref_length), mapq,
              "tp:A:P", "cs:Z:" + cs]
    return "\t".join(map(str, fields)) + "\n"


class AlignasmCLI(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="alignasm-test-")
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def run_case(self, label, text, options=(), report=True, alternate=None, success=True):
        directory = self.root / label
        directory.mkdir()
        paf = directory / "input.paf"
        paf.write_text(text)
        command = [str(BINARY), str(paf), *map(str, options)]
        if alternate is not None:
            alt = directory / "alternate.paf"
            alt.write_text(alternate)
            command += ["--alt", str(alt)]
        if report:
            command += ["--score-report", str(directory / "scores.tsv")]
        process = subprocess.run(command, text=True, capture_output=True, timeout=60)
        if success:
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        else:
            self.assertNotEqual(process.returncode, 0)
        return directory, process

    def paf(self, directory, suffix="aln.paf"):
        return [line.split("\t") for line in (directory / f"input.{suffix}").read_text().splitlines()]

    def scores(self, directory):
        with (directory / "scores.tsv").open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_real_unitig_selects_primary_chr1_pair(self):
        source = (FIXTURES / "utg012958l.paf").read_text()
        directory, _ = self.run_case("real", source, ["--write-all"])
        raw = [line.split("\t") for line in source.splitlines()]
        self.assertEqual([r[:12] for r in self.paf(directory)], [raw[1][:12], raw[0][:12]])
        score = self.scores(directory)[0]
        self.assertEqual(score["total_scaled"], str(2004 * 60000))
        self.assertEqual(score["mapq_deficit"], "0")
        self.assertEqual(score["unmapped_query_bp"], "4")
        self.assertEqual(score["non_collinear_joins"], "1")
        self.assertEqual(score["paths_examined"], "1547")
        self.assertEqual(score["enumeration_limit_reached"], "0")
        fast, _ = self.run_case("real_fast", source)
        self.assertEqual((directory / "input.aln.paf").read_bytes(), (fast / "input.aln.paf").read_bytes())
        self.assertEqual(self.scores(fast)[0]["paths_examined"], "1")

    def test_legacy_matches_frozen_golden_outputs(self):
        source = (FIXTURES / "utg012958l.paf").read_text()
        directory, _ = self.run_case("legacy", source, ["--scoring", "legacy", "--write-all"], report=False)
        for suffix in ["aln.paf", "aln.alt.paf", "aln.all.paf"]:
            self.assertEqual((directory / f"input.{suffix}").read_bytes(),
                             (FIXTURES / f"utg012958l.legacy.{suffix}").read_bytes(), suffix)

    def test_high_quality_translocation_is_preserved(self):
        source = (row("trans", 2004, 0, 1000, "chr1", 10000, 60)
                  + row("trans", 2004, 1004, 2004, "chr7", 20000, 60)
                  + row("trans", 2004, 1000, 2004, "chr1", 20000, 0))
        directory, _ = self.run_case("trans", source)
        selected = self.paf(directory)
        self.assertEqual([r[5] for r in selected], ["chr1", "chr7"])
        self.assertEqual([r[11] for r in selected], ["60", "60"])
        self.assertEqual(self.scores(directory)[0]["total_scaled"], str(2004 * 60000))

    def test_overlap_clipping_accounts_once_on_both_strands(self):
        for strand in ["+", "-"]:
            positions = [1000, 1090] if strand == "+" else [1100, 1000]
            source = (row("overlap", 200, 0, 100, "chr1", positions[0], 60, strand)
                      + row("overlap", 200, 90, 200, "chr1", positions[1], 60, strand))
            directory, _ = self.run_case("overlap" + strand, source)
            selected = self.paf(directory)
            self.assertEqual([(int(r[2]), int(r[3])) for r in selected], [(0, 90), (90, 200)])
            self.assertEqual([int(r[9]) for r in selected], [90, 110])
            self.assertIn("cs:Z::90", selected[0])
            score = self.scores(directory)[0]
            self.assertEqual(score["total_scaled"], "0")
            self.assertEqual(score["mapq_deficit"], "0")
            self.assertEqual(score["pieces"], "2")

    def test_reference_span_does_not_override_ties(self):
        source = (row("span", 1000, 0, 1000, "chr1", 0, 0)
                  + row("span", 1000, 0, 1000, "chr2", 0, 0,
                        cs=":500-" + "a" * 500 + ":500", ref_length=1500))
        directory, _ = self.run_case("span", source, ["--write-all"])
        self.assertEqual(self.paf(directory)[0][5], "chr1")
        self.assertEqual(self.paf(directory, "aln.all.paf")[0][5], "chr2")
        self.assertEqual([r["path_kind"] for r in self.scores(directory)], ["selected", "tied"])
        fast, _ = self.run_case("span_fast", source)
        self.assertEqual(self.paf(directory), self.paf(fast))

    def test_zero_weight_disables_quality_tie_break(self):
        source = (row("zero", 1000, 0, 1000, "chr1", 0, 0)
                  + row("zero", 1000, 0, 1000, "chr2", 0, 60))
        weighted, _ = self.run_case("weighted", source)
        disabled, _ = self.run_case("disabled", source, ["--mapq-loss-per-kb", "0"])
        self.assertEqual(self.paf(weighted)[0][5], "chr2")
        self.assertEqual(self.paf(disabled)[0][5], "chr1")
        score = self.scores(disabled)[0]
        self.assertEqual(score["mapq_deficit"], "60000")
        self.assertEqual(score["mapq_loss_scaled"], "0")

    def test_unknown_and_capped_mapq_are_preserved_in_paf(self):
        for mapq, alternative in [(255, 0), (254, 60)]:
            source = (row("quality", 100, 0, 100, "chr1", 0, mapq)
                      + row("quality", 100, 0, 100, "chr2", 0, alternative))
            directory, _ = self.run_case("mapq" + str(mapq), source)
            self.assertEqual(self.paf(directory)[0][11], str(mapq))
            self.assertEqual(self.scores(directory)[0]["total_scaled"], "0")

    def test_single_alignment_keeps_terminal_gap_cost(self):
        directory, _ = self.run_case("single", row("single", 1000, 100, 900, "chr1", 2000, 7))
        score = self.scores(directory)[0]
        self.assertEqual(score["query_cost"], "400")
        self.assertEqual(score["unmapped_query_bp"], "200")
        self.assertEqual(score["mapq_deficit"], "0")
        self.assertEqual(score["pieces"], "1")

    def test_omitted_candidate_still_contributes_to_envelope(self):
        source = (row("omit", 10000, 0, 10000, "chr1", 0, 0)
                  + row("omit", 10000, 0, 1000, "chr2", 5000, 60))
        directory, _ = self.run_case("omit", source)
        self.assertEqual(self.paf(directory)[0][5], "chr1")
        self.assertEqual(self.scores(directory)[0]["mapq_loss_cost"], "10.000000")

    def test_alt_input_filter_and_fallback_define_envelope(self):
        primary = row("alt", 1000, 0, 1000, "chr1", 0, 0)
        alternative = (row("alt:1-1000", 1000, 0, 600, "chr2", 0, 60)
                       + row("alt:1-1000", 1000, 600, 1000, "chr2", 600, 60))
        directory, _ = self.run_case("alt_filter", primary, alternate=alternative)
        self.assertEqual(self.scores(directory)[0]["mapq_deficit"], "36000")
        directory, _ = self.run_case("alt_fallback", primary,
                                     alternate=row("alt:1-1000", 1000, 0, 400, "chr2", 0, 60))
        self.assertEqual(self.scores(directory)[0]["mapq_deficit"], "24000")
        directory, _ = self.run_case("alt_full", primary,
                                     alternate=row("alt:1-1000", 1000, 0, 1000, "chr2", 0, 60))
        self.assertEqual(self.paf(directory)[0][5], "chr2")
        self.assertIn("xi:Z:A_0", self.paf(directory)[0])
        directory, _ = self.run_case("alt_empty", primary, alternate="")
        self.assertEqual(self.scores(directory)[0]["mapq_deficit"], "0")

    def test_thread_and_enumeration_determinism(self):
        source = "".join(row(f"tie{i}", 100, 0, 100, chromosome, 0, 0)
                         for i in range(12) for chromosome in ["chr1", "chr2", "chr3"])
        serial, _ = self.run_case("serial", source, ["--write-all", "-t", "1"])
        for number in range(2):
            parallel, _ = self.run_case("parallel" + str(number), source, ["--write-all", "-t", "4"])
            for name in ["input.aln.paf", "input.aln.alt.paf", "input.aln.all.paf", "scores.tsv"]:
                self.assertEqual((serial / name).read_bytes(), (parallel / name).read_bytes(), name)

    def test_auxiliary_limit_does_not_change_optimum(self):
        source = "".join(row("many", 140, layer * 10, (layer + 1) * 10, "chr1", layer * 10, 0)
                         for layer in range(14) for _ in range(2))
        limited, process = self.run_case("limited", source, ["--write-all"])
        fast, _ = self.run_case("many_fast", source)
        self.assertEqual((limited / "input.aln.paf").read_bytes(), (fast / "input.aln.paf").read_bytes())
        scores = self.scores(limited)
        self.assertEqual(len(scores), 10000)
        self.assertEqual(scores[0]["paths_examined"], "10000")
        self.assertEqual(scores[0]["enumeration_limit_reached"], "1")
        self.assertIn("10000-path", process.stderr)

    def test_alternative_can_have_zero_total_cost_increase(self):
        source = (row("zero_up", 200000, 0, 100000, "chr1", 0, 60)
                  + row("zero_up", 200000, 100000, 200000, "chr2", 0, 60)
                  + row("zero_up", 200000, 0, 200000, "chr3", 0, 0))
        directory, _ = self.run_case("zero_up", source, ["--write-all"])
        self.assertEqual([r[11] for r in self.paf(directory)], ["60", "60"])
        self.assertEqual(self.paf(directory, "aln.alt.paf")[0][5], "chr3")
        scores = {r["path_kind"]: r for r in self.scores(directory)}
        self.assertEqual(scores["selected"]["total_scaled"], scores["alternative"]["total_scaled"])
        self.assertEqual(scores["alternative"]["non_collinear_joins"], "0")

    def test_invalid_options_and_overflow_fail_clearly(self):
        single = row("invalid", 100, 0, 100, "chr1", 0, 60)
        for index, options in enumerate([
            ["--scoring", "invalid"], ["--sv-cost", "0"], ["--sv-cost", "-1"],
            ["--mapq-loss-per-kb", "-1"], ["--mapq-loss-per-kb", "0.5"],
            ["--scoring", "legacy", "--sv-cost", "2000"],
            ["--scoring", "legacy", "--mapq-loss-per-kb", "10"],
        ]):
            self.run_case("invalid" + str(index), single, options, report=False, success=False)
        self.run_case("legacy_report", single, ["--scoring", "legacy"], success=False)
        source = (row("overflow", 200, 0, 100, "chr1", 0, 60)
                  + row("overflow", 200, 100, 200, "chr2", 0, 60))
        _, process = self.run_case("overflow", source, ["--sv-cost", str(2**63 - 1)], success=False)
        self.assertIn("64-bit", process.stderr)
        self.assertIn("overflow", process.stderr)

    def test_score_report_cannot_overwrite_input(self):
        path = self.root / "protected.paf"
        source = row("protected", 100, 0, 100, "chr1", 0, 60)
        path.write_text(source)
        process = subprocess.run([str(BINARY), str(path), "--score-report", str(path)],
                                 text=True, capture_output=True, timeout=10)
        self.assertNotEqual(process.returncode, 0)
        self.assertEqual(path.read_text(), source)


if __name__ == "__main__":
    unittest.main()
