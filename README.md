# Alignasm : Assembly tool for de novo genome alignment

Alignasm is part of the ACCtools pipeline, a project designed to resolve the complex karyotypes of cancer using graph-based approaches. Alignasm takes as input a PAF format file, where a *de novo* assembled genome has been aligned to a reference. Then, alignasm performs graph-based analysis to infer which references the contigs are likely derived from.

## Install
Building requires CMake 3.27 or newer and a C++20 compiler. Use a compatible
C++ compiler/runtime for alignasm and its vcpkg dependencies.

```sh
git clone https://github.com/ACCtools/alignasm.git
cd alignasm

git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh

mkdir build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build

```

## Quick Start
```sh
alignasm <input.paf>
```

When provided with `input.paf`, alignasm generates `input.aln.paf`, which
contains the selected alignment paths.

**Quality scoring is the default.** To inspect the cost of each emitted path:

```sh
alignasm input.paf --score-report input.scores.tsv
```

For each query position, the baseline is the highest MAPQ among all admitted
input alignments covering that position. This includes `--alt` records after
the existing coverage filter. MAPQ is capped at 60 for scoring; 255 means
missing and contributes zero. The original MAPQ is preserved in the output.

The selected alignment supplies the quality at a mapped query position. An
unmapped position supplies zero, so dropping a high-quality interval still
incurs a quality loss. Positions without any candidate have a zero baseline.
Overlapping alignments are clipped at compatible `cs` match boundaries, and
each query position is counted once. Quality mode retains up to three cut
states per overlapping pair: the earliest valid cut, its reverse-query mirror
(the latest valid cut), and the best two-piece cut under the same scoring key.
All retained states participate in the joint path optimization; choosing a
pairwise-best cut alone can make adjacent cuts incompatible in a three-piece
chain. Equivalent cuts use a canonical, unordered reference-breakend key for
deterministic ordering. Legacy mode retains its original first-feasible cut.
Different complete paths may still tie on the whole comparison key, so this
does not establish unique breakpoint coordinates. `--write-all` exposes tied
paths within the documented enumeration limit.

```text
D = sum over query positions (baseline MAPQ - selected MAPQ)
MAPQ loss cost = K * D / 60000
total cost = reference connection cost + unmapped query cost + MAPQ loss cost
```

| Parameter or cost | Quality-mode default |
| --- | --- |
| `--sv-cost B` | 2000; positive integer |
| `--mapq-loss-per-kb K` | 10; nonnegative integer |
| Same chromosome, same strand | `min(abs(signed reference gap), B)` |
| Strand switch on the same chromosome | `B` |
| Different chromosomes | `B` |
| Internal unmapped query interval | 1 per bp |
| Unmapped query prefix or suffix | 2 per bp |

With K=10, replacing MAPQ 60 by MAPQ 0 over 24,800 bp costs 248. Replacing
MAPQ 7 by MAPQ 0 costs proportionally less. This is a length-weighted selection
penalty, not a claim that mapping errors at individual bases are independent.
K=0 disables both the MAPQ cost and its tie-break.

All comparisons use the exact integer numerator
`60000 * (reference cost + query cost) + K * D`. There is no per-piece rounding.
Equal total scores are ordered by lower MAPQ loss cost, fewer unmapped query
bases, fewer non-collinear joins, and fewer alignment pieces. A join counts
once when chromosome/strand differs or its signed reference gap differs from
its query gap. This is an alignment-boundary metric, not a count of complete
biological SV events. Remaining ties use a fixed candidate/enumeration order
for a fixed input row order, including across thread counts.

The quality mode optimizes this one objective throughout the existing DAG.
It does not run the legacy query-first gap upgrade or the final
query-plus-reference-span maximization. The reconstructed path is independently
rescored and checked against its graph distance. Invalid score arithmetic or
signed 64-bit overflow is reported as an error. Cut candidates are endpoints
of adjacent matched-block overlaps, or nearest matched-block gaps when an
adjacent cut is unavailable. The earliest, latest, and two-piece-optimal
candidates define a bounded, reverse-complement-closed set; the search does
not include every possible cut coordinate. Admitted alignments, `--alt`
filtering, mandatory overlap components, and `--non_skip_linkable` restrictions
are unchanged. A low-quality isolated component can therefore still force a
join; this cut correction does not change that separate graph restriction.

**Optional confident-coverage priority:**

```sh
alignasm input.paf --scoring coverage --coverage-mapq 20 --score-report scores.tsv
```

This mode addresses a tradeoff in the default scalar objective: omitting an
801 bp MAPQ-60 alignment can cost less than one extra connection at B=2000.
Coverage mode first minimizes the number of query bases that have an admitted
alignment at or above `--coverage-mapq` but are not covered by a selected
alignment at or above that threshold. It then uses the entire quality-mode
comparison key. The default threshold is 20; valid values are 1 through 60.
MAPQ 255 is missing and never qualifies; other values above 60 are capped.
`--coverage-mapq` is accepted only with `--scoring coverage`.

Coverage is the union of PAF query intervals, including internal alignment
gaps; it is not a count of matching bases or supporting molecules. Overlapping
or duplicated records do not multiply coverage. No input locus or individual
record is mandatory: an equally confident continuous alignment can explain
the same sequence with fewer connections. Bases for which no qualifying
alignment exists have no coverage priority. Low-MAPQ alignments retain the
usual costs. `--mapq-loss-per-kb 0` disables the secondary MAPQ penalty, not
the leading coverage criterion.

Both the main path and overlap cuts use this additive, lexicographic objective;
the recovered path is independently rescored. The admissible graph, including
mandatory components and its incomplete set of overlap cuts, is unchanged.
If qualifying intervals cannot coexist on a legal path, the reported minimum
`confident_coverage_loss_bp` may be positive. Alternative paths with fewer
joins must preserve this minimum before their scalar cost is compared.

This is an opt-in alignment policy, not a calibrated SV confidence model.
High MAPQ does not rule out incorrect reference placement or a chimeric
assembly, and even a very short qualifying interval can change the selected
path. Independent junction evidence and a broader accuracy benchmark are
needed before adopting this policy as a pipeline default. `quality` remains
the default, and historical `legacy` selection is unchanged.

**Legacy compatibility:**

```sh
alignasm input.paf --scoring legacy
```

This retains the scoring, gap upgrades, and final selection of commit
`b234c268`, including the old reference-span tie-break. Explicit `--sv-cost`,
`--mapq-loss-per-kb`, and `--score-report` options are rejected in legacy mode.

**Optional alternative outputs:** `--write-all` additionally writes
`input.aln.all.paf` and `input.aln.alt.paf`. These can be large and are disabled
by default.

- In quality mode, `.aln.all.paf` contains additional paths tied on the entire
  comparison key above. Their query names retain the `.1`, `.2`, ... suffixes.
  Different compatible cut states can contribute additional tied paths, so
  enumeration counts need not match the earlier first-feasible-cut policy.
  `.aln.alt.paf` contains a path with fewer non-collinear joins minimizing
  `increase in total cost / decrease in joins`, with the common comparison key
  breaking ties. A zero cost increase is allowed.
- Coverage mode uses the same outputs with confident-coverage loss preceding
  the quality comparison key. A path sacrificing confident coverage is not
  emitted as a simpler `.aln.alt.paf` alternative.
- Auxiliary paths are considered within the first 10,000 enumerated paths.
  Reaching that limit is reported; it does not affect the exact selected
  optimum within the graph. A limit notification alone does not establish
  whether more paths exist.
- In legacy mode, `.aln.all.paf` retains the previous maximum-coverage ties at
  the same legacy score and anomaly count; `.aln.alt.paf` retains its previous
  behavior.

The optional score TSV records query and path kind/number, B and K,
`total_scaled` and `score_scale=60000`, reference/query costs, raw MAPQ
deficit, scaled MAPQ loss, unmapped bp, joins, pieces, paths examined, and
whether the auxiliary enumeration limit was reached. It also records
`coverage_mapq` and `confident_coverage_loss_bp` (both `.` outside coverage
mode, where confident coverage is not measured). In coverage mode `total_cost` is the secondary
scalar score, not the full path ranking. Decimal display columns
are rounded to six places; the integer columns are authoritative. The TSV
path must differ from PAF input/output paths. PAF filenames, basic fields,
`xi` provenance and clipped `cs` representation keep their existing format.

**Tests:** CTest runs C++ score/graph tests and Python standard-library CLI
tests, including the frozen H1437 fixture and legacy outputs.

```sh
cmake --build build --target alignasm quality_scoring_test
ctest --test-dir build --output-on-failure
```

Tests are enabled by default and need Python 3. Configure with
`-DBUILD_TESTING=OFF` when building only the executable.

It is recommended to have more than 512 GB of RAM available for running alignasm.  
For a more detailed pipeline on how to use Alignasm, refer to [ACCtools-pipeline](https://github.com/ACCtools/ACCtools-pipeline).

## Author

Kyungmo Ku <pentagon03.codes@gmail.com>  
Hyunwoo Ryu <wowo0118@korea.ac.kr>
