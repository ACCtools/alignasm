# H1437 regression fixture

`utg012958l.paf` contains the 45 records for this unitig extracted, in original
order, from `/Data/hyunwoo/00_skype_run_data/H1437/20_alignasm/H1437.utg.paf`.
Its SHA256 is `2dabf77dc6c3e9174db5c466c5ccb0e634649c8d216765d5c1fb0b3116e0ec47`.

The three `utg012958l.legacy.*.paf` files were generated from this exact subset
with commit `b234c268f1bb748ed895383327c373e1f7c70120`, using `--write-all`.
The baseline was compiled with GCC 11.5, CMake 3.31.10, and the pre-existing
CSV/argparse/vcpkg dependencies. `xi` indices therefore refer to this 45-row
fixture, not row numbers in the full sample input.

The quality-mode regression expects the original chr1 alignments with MAPQ 7
and 60, retaining the four-base query gap, at total cost 2,004. This checks the
specified scoring behavior; it is not a general biological truth set.
