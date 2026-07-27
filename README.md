# Alignasm : Assembly tool for de novo genome alignment

Alignasm is part of the ACCtools pipeline, a project designed to resolve the complex karyotypes of cancer using graph-based approaches. Alignasm takes as input a PAF format file, where a *de novo* assembled genome has been aligned to a reference. Then, alignasm performs graph-based analysis to infer which references the contigs are likely derived from.

## Install
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

Use `--write-all` to additionally generate `input.aln.alt.paf`, which represents
the alternative paths, and `input.aln.all.paf`, which contains the extra
maximum-coverage paths tied with the selected path for the same alignment score
and anomaly count. These optional outputs are disabled by default because they
can be very large.
It is recommended to have more than 512 GB of RAM available for running alignasm.  
For a more detailed pipeline on how to use Alignasm, refer to [ACCtools-pipeline](https://github.com/ACCtools/ACCtools-pipeline).

## Author

Kyungmo Ku <pentagon03.codes@gmail.com>  
Hyunwoo Ryu <wowo0118@korea.ac.kr>
