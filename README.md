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
```
alignasm <input.paf>
```

When provided with input.paf as input, it generates two serialized outputs: `input.aln.paf`, which contains the aligned data, and `input.aln.alt.paf`, which represents the alternative paths.
It is recommended to have more than 512 GB of RAM available for running alignasm.  
For a more detailed pipeline on how to use Alignasm, refer to [ACCtools-pipeline](https://github.com/ACCtools/ACCtools-pipeline).

## Author

Kyungmo Ku <kyungmoku7141@gmail.com>  
Hyunwoo Ryu <wowo0118@korea.ac.kr>