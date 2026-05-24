# OptiType reference data

Files consumed by the pipeline:

| File | Used by | Content |
|---|---|---|
| `hla_reference_dna.fasta` | YARA (DNA mode) | per allele: `intron1 + exon2 + intron2 + exon3 + intron3` (contiguous) |
| `hla_reference_rna.fasta` | YARA (RNA mode) | per allele: spliced `exon2 + exon3` |
| `alleles.h5`              | OPTITYPE        | pandas HDFStore with `table` + `features` (the only keys OptiType reads at typing time) |

All three are restricted to canonical class-I + pseudogene loci
`{A, B, C, E, F, G, H, J, K, L}` and deduplicated by sequence content (one
FASTA record per unique sequence; rep = lowest-ID allele in the cluster).
Provenance — IMGT release, source URL, per-file MD5s — is recorded in
[`assets/software_meta.json`](../../assets/software_meta.json).

## Regenerating from a new IMGT release

```bash
# default IMGT_TAG is v3.63.1-alpha
bash data/references/regenerate.sh [IMGT_TAG]
```

Run from the repo root. Requirements:

- `mamba` (or `conda`) — used to create the pinned conda env
  `optitype-refgen-1.3.5` on first run (Python 2.7 + biopython 1.70 +
  pandas 0.24 + pytables 3.5; matches upstream OptiType v1.3.5)
- `git`, `unzip`
- ~6 GB free disk for the IMGTHLA clone + intermediate output

The script is idempotent: re-runs reuse the existing conda env and a cached
IMGTHLA clone (set `WORKDIR=/path/to/keep` to keep the clone across runs).

Expected wall time: ~5 min on the first run (mostly the IMGT EMBL parse +
HDF5 write), seconds on re-runs that hit the conda-env cache.

After running, copy the three printed MD5s into `assets/software_meta.json`
under the `optitype` block and bump `imgt_release` / `imgt_release_date` /
`imgt_source` accordingly.

## What the build does

`regenerate.sh` is a thin wrapper. `regenerate.py` does the work:

1. Import OptiType v1.3.5's `hlatyper.py` from the conda env. Apply two
   small in-memory patches that fix a pandas-API regression
   (`irow` → `iloc`) and two `O(n²)` hot loops that hurt on the larger
   modern IMGT releases.
2. Run `hlatyper.create_allele_dataframes(hla.dat, hla_gen.fasta, hla_nuc.fasta)`.
3. Restrict to the canonical 11 loci listed above.
4. For each allele, extract its DNA and RNA typing windows.
5. Deduplicate each FASTA by sequence content.
6. Write `alleles.h5` containing only `table` + `features`.

## Why this scope (and not raw IMGT)

OptiType only types class I on exons 2+3 (Szolek et al. 2014). The pipeline
ships only what's needed:

- the FASTAs are **trimmed** to the typing window — shipping full-gene
  FASTAs would multiply yara mapping and OptiType ILP runtime several-fold
  for no typing benefit;
- `alleles.h5` is **stripped** to `table` + `features` only — the bundled
  `sequences` and `feature_sequences` tables (~180 MB pre-trim) are never
  read by OptiType at typing time.
