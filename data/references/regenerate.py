#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Build OptiType reference (alleles.h5 + hla_reference_{dna,rna}.fasta) from an
IPD-IMGT/HLA release. Invoked by regenerate.sh; can also be run standalone.

Run on Python 2.7 inside the optitype-refgen-1.3.5 conda env. Imports OptiType
v1.3.5's hlatyper.py from the active env, applies in-memory patches for two
pandas-API regressions (irow -> iloc) and two O(N^2) hot loops, then runs:

  1. hlatyper.create_allele_dataframes(hla.dat, hla_gen.fasta, hla_nuc.fasta)
  2. restrict to canonical class-I + pseudogene loci {A,B,C,E,F,G,H,J,K,L,V}
  3. extract per-allele typing-window subsequences:
       DNA = intron1 + exon2 + intron2 + exon3 + intron3 (contiguous span)
       RNA = exon2 + exon3 (spliced)
  4. dedup FASTA records by sequence content (rep = lowest-ID member)
  5. write alleles.h5 keeping only `table` + `features` (the only tables
     OptiType reads at typing time)
"""
from __future__ import print_function
import argparse
import imp
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import OrderedDict

CANONICAL_LOCI = {"A","B","C","E","F","G","H","J","K","L","V"}


# ---------------------------------------------------------------------------
# hlatyper.py loading + monkeypatches
# ---------------------------------------------------------------------------

def load_patched_hlatyper():
    """
    Locate the OptiType v1.3.5 hlatyper.py from the active conda env, copy it
    to a temp file, apply two byte-level patches, and import.
    """
    opti = subprocess.check_output(["which", "OptiTypePipeline.py"]).strip()
    if not opti:
        sys.exit("OptiTypePipeline.py not on PATH; activate the optitype-refgen-1.3.5 env")
    src = os.path.join(os.path.dirname(opti.decode() if isinstance(opti, bytes) else opti),
                        "hlatyper.py")
    tmp = tempfile.NamedTemporaryFile(prefix="hlatyper_", suffix=".py", delete=False)
    with open(src, "rb") as fin:
        body = fin.read().decode("utf-8")
    # Patch 1: irow(0) was removed in pandas >= 0.20.
    body = body.replace(".irow(0)", ".iloc[0]")
    # Patch 2: vectorize the per-allele sanity-check loop (was O(N) with
    #          chained-assignment that triggers SettingWithCopyWarning).
    body = body.replace(
        "    for allele, features in joined.groupby('id'):\n"
        "        row = features.iloc[0]  # first row of the features subtable. Contains all allele information because of the join\n"
        "        sum_features_length = features['length'].sum()\n"
        "        sum_exons_length = features.loc[features['feature']=='exon']['length'].sum()\n"
        "        if row['len_gen']>0 and row['len_gen'] != sum_features_length:\n"
        "            if VERBOSE:\n"
        "                print(\"\\tFeature lengths don't add up to gen sequence length\", allele, row['len_gen'], sum_features_length, row['type'])\n"
        "            table.loc[allele]['flags'] += 4\n"
        "        if row['len_nuc']>0 and row['len_nuc'] != sum_exons_length:\n"
        "            if VERBOSE:\n"
        "                print(\"\\tExon lengths don't add up to nuc sequence length\", allele, row['len_nuc'], sum_exons_length, row['type'])\n"
        "            table.loc[allele]['flags'] += 8",
        "    grp_len_sum  = joined.groupby('id')['length'].sum()\n"
        "    grp_exon_sum = joined[joined['feature']=='exon'].groupby('id')['length'].sum()\n"
        "    head_rows    = joined.drop_duplicates(subset='id', keep='first').set_index('id')\n"
        "    for allele in head_rows.index:\n"
        "        len_gen = head_rows.at[allele, 'len_gen']\n"
        "        len_nuc = head_rows.at[allele, 'len_nuc']\n"
        "        atype   = head_rows.at[allele, 'type']\n"
        "        sum_features_length = grp_len_sum.get(allele, 0)\n"
        "        sum_exons_length    = grp_exon_sum.get(allele, 0)\n"
        "        if len_gen > 0 and len_gen != sum_features_length:\n"
        "            if VERBOSE:\n"
        "                print(\"\\tFeature lengths don't add up to gen sequence length\", allele, len_gen, sum_features_length, atype)\n"
        "            table.at[allele, 'flags'] = table.at[allele, 'flags'] + 4\n"
        "        if len_nuc > 0 and len_nuc != sum_exons_length:\n"
        "            if VERBOSE:\n"
        "                print(\"\\tExon lengths don't add up to nuc sequence length\", allele, len_nuc, sum_exons_length, atype)\n"
        "            table.at[allele, 'flags'] = table.at[allele, 'flags'] + 8"
    )
    # Patch 3: replace the O(N^2) boolean-mask scan in the feature-sequences
    #          loop with a dict lookup.
    body = body.replace(
        "    for i_id, i_features in all_features.groupby('id'):\n"
        "        seq = sequences.loc[(sequences['id']==i_id) & (sequences['source']=='dat')].iloc[0]['sequence']\n"
        "        for ft_idx, feature in i_features.iterrows():\n"
        "            ft_seq = seq[feature['start']:feature['end']]\n"
        "            all_ft_counter += 1\n"
        "            if ft_seq not in ft_seq_lookup:\n"
        "                ft_seq_lookup[ft_seq] = ft_counter\n"
        "                all_features.loc[ft_idx, 'seq_id'] = ft_counter\n"
        "                ft_counter += 1\n"
        "            else:\n"
        "                all_features.loc[ft_idx, 'seq_id'] = ft_seq_lookup[ft_seq]",
        "    _dat = sequences[sequences['source']=='dat']\n"
        "    dat_seq_by_id = dict(zip(_dat['id'].tolist(), _dat['sequence'].tolist()))\n"
        "    seq_id_col_pos = all_features.columns.get_loc('seq_id')\n"
        "    for i_id, i_features in all_features.groupby('id'):\n"
        "        seq = dat_seq_by_id[i_id]\n"
        "        for ft_idx, feature in i_features.iterrows():\n"
        "            ft_seq = seq[feature['start']:feature['end']]\n"
        "            all_ft_counter += 1\n"
        "            if ft_seq not in ft_seq_lookup:\n"
        "                ft_seq_lookup[ft_seq] = ft_counter\n"
        "                all_features.iat[ft_idx, seq_id_col_pos] = ft_counter\n"
        "                ft_counter += 1\n"
        "            else:\n"
        "                all_features.iat[ft_idx, seq_id_col_pos] = ft_seq_lookup[ft_seq]"
    )
    tmp.write(body.encode("utf-8"))
    tmp.close()
    return imp.load_source("hlatyper", tmp.name)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def find_input(imgt_dir, name):
    """IMGTHLA layout changed in 3.64 (fastas at repo root, hla_gen zipped)."""
    for cand in (os.path.join(imgt_dir, "fasta", name),
                 os.path.join(imgt_dir, name)):
        if os.path.exists(cand):
            return cand
    zp = os.path.join(imgt_dir, name + ".zip")
    if os.path.exists(zp):
        with zipfile.ZipFile(zp) as z:
            z.extractall(imgt_dir)
        if os.path.exists(os.path.join(imgt_dir, name)):
            return os.path.join(imgt_dir, name)
    sys.exit("could not locate %s under %s" % (name, imgt_dir))


def ensure_hla_dat(imgt_dir):
    p = os.path.join(imgt_dir, "hla.dat")
    if os.path.exists(p): return p
    zp = os.path.join(imgt_dir, "hla.dat.zip")
    if not os.path.exists(zp):
        sys.exit("hla.dat[.zip] not found in %s" % imgt_dir)
    with zipfile.ZipFile(zp) as z: z.extractall(imgt_dir)
    return p


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--imgt-dir",   required=True, help="path to IMGTHLA git checkout")
    ap.add_argument("--output-dir", required=True, help="output dir for FASTAs + h5")
    args = ap.parse_args()

    out = args.output_dir
    if not os.path.isdir(out):
        os.makedirs(out)

    import pandas as pd
    ht = load_patched_hlatyper()
    ht.VERBOSE = True

    print("[1/4] hlatyper.create_allele_dataframes ...")
    table, features, sequences, _ = ht.create_allele_dataframes(
        ensure_hla_dat(args.imgt_dir),
        find_input(args.imgt_dir, "hla_gen.fasta"),
        find_input(args.imgt_dir, "hla_nuc.fasta"),
    )

    print("[2/4] restrict to canonical loci %s" % sorted(CANONICAL_LOCI))
    t = table[table["locus"].isin(CANONICAL_LOCI)].copy()
    keep_ids = set(t["id"])
    features = features[features["id"].isin(keep_ids)].copy()
    sequences = sequences[sequences["id"].isin(keep_ids)].copy()
    print("    %d alleles after locus filter" % len(t))

    print("[3/4] extract typing-window subseqs + dedup ...")
    dat_by_id = dict(zip(sequences[sequences["source"]=="dat"]["id"].tolist(),
                         sequences[sequences["source"]=="dat"]["sequence"].tolist()))
    features_by_id = {aid: g for aid, g in features.groupby("id")}

    records = []  # (aid, atype, dna_or_None, rna)
    for _, row in t.iterrows():
        aid, atype = row["id"], row["type"]
        f = features_by_id.get(aid); dat = dat_by_id.get(aid)
        if f is None or dat is None: continue
        e2 = f[(f["feature"]=="exon") & (f["number"]==2)]
        e3 = f[(f["feature"]=="exon") & (f["number"]==3)]
        if not (len(e2) and len(e3)): continue
        rna = (dat[int(e2.iloc[0]["start"]):int(e2.iloc[0]["end"])] +
               dat[int(e3.iloc[0]["start"]):int(e3.iloc[0]["end"])])
        dna = None
        i1 = f[(f["feature"]=="intron") & (f["number"]==1)]
        i3 = f[(f["feature"]=="intron") & (f["number"]==3)]
        if len(i1) and len(i3):
            dna = dat[int(i1.iloc[0]["start"]):int(i3.iloc[0]["end"])]
        records.append((aid, atype, dna, rna))

    def dedup(field_index):
        by_seq = {}
        for rec in records:
            seq = rec[field_index]
            if seq is None: continue
            prev = by_seq.get(seq)
            if prev is None or rec[0] < prev[0]:
                by_seq[seq] = rec
        return by_seq

    dna_reps = dedup(2)
    rna_reps = dedup(3)
    print("    unique DNA: %d   unique RNA: %d" % (len(dna_reps), len(rna_reps)))

    def write_fasta(path, items, field_index):
        with open(path, "w") as fout:
            for rec in sorted(items.values(), key=lambda r: r[0]):
                aid, atype, _, _ = rec
                seq = rec[field_index]
                fout.write(">%s HLA-%s\n" % (aid, atype))
                for i in range(0, len(seq), 80):
                    fout.write(seq[i:i+80] + "\n")
    write_fasta(os.path.join(out, "hla_reference_dna.fasta"), dna_reps, 2)
    write_fasta(os.path.join(out, "hla_reference_rna.fasta"), rna_reps, 3)

    print("[4/4] write alleles.h5 (table + features only) ...")
    rep_ids = {r[0] for r in dna_reps.values()} | {r[0] for r in rna_reps.values()}
    table_kept    = t[t["id"].isin(rep_ids)].copy()
    features_kept = features[features["id"].isin(rep_ids)].copy()
    features_kept = features_kept.drop(
        columns=[c for c in ("seq_id","order") if c in features_kept.columns]
    )

    h5 = os.path.join(out, "alleles.h5")
    if os.path.exists(h5): os.remove(h5)
    s = pd.HDFStore(h5, complevel=9, complib="zlib")
    s["table"] = table_kept
    s["features"] = features_kept
    s.close()
    print("    h5 rows: table=%d, features=%d" % (len(table_kept), len(features_kept)))


if __name__ == "__main__":
    main()
