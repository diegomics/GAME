#!/usr/bin/env python3

# Sex chromosome and sample sex configuration handling
# by Diego De Panis, 2026
# This script is part of the GAME pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

"""
sex_helpers.py

Provides parsing, validation, and accessors for two config fields:

  - asm_sex_chr      (assembly-level): names of contigs that are sex
                                   chromosomes in this reference assembly
  - sample_sex_karyo (sample-level):   karyotype of this sample, drawn from
                                   a controlled vocabulary

Empty-like values ("", None, "None", "NA", "-", etc.) are normalized
silently to [] / "unknown".
"""

import re
import sys


# -------------------------------------------------------------------------------
# EMPTY-VALUE NORMALIZATION
# -------------------------------------------------------------------------------
# Many YAML inputs encode "no value" in many ways. Treat them all the same.

_EMPTY_STRINGS = {"none", "na", "n/a", "-", ".", "null", "nan"}


def _is_empty(v):
    """True if v should be treated as 'no value provided'."""
    if v is None:
        return True
    if isinstance(v, str) and (v.strip() == "" or v.strip().lower() in _EMPTY_STRINGS):
        return True
    return False


# -------------------------------------------------------------------------------
# CONTROLLED VOCABULARY FOR sample_sex
# -------------------------------------------------------------------------------
# HETEROGAMETIC means that all contigs listed in sex_chr should be
# called as haploid for this sample. HOMOGAMETIC means that those
# contigs are diploid (no --haploid_contigs flag). "unknown" triggers the
# coverage-based inference path.

HETEROGAMETIC = {"ZW", "XY", "X0", "Z0", "XY1Y2", "X1X2Y"}
HOMOGAMETIC   = {"ZZ", "XX", "X1X1X2X2", "Z1Z1Z2Z2"}
_VALID_SEX    = HETEROGAMETIC | HOMOGAMETIC | {"unknown"}


def parse_sample_sex(value):
    """
    Normalize a raw sample_sex value to the controlled vocabulary.

    Returns one of the entries in HETEROGAMETIC, HOMOGAMETIC, or "unknown".
    Empty-like values become "unknown" silently.
    Raises ValueError if value is non-empty but not in the vocabulary.
    """
    if _is_empty(value):
        return "unknown"
    s = str(value).strip().upper()
    if s in _VALID_SEX:
        return s
    # Tolerate a few common spellings of "unknown"
    if s in {"UNKNOWN", "UNK", "?"}:
        return "unknown"
    raise ValueError(
        f"Invalid sample_sex value: {value!r}. "
        f"Allowed: heterogametic {sorted(HETEROGAMETIC)}, "
        f"homogametic {sorted(HOMOGAMETIC)}, "
        f"or 'unknown' / empty."
    )


def is_heterogametic(sex):
    """True if this sample's sex contigs should be called haploid."""
    return sex in HETEROGAMETIC


def is_homogametic(sex):
    """True if this sample's sex contigs should be called diploid."""
    return sex in HOMOGAMETIC


# -------------------------------------------------------------------------------
# PARSING sex_chr
# -------------------------------------------------------------------------------
# Accepted forms:
#   - missing / None / "" / "NA" / "-"           -> []
#   - "CM098573.1, CM098574.1"                   -> ["CM098573.1", "CM098574.1"]
#   - ["CM098573.1", "CM098574.1"]               -> same
#   - "CM098573.1; CM098574.1" or whitespace     -> same (commas, semis, spaces)

_SEX_CHR_SPLIT = re.compile(r"[,\s;]+")


def parse_sex_chr(value):
    """Normalize a raw sex_chr value to a list of contig names. May be empty."""
    if _is_empty(value):
        return []
    if isinstance(value, list):
        return [str(x).strip() for x in value if not _is_empty(x) and str(x).strip()]
    return [t.strip() for t in _SEX_CHR_SPLIT.split(str(value)) if t.strip()]


# -------------------------------------------------------------------------------
# CONFIG ACCESSORS
# -------------------------------------------------------------------------------
# return safe defaults instead of raising on missing keys

def get_sex_chr(samples_config, sp, asm):
    """Return the list of sex contigs declared for this assembly. May be []."""
    try:
        raw = samples_config["sp_name"][sp]["asm_id"][asm].get("asm_sex_chr")
    except (KeyError, TypeError, AttributeError):
        return []
    return parse_sex_chr(raw)


def get_sample_sex(samples_config, sp, asm, sid):
    """Return the parsed sample_sex value. Always a valid vocabulary string."""
    try:
        raw = (samples_config["sp_name"][sp]["asm_id"][asm]
               ["sample_id"][sid].get("sample_sex_karyo"))
    except (KeyError, TypeError, AttributeError):
        return "unknown"
    try:
        return parse_sample_sex(raw)
    except ValueError:
        # I think the validation shooould halt the pipeline before we get here.
        # If we somehow do (e.g. helper called outside the validator path),
        # degrade safely to "unknown" rather than crashing in a rule.
        return "unknown"


# -------------------------------------------------------------------------------
# STARTUP VALIDATOR
# -------------------------------------------------------------------------------
# TO CHECK!! Called once from the Snakefile right after samples_config is loaded.
# Walks the whole config, collects ALL errors before exiting (friendlier than
# halting on the first one when a config has many samples). Info-level notes
# are printed but do not halt.

def validate_sex_config(samples_config, is_main=True):
    """
    Validate sex_chr and sample_sex across the whole config.

    Halts the pipeline (sys.exit(1)) on:
      - sample_sex value outside the controlled vocabulary

    Emits an info message (does not halt) on:
      - sample with heterogametic sample_sex but no sex_chr on its assembly
        (typical case: draft reference with unresolved sex chromosomes)
    """
    errors = []
    info_msgs = []

    for sp, sd in (samples_config.get("sp_name") or {}).items():
        if not sd or "asm_id" not in sd:
            continue
        for asm, ad in (sd["asm_id"] or {}).items():
            if not ad:
                continue

            sex_chr_list = get_sex_chr(samples_config, sp, asm)

            for sid, sdata in (ad.get("sample_id") or {}).items():
                if sid in (None, "None"):
                    continue
                raw_sex = (sdata or {}).get("sample_sex_karyo")
                try:
                    sex = parse_sample_sex(raw_sex)
                except ValueError as e:
                    errors.append(
                        f"  sp={sp}, asm={asm}, sample={sid}: {e}"
                    )
                    continue

                # Info note: declared heterogametic but no sex_chr on the asm
                if is_heterogametic(sex) and not sex_chr_list:
                    info_msgs.append(
                        f"  sp={sp}, asm={asm}, sample={sid}: "
                        f"sample_sex_karyo={sex} but no sex_chr declared on this "
                        f"assembly (reference may not have resolved sex "
                        f"chromosomes). Sex contigs cannot be flagged as "
                        f"haploid; sample_sex is preserved for downstream "
                        f"analyses."
                    )

    if info_msgs and is_main:
        print("[GAME] asm_sex_chr / sample_sex_karyo notes:")
        for m in info_msgs:
            print(m)

    if errors:
        print("❌ Invalid sex configuration:", file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(1)


# -------------------------------------------------------------------------------
# HAPLOID-CONTIG RESOLVER
# -------------------------------------------------------------------------------
# Per-sample dispatcher: combines sex_chr (assembly-level) and sample_sex
# (sample-level) into the list of contigs that DeepVariant should call as
# haploid for this sample, plus a `source` tag explaining the decision.
#
#  For the `inferred` case:
# (sex_chr populated but sample_sex == "unknown"), the actual coverage-based
# classification lives in the standalone inference script and its output
# (sex_inference.tsv) is read by the caller (the DeepVariant input function)
#
# Return: (haploid_contigs, source)
#   haploid_contigs : list[str] or None
#                     - list of contig names to pass to --haploid_contigs
#                     - None  -> caller must consult sex_inference.tsv
#                     - []    -> no flag should be passed
#   source          : one of:
#                     - "declared"          haploid set determined directly
#                                           from sample_sex (heterogametic or
#                                           homogametic). Includes the case
#                                           where the set is empty because
#                                           the sample is homogametic.
#                     - "inferred"          sample_sex is unknown; caller must
#                                           read sex_inference.tsv to get the
#                                           per-contig classifications.
#                     - "missing_sex_chr"   sample is heterogametic but the
#                                           assembly has no sex_chr declared.
#                                           No flag passed; info already logged
#                                           by validate_sex_config().
#                     - "none"              no info available; nothing to do.

def resolve_haploid_contigs(samples_config, sp, asm, sid):
    """
    Return (haploid_contigs, source) for one sample.

    See module docstring for the 6-case matrix this implements.
    """
    sex_chr = get_sex_chr(samples_config, sp, asm)
    sex     = get_sample_sex(samples_config, sp, asm, sid)

    if not sex_chr:
        # No sex contigs declared on this assembly.
        if is_heterogametic(sex):
            return [], "missing_sex_chr"
        return [], "none"

    # sex_chr is populated from here on.
    if sex == "unknown":
        return None, "inferred"
    if is_heterogametic(sex):
        return list(sex_chr), "declared"
    # Homogametic: sex contigs exist but none are haploid for this sample.
    return [], "declared"
