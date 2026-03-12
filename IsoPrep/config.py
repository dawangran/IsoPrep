"""Configuration dataclasses for tools and pipeline defaults."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ToolPaths:
    """Executable paths used by the pipeline."""

    cutadapt: str = "cutadapt"
    tsoclip: str = "tsoclip"
    seqkit: str = "seqkit"
    minimap2: str = "minimap2"
    samtools: str = "samtools"
    pisa: str = "PISA"
    python: str = "python"


@dataclass
class AdapterModel:
    """Adapter/TSO model and read slicing configuration."""

    model: str = (
        "CGACATGGCTACGATCCGACTTTCTGCGNNNNNNNNNNCCTTCCNNNNNNNNNNCGATGNNNNNNNNNNTTTTTTTTTTT;"
        "max_error_rate=0.2;min_overlap=70"
    )
    tso: str = "CCCCTCTGCGTTGATACCACTGCTT"
    r1_slice_1: str = "1:80"
    r1_slice_2: str = "70:-1"
    cb_umi_model: str = "BBBBBBBBBBBBBBBBBBBBUUUUUUUUUU"


@dataclass
class ShardedParams:
    """Legacy sharding-related default parameters."""

    shards: int = 32
    threads: int = 32
    no_gene: bool = True
    fast_h1: bool = True
    ham: int = 1
    ratio: float = 2.0
    jitter: int = 10
    locus_bin: int = 1000


@dataclass
class Defaults:
    """General pipeline threshold defaults."""

    min_read1_len: int = 300
    min_retain_len: int = 100
