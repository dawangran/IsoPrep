
from dataclasses import dataclass

@dataclass
class ToolPaths:
    cutadapt: str = "cutadapt"
    tsoclip : str = "tsoclip"
    seqkit: str = "seqkit"
    minimap2: str = "minimap2"
    samtools: str = "samtools"
    pisa: str = "PISA"
    python: str = "python"

@dataclass
class AdapterModel:
    model: str = ("CGACATGGCTACGATCCGACTTTCTGCGNNNNNNNNNNCCTTCCNNNNNNNNNNCGATGNNNNNNNNNNTTTTTTTTTTT;max_error_rate=0.2;min_overlap=70")
    tso: str = "CCCCTCTGCGTTGATACCACTGCTT"
    # R1 slicing windows
    r1_slice_1: str = "1:80"
    r1_slice_2: str = "70:-1"
    # add_cb_umi_db model (B=CB, U=UMI) — adjust to your design
    cb_umi_model: str = "BBBBBBBBBBBBBBBBBBBBUUUUUUUUUU"

@dataclass
class ShardedParams:
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
    min_read1_len: int = 300
    min_retain_len: int = 100
