
import re, subprocess
from pathlib import Path
from sclrtoolkit.logging import setup_logger

logger = setup_logger(__name__)

def log(msg: str):
    logger.info(msg)

def run_cmd(cmd: str, cwd: Path = None, check: bool = True):
    logger.info(f"$ {cmd}")
    proc = subprocess.run(cmd, shell=True, cwd=cwd)
    if check and proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {cmd}")

def safe_mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def sample_name_from_fastq(fq: Path) -> str:
    name = fq.name
    name = re.sub(r"\.fastq\.gz$|\.fq\.gz$|\.fastq$|\.fq$", "", name, flags=re.I)
    return name

def symlink_force(src: Path, dst: Path):
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())
