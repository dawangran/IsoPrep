"""Utility helpers for filesystem and subprocess operations."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Optional

from .logging import setup_logger

logger = setup_logger(__name__)


def log(msg: str) -> None:
    """Write a convenience INFO log message."""
    logger.info(msg)


def run_cmd(cmd: str, cwd: Optional[Path] = None, check: bool = True) -> None:
    """Execute a shell command and raise RuntimeError on failure when enabled."""
    logger.info("$ %s", cmd)
    proc = subprocess.run(cmd, shell=True, cwd=cwd)
    if check and proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {cmd}")


def safe_mkdir(path: Path) -> None:
    """Create a directory recursively if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)


def sample_name_from_fastq(fq: Path) -> str:
    """Infer sample name by removing FASTQ-related extensions."""
    return re.sub(r"\.fastq\.gz$|\.fq\.gz$|\.fastq$|\.fq$", "", fq.name, flags=re.I)


def symlink_force(src: Path, dst: Path) -> None:
    """Create or replace a symlink from ``dst`` to ``src``."""
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())
