
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sclrtoolkit.logging — Unified logging helper."""
import logging, os
from logging import Logger
from typing import Optional

def setup_logger(name: Optional[str] = None) -> Logger:
    level_name = os.getenv("SCLR_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    show_time = os.getenv("SCLR_LOG_TIME", "1") == "1"
    fmt = "[%(asctime)s] %(levelname)s: %(message)s" if show_time else "%(levelname)s: %(message)s"
    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(level=level, format=fmt)
    logger = logging.getLogger(name if name else __name__)
    logger.setLevel(level)
    return logger
