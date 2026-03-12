#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Centralized logging setup for IsoPrep."""

from __future__ import annotations

import logging
import os
from logging import Logger
from typing import Optional


def setup_logger(name: Optional[str] = None) -> Logger:
    """Create a logger configured through environment variables.

    Environment variables:
    - ``SCLR_LOG_LEVEL``: logging level, default ``INFO``.
    - ``SCLR_LOG_TIME``: ``1`` to include timestamps, ``0`` to omit.
    """
    level_name = os.getenv("SCLR_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    show_time = os.getenv("SCLR_LOG_TIME", "1") == "1"
    fmt = "[%(asctime)s] %(levelname)s: %(message)s" if show_time else "%(levelname)s: %(message)s"

    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(level=level, format=fmt)

    logger = logging.getLogger(name or __name__)
    logger.setLevel(level)
    return logger
