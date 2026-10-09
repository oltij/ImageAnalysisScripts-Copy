"""Crash-safe TIFF writes and structural completeness checks.

Pipeline TIFFs are written to a temporary file in the destination directory,
validated, flushed, and atomically renamed into place.  A process interruption
therefore leaves either the previous complete file or no final file, never a
partially written final pathname.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import tifffile


def validate_tiff(
    path: Path,
    expected_shape: tuple[int, ...] | None = None,
    expected_dtype: Any | None = None,
) -> tuple[tuple[int, ...], np.dtype[Any]]:
    """Validate TIFF structure and return its first-series shape and dtype.

    The strip/tile bounds check catches a common interrupted-write case where
    TIFF metadata exists but the pixel payload ends before its declared byte
    range.  Expected shape and dtype checks prevent reuse of a complete but
    incompatible intermediate.
    """
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"TIFF does not exist: {path}")
    file_size = path.stat().st_size
    if file_size <= 0:
        raise ValueError(f"TIFF is empty: {path}")

    try:
        with tifffile.TiffFile(path) as tif:
            if not tif.series:
                raise ValueError("contains no image series")
            series = tif.series[0]
            shape = tuple(int(value) for value in series.shape)
            dtype = np.dtype(series.dtype)
            if not shape or any(value <= 0 for value in shape):
                raise ValueError(f"has invalid image shape {shape}")
            if expected_shape is not None and shape != tuple(expected_shape):
                raise ValueError(
                    f"has shape {shape}, expected {tuple(expected_shape)}"
                )
            if expected_dtype is not None and dtype != np.dtype(expected_dtype):
                raise ValueError(
                    f"has dtype {dtype}, expected {np.dtype(expected_dtype)}"
                )
            if not tif.pages:
                raise ValueError("contains no TIFF pages")
            for page_number, page in enumerate(tif.pages):
                offsets = tuple(int(value) for value in page.dataoffsets)
                bytecounts = tuple(int(value) for value in page.databytecounts)
                if not offsets or len(offsets) != len(bytecounts):
                    raise ValueError(
                        f"page {page_number} has invalid strip/tile tables"
                    )
                for offset, bytecount in zip(offsets, bytecounts):
                    if offset < 0 or bytecount <= 0 or offset + bytecount > file_size:
                        raise ValueError(
                            f"page {page_number} pixel data extends beyond end of file"
                        )
    except (OSError, tifffile.TiffFileError) as error:
        raise ValueError(f"TIFF cannot be read: {path}: {error}") from error
    return shape, dtype


def is_complete_tiff(
    path: Path,
    expected_shape: tuple[int, ...] | None = None,
    expected_dtype: Any | None = None,
) -> bool:
    """Return whether a TIFF is structurally complete and compatible."""
    try:
        validate_tiff(path, expected_shape, expected_dtype)
    except (OSError, ValueError):
        return False
    return True


def _temporary_tiff(destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    prefix = f".{destination.name}."
    suffix = ".tmp.tif"
    # SIGKILL cannot run a finally block. Remove abandoned private temporary
    # files from earlier attempts before allocating the next one.
    for candidate in destination.parent.iterdir():
        if (
            candidate.name.startswith(prefix)
            and candidate.name.endswith(suffix)
            and (candidate.is_file() or candidate.is_symlink())
        ):
            candidate.unlink()
    descriptor, name = tempfile.mkstemp(
        prefix=prefix,
        suffix=suffix,
        dir=destination.parent,
    )
    os.close(descriptor)
    return Path(name)


def _flush_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def atomic_tiff_write(destination: Path, data: Any, **kwargs: Any) -> None:
    """Write and validate a TIFF before atomically replacing destination."""
    destination = Path(destination)
    temporary = _temporary_tiff(destination)
    array = np.asanyarray(data)
    try:
        tifffile.imwrite(temporary, data, **kwargs)
        validate_tiff(temporary, tuple(array.shape), array.dtype)
        _flush_file(temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_tiff_copy(source: Path, destination: Path) -> None:
    """Validate and atomically publish a copy of an externally written TIFF."""
    source = Path(source)
    destination = Path(destination)
    shape, dtype = validate_tiff(source)
    temporary = _temporary_tiff(destination)
    try:
        shutil.copy2(source, temporary)
        validate_tiff(temporary, shape, dtype)
        _flush_file(temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
