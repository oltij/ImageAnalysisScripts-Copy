#!/usr/bin/env python3

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import tifffile


def write_cellprofiler_binary_mask(source, destination):
    """Write a nonempty 2-D mask as an 8-bit 0/255 CellProfiler input."""
    source = Path(source)
    destination = Path(destination)
    try:
        original = np.squeeze(tifffile.imread(source))
    except Exception as error:
        raise ValueError(f"Could not read organoid mask TIFF:\n{source}") from error
    if original.ndim != 2:
        raise ValueError(
            f"Organoid mask must be a 2-D TIFF; got shape {original.shape}:\n{source}"
        )
    if np.issubdtype(original.dtype, np.floating) and not np.all(
        np.isfinite(original)
    ):
        raise ValueError(f"Organoid mask contains non-finite values:\n{source}")

    binary = original > 0
    foreground_pixels = int(np.count_nonzero(binary))
    if foreground_pixels == 0:
        raise ValueError(f"Organoid mask contains no foreground pixels:\n{source}")

    encoded = np.where(binary, 255, 0).astype(np.uint8)
    tifffile.imwrite(
        destination,
        encoded,
        photometric="minisblack",
        metadata={"axes": "YX"},
    )

    # Verify the exact bytes CellProfiler will read, not only the source array.
    written = np.squeeze(tifffile.imread(destination))
    values = set(int(value) for value in np.unique(written))
    if (
        written.shape != encoded.shape
        or written.dtype != np.uint8
        or not values.issubset({0, 255})
        or not np.array_equal(written, encoded)
    ):
        raise RuntimeError(
            f"Temporary CellProfiler mask failed 8-bit binary validation:\n{destination}"
        )
    return foreground_pixels


def run_cellprofiler(pipeline_path, input_file, output_dir, organoid_mask=None):

    pipeline_path = Path(pipeline_path).resolve()
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    organoid_mask = Path(organoid_mask).resolve() if organoid_mask else None

    # ============================================================
    # CHECK INPUTS
    # ============================================================

    if not pipeline_path.exists():
        raise FileNotFoundError(
            f"CellProfiler pipeline not found:\n{pipeline_path}"
        )

    if not input_file.exists():
        raise FileNotFoundError(
            f"Input image not found:\n{input_file}"
        )

    if organoid_mask is not None and not organoid_mask.exists():
        raise FileNotFoundError(
            f"Organoid mask not found:\n{organoid_mask}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("CELLPROFILER HEADLESS RUN")
    print("=" * 70)
    print(f"Pipeline:   {pipeline_path}")
    print(f"Input file: {input_file}")
    if organoid_mask is not None:
        print(f"Organoid mask: {organoid_mask}")
    print(f"Output dir: {output_dir}")
    print("=" * 70)

    # ============================================================
    # TEMPORARY INPUT DIRECTORY
    # ============================================================

    with tempfile.TemporaryDirectory(
        prefix="cellprofiler_input_"
    ) as temp_dir_string:

        temp_dir = Path(temp_dir_string)

        temporary_input = temp_dir / (
            "__analysis_fluorescence__.tif"
            if organoid_mask is not None
            else input_file.name
        )

        # Symlink instead of copying the large microscopy TIFF
        temporary_input.symlink_to(input_file)
        temporary_mask = None
        if organoid_mask is not None:
            temporary_mask = temp_dir / "__hoechst_organoid_mask__.tif"
            foreground_pixels = write_cellprofiler_binary_mask(
                organoid_mask, temporary_mask
            )

        print("\nTemporary input directory:")
        print(f"  {temp_dir}")

        print("\nLinked input image:")
        print(f"  {temporary_input}")
        if temporary_mask is not None:
            print("\nPrepared Hoechst organoid mask for CellProfiler:")
            print(f"  {temporary_mask}")
            print("  encoding: uint8, background=0, foreground=255")
            print(f"  foreground pixels: {foreground_pixels}")

        # ========================================================
        # RUN CELLPROFILER
        # ========================================================

        command = [
            "cellprofiler",
            "-c",
            "-r",
            "-p", str(pipeline_path),
            "-i", str(temp_dir),
            "-o", str(output_dir),
            "-L", "INFO",
        ]

        print("\nRunning:")
        print(
            " ".join(
                f'"{item}"' if " " in item else item
                for item in command
            )
        )

        print()

        result = subprocess.run(command)

        if result.returncode != 0:
            raise RuntimeError(
                f"CellProfiler exited with code "
                f"{result.returncode}"
            )

        # ========================================================
        # COLLECT FILES WRITTEN BESIDE THE INPUT IMAGE
        # ========================================================
        #
        # Some CellProfiler pipelines have SaveImages and/or
        # ExportToSpreadsheet configured to save into the input
        # image's directory.
        #
        # Since our input directory is temporary, copy all outputs
        # out of it before TemporaryDirectory deletes the directory.
        #
        # The original input symlink itself is skipped.
        # ========================================================

        print()
        print("=" * 70)
        print("COLLECTING CELLPROFILER OUTPUTS")
        print("=" * 70)

        copied_files = []

        for item in temp_dir.iterdir():

            # Do NOT copy the input-image symlink
            if item == temporary_input or item == temporary_mask:
                continue

            destination = output_dir / item.name

            # ----------------------------------------------------
            # Directory output
            # ----------------------------------------------------
            if item.is_dir():

                if destination.exists():
                    shutil.rmtree(destination)

                shutil.copytree(
                    item,
                    destination
                )

                copied_files.append(destination)

                print(
                    f"Copied directory:\n"
                    f"  {item}\n"
                    f"      -> {destination}"
                )

            # ----------------------------------------------------
            # File output
            # ----------------------------------------------------
            elif item.is_file():

                shutil.copy2(
                    item,
                    destination
                )

                copied_files.append(destination)

                print(
                    f"Copied:\n"
                    f"  {item}\n"
                    f"      -> {destination}"
                )

        # ========================================================
        # FINAL REPORT
        # ========================================================

        print()
        print("=" * 70)
        print("CELLPROFILER COMPLETED SUCCESSFULLY")
        print("=" * 70)

        print(f"\nFinal output directory:\n{output_dir}")

        print("\nFiles now present:")

        final_files = sorted(output_dir.iterdir())

        if not final_files:
            print("  [NO FILES FOUND]")

        else:
            for file in final_files:
                print(f"  {file.name}")

        print()

        if copied_files:
            print(
                f"Recovered {len(copied_files)} "
                f"output item(s) written beside the input image."
            )


# ================================================================
# COMMAND-LINE INTERFACE
# ================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Run CellProfiler headlessly on one image and collect "
            "all generated outputs into a specified directory."
        )
    )

    parser.add_argument(
        "--pipeline",
        required=True,
        help="Path to CellProfiler .cppipe pipeline.",
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to input TIFF/image.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Destination directory for all CellProfiler outputs.",
    )

    parser.add_argument(
        "--organoid-mask",
        help=(
            "Optional aligned Hoechst OrganoidMask.tiff. When supplied, the "
            "driver supplies a temporary uint8 0/255 binary copy under the "
            "name __hoechst_organoid_mask__; the pipeline must assign it as "
            "HoechstOrganoidMask. The source mask is not modified."
        ),
    )

    args = parser.parse_args()

    run_cellprofiler(
        pipeline_path=args.pipeline,
        input_file=args.input,
        output_dir=args.output,
        organoid_mask=args.organoid_mask,
    )
