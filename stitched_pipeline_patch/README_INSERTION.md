# Additional entry point: begin from Imaris-stitched TIFFs

For data already stitched independently in Imaris, use the optional `run_stitched_pipeline.py` controller in this repository. It preserves the eight-step order starting at **initial segmentation** (Step 2) and skips `.ims` export, BaSiC and BigStitcher (Step 1). All existing scientific scripts remain unchanged.

See [STITCHED_WORKFLOW.md](STITCHED_WORKFLOW.md) for execution, calibration, data validation, output layout, provenance, error recovery, and scientific caveats, and [SCRIPT_HANDOFF_REFERENCE.md](SCRIPT_HANDOFF_REFERENCE.md) for a per-script input/output inventory. Edit `stitched_pipeline.example.json` into a project-specific config; run `python run_stitched_pipeline.py --config stitched_pipeline.json --dry-run` to inspect its configured plan, then remove `--dry-run` to execute. Do not use historical pixel-size and filter/overlap values without checking your own data.
