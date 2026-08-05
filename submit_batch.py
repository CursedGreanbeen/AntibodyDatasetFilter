#!/usr/bin/env python3
"""
Submit a batch of FASTA files for AlphaFold folding.
Runs tasks in parallel and saves task IDs to a file.
"""

import argparse
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Optional


# Conda environment activation and script path
CONDA_ACTIVATE = "source ~/miniconda3/etc/profile.d/conda.sh && conda activate folding_env &&"
ALPHAFOLD_SCRIPT = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/scripts/run_alphafold.py")
MAX_WORKERS = 7  # Maximum parallel submissions


def run_command(cmd: str, capture: bool = True) -> tuple[str, int]:
    """Run a shell command and return (stdout, return_code)."""
    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=capture,
        text=True,
        executable="/bin/bash"
    )
    return result.stdout.strip(), result.returncode


def parse_task_id(output: str) -> Optional[str]:
    """Extract TASK_ID from script output."""
    import re
    match = re.search(r'TASK_ID:\s*(\S+)', output)
    return match.group(1) if match else None


def submit_single(fasta_path: Path, memory_gb: int = 128) -> tuple[Path, Optional[str], Optional[str]]:
    """
    Submit a single FASTA file for folding.
    Returns (fasta_path, task_id, error_message).
    """
    stdout, returncode = run_command(
        f"{CONDA_ACTIVATE} python3 {ALPHAFOLD_SCRIPT} -i {fasta_path} --memory={memory_gb}",
        capture=True
    )

    if returncode != 0:
        # Also capture stderr for debugging
        full_cmd = f"{CONDA_ACTIVATE} python3 {ALPHAFOLD_SCRIPT} -i {fasta_path}"
        result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash")
        stderr_output = result.stderr[:500] if result.stderr else "No stderr"
        return fasta_path, None, f"Exit code {returncode}: {stderr_output}"

    task_id = parse_task_id(stdout)
    if not task_id:
        return fasta_path, None, f"Could not parse TASK_ID. Output: {stdout[:200]}"

    return fasta_path, task_id, None


def main():
    parser = argparse.ArgumentParser(description="Submit batch of FASTA files for AlphaFold")
    parser.add_argument(
        "input",
        type=Path,
        help="Directory containing FASTA files or a single FASTA file"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output file for task IDs (default: batch_tasks_YYYYMMDD_HHMMSS.txt)"
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=MAX_WORKERS,
        help=f"Maximum parallel submissions (default: {MAX_WORKERS})"
    )
    parser.add_argument(
        "--memory",
        type=int,
        default=128,
        help="Memory in GB for each AlphaFold task (default: 128)"
    )

    args = parser.parse_args()

    input_path = args.input

    # Collect FASTA files
    if input_path.is_file():
        if input_path.suffix != '.fasta':
            print(f"ERROR: {input_path} is not a FASTA file")
            sys.exit(1)
        fasta_files = [input_path]
    elif input_path.is_dir():
        fasta_files = sorted(input_path.glob("*.fasta"))
        if not fasta_files:
            print(f"ERROR: No FASTA files found in {input_path}")
            sys.exit(1)
    else:
        print(f"ERROR: {input_path} does not exist")
        sys.exit(1)

    print(f"Found {len(fasta_files)} FASTA file(s)")

    # Determine output file
    if args.output:
        output_file = args.output
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = Path(f"batch_tasks_{timestamp}.txt")

    # Submit all tasks
    print(f"\nSubmitting tasks (max {args.max_workers} parallel)...")
    results: list[tuple[Path, Optional[str], Optional[str]]] = []

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {executor.submit(submit_single, f, args.memory): f for f in fasta_files}

        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            fasta_path, task_id, error = result

            if task_id:
                print(f"  OK: {fasta_path.name} -> {task_id[:12]}...")
            else:
                print(f"  FAILED: {fasta_path.name} -> {error}")

    # Summary
    success_count = sum(1 for _, tid, _ in results if tid)
    fail_count = len(results) - success_count

    print(f"\n{'='*60}")
    print(f"Summary: {success_count} submitted, {fail_count} failed")
    print(f"{'='*60}")

    if success_count == 0:
        print("No tasks were submitted successfully. Exiting.")
        sys.exit(1)

    # Save task IDs
    with open(output_file, 'w') as f:
        for fasta_path, task_id, _ in results:
            if task_id:
                f.write(f"{task_id} {fasta_path}\n")

    print(f"Task IDs saved to: {output_file}")
    print(f"\nTo download results, run:")
    print(f"  python3 download_results.py {output_file}")


if __name__ == "__main__":
    main()
