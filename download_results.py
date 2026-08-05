#!/usr/bin/env python3
"""
Download results from completed AlphaFold batch tasks.
Reads task IDs from a batch file and downloads results using oneq get-results.

Results are saved to the specified output directory (default: current directory).
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional


# OneQ results are typically saved to:
# - Current working directory if no output specified
# - Or to the directory specified with -o / --output flag
DEFAULT_OUTPUT_DIR = Path.cwd()


def run_command(cmd: str, capture: bool = True) -> tuple[str, int]:
    """Run a shell command and return (stdout, return_code)."""
    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=capture,
        text=True
    )
    if capture:
        return result.stdout.strip(), result.returncode
    else:
        return "", result.returncode


def parse_batch_file(batch_path: Path) -> list[tuple[str, str]]:
    """
    Parse batch tasks file and return list of (task_id, fasta_path) tuples.
    """
    tasks = []
    content = batch_path.read_text()

    for line in content.strip().splitlines():
        parts = line.strip().split()
        if len(parts) >= 2:
            task_id = parts[0]
            fasta_path = parts[1]
            tasks.append((task_id, fasta_path))

    return tasks


def download_task_results(task_id: str, output_dir: Path) -> tuple[str, bool]:
    """
    Download results for a single task using oneq get-results.
    Returns (message, success).
    """
    # Use oneq get-results command: oneq get-results TASK_ID SAVE_PATH
    cmd = f"oneq get-results {task_id} {output_dir}"

    stdout, returncode = run_command(cmd, capture=False)

    if returncode != 0:
        return f"Failed to download task {task_id}: exit code {returncode}", False

    return f"Downloaded results for task {task_id[:12]}... to {output_dir}", True


def check_task_status(task_id: str) -> tuple[str, bool]:
    """
    Check if a task is completed using oneq status.
    Returns (status_message, is_completed).
    """
    cmd = f"oneq status {task_id}"
    stdout, returncode = run_command(cmd, capture=True)

    if returncode != 0:
        return f"Could not get status for {task_id}: {stdout}", False

    # Check if task is completed (oneq typically shows "COMPLETED", "RUNNING", "FAILED", etc.)
    is_completed = "COMPLETED" in stdout.upper() or "FINISHED" in stdout.upper()
    return stdout, is_completed


def main():
    parser = argparse.ArgumentParser(
        description="Download results from completed AlphaFold batch tasks"
    )
    parser.add_argument(
        "batch_file",
        type=Path,
        help="Batch tasks file (e.g., batch_tasks_20260717_163616.txt)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory for results (default: current directory)"
    )
    parser.add_argument(
        "--check-status",
        action="store_true",
        help="Check task status before downloading (skip completed tasks)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force download even if results already exist"
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Start from task N (0-based index, default: 0)"
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="End at task N (exclusive, default: all tasks)"
    )

    args = parser.parse_args()

    # Validate batch file
    if not args.batch_file.exists():
        print(f"ERROR: Batch file not found: {args.batch_file}")
        sys.exit(1)

    # Parse batch file
    tasks = parse_batch_file(args.batch_file)

    if not tasks:
        print(f"ERROR: No tasks found in {args.batch_file}")
        sys.exit(1)

    print(f"Found {len(tasks)} tasks in batch file")

    # Apply slicing if specified
    if args.start > 0:
        tasks = tasks[args.start:]
    if args.end is not None:
        tasks = tasks[:args.end]

    print(f"Processing {len(tasks)} tasks (start={args.start}, end={args.end})\n")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Download results
    results: list[tuple[str, str, bool]] = []  # (task_id, message, success)

    for i, (task_id, fasta_path) in enumerate(tasks):
        print(f"[{i+1}/{len(tasks)}] Task {task_id[:12]}... ({fasta_path})")

        # Check status if requested
        if args.check_status:
            status, is_completed = check_task_status(task_id)
            if not is_completed:
                print(f"  Skipping: Task not completed yet")
                results.append((task_id, "Skipped (not completed)", False))
                continue

        # Download results
        message, success = download_task_results(task_id, args.output_dir)
        print(f"  {message}")
        results.append((task_id, message, success))

    # Summary
    success_count = sum(1 for _, _, s in results if s)
    fail_count = len(results) - success_count

    print(f"\n{'='*60}")
    print(f"Summary: {success_count} downloaded, {fail_count} failed/skipped")
    print(f"{'='*60}")

    # Show failed tasks
    if fail_count > 0:
        print("\nFailed/Skipped tasks:")
        for task_id, message, success in results:
            if not success:
                print(f"  - {task_id}: {message}")

    # Exit with error if any failures
    sys.exit(0 if fail_count == 0 else 1)


if __name__ == "__main__":
    main()
