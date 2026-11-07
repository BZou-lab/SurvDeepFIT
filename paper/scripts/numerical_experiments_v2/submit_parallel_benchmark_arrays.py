"""Submit arrays for an already snapshotted and submitted setup job."""
import json
from pathlib import Path
import subprocess
import sys

record_path = Path(sys.argv[1]).resolve() / "submission.json"
record = json.loads(record_path.read_text())
run = Path(record["run_root"])
root = Path(record["benchmark_root"])
project = record["project_root"]
jobs = record["jobs"]


def submit(key, options, script, arguments):
    if key in jobs:
        print(f"{key}: already submitted as {jobs[key]}", flush=True)
        return jobs[key]
    command = ["sbatch", "--parsable", "--partition=general", *options,
               f"--output={run}/logs/{key}_%A_%a.out",
               f"--error={run}/logs/{key}_%A_%a.err",
               str(root / "slurm" / script), project, str(root), str(run), *arguments]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    job_id = result.stdout.strip().split(";")[0]
    if not job_id.isdigit():
        raise RuntimeError(f"Unexpected sbatch output: {result.stdout}")
    jobs[key] = job_id
    record.setdefault("submission_commands", {})[key] = command
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    print(f"{key}: {job_id}", flush=True)
    return job_id


arrays = []
for workers in (1, 2, 4):
    arrays.append(submit(f"w{workers}", [f"--job-name=surv-par-w{workers}",
        "--array=1-12%4", f"--cpus-per-task={workers}", "--mem=16G",
        "--time=24:00:00", f"--dependency=afterok:{jobs['setup']}",
        "--kill-on-invalid-dep=yes"], "run_parallel_benchmark.slurm", [str(workers)]))
submit("summary", ["--job-name=surv-par-summary", "--cpus-per-task=1", "--mem=4G",
    "--time=00:30:00", "--dependency=afterany:" + ":".join(arrays)],
    "summarize_parallel_benchmark.slurm", [])
