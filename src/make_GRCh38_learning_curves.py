#!/usr/bin/env python3
"""Build reproducible, strictly nested FASTA learning-curve datasets.

For each requested outer set, the positive and negative full training pools are
shuffled once and every percentage is written as a prefix of that fixed order.
Validation FASTAs are deliberately ignored, and no negatives are generated.
"""

from __future__ import annotations

import argparse
import os
import random
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence


PERCENTAGES_DESCENDING = (100, 80, 60, 40, 20, 10)
PERCENTAGES_ASCENDING = tuple(reversed(PERCENTAGES_DESCENDING))
DEFAULT_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = (
    DEFAULT_REPOSITORY_ROOT / "data/processed/evn/GRCh38/training"
)
DEFAULT_OUTPUT_DIR = (
    DEFAULT_REPOSITORY_ROOT / "data/processed/evn/GRCh38/learning_curves"
)


@dataclass(frozen=True)
class FastaRecord:
    """One complete FASTA record preserved exactly as input text."""

    identifier: str
    lines: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create nested 10/20/40/60/80/100 percent training FASTAs "
            "for GRCh38 outer sets."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=(
            "Directory containing outer1 through outer4, each with "
            "train_pos.fa and train_neg.fa."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination directory for outer-specific learning-curve FASTAs.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Base random seed (default: 42).",
    )
    parser.add_argument(
        "--outers",
        type=int,
        nargs="+",
        choices=(1, 2, 3, 4),
        default=(1, 2, 3, 4),
        help="Outer sets to process in one run (default: 1 2 3 4).",
    )
    return parser.parse_args()


def fasta_identifier(header_line: str, path: Path, line_number: int) -> str:
    """Return the FASTA ID (the first token following '>')."""
    identifier = header_line[1:].strip().split(maxsplit=1)[0]
    if not identifier:
        raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
    return identifier


def read_fasta_records(path: Path) -> list[FastaRecord]:
    """Read complete multiline FASTA records and reject duplicate IDs."""
    if not path.is_file():
        raise FileNotFoundError(f"Expected FASTA does not exist: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Expected FASTA is empty: {path}")

    records: list[FastaRecord] = []
    seen_ids: set[str] = set()
    current_lines: list[str] | None = None
    current_id: str | None = None

    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if current_lines is not None:
                    if len(current_lines) == 1:
                        raise ValueError(
                            f"{path}: FASTA record {current_id!r} has no sequence"
                        )
                    records.append(
                        FastaRecord(current_id, tuple(current_lines))  # type: ignore[arg-type]
                    )

                current_id = fasta_identifier(line, path, line_number)
                if current_id in seen_ids:
                    raise ValueError(
                        f"{path}:{line_number}: duplicate FASTA ID {current_id!r}"
                    )
                seen_ids.add(current_id)
                current_lines = [line]
            else:
                if current_lines is None:
                    if line.strip():
                        raise ValueError(
                            f"{path}:{line_number}: sequence text before first header"
                        )
                    continue
                current_lines.append(line)

    if current_lines is not None:
        if len(current_lines) == 1:
            raise ValueError(
                f"{path}: FASTA record {current_id!r} has no sequence"
            )
        records.append(
            FastaRecord(current_id, tuple(current_lines))  # type: ignore[arg-type]
        )

    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def subset_sizes(record_count: int) -> dict[int, int]:
    """Return safe prefix sizes using floor division.

    For percentage p, the size is floor(record_count * p / 100). This is
    deterministic, monotonic, and can never exceed the available record count.
    """
    sizes = {
        percentage: (record_count * percentage) // 100
        for percentage in PERCENTAGES_DESCENDING
    }
    sizes[100] = record_count

    ascending = [sizes[p] for p in PERCENTAGES_ASCENDING]
    if ascending != sorted(ascending):
        raise AssertionError(f"Non-monotonic subset sizes: {sizes}")
    if any(size < 0 or size > record_count for size in sizes.values()):
        raise AssertionError(f"Invalid subset size for N={record_count}: {sizes}")
    return sizes


def class_seed(base_seed: int, outer: int, class_name: str) -> int:
    """Derive independent integer RNG streams for each outer and class."""
    class_offset = 0 if class_name == "pos" else 1
    return base_seed + (outer * 2) + class_offset


def write_records_atomic(path: Path, records: Sequence[FastaRecord]) -> None:
    """Write complete FASTA records atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            for record in records:
                for line in record.lines:
                    handle.write(line)
                # Prevent the next header from joining a final unterminated line.
                if record.lines and not record.lines[-1].endswith(("\n", "\r")):
                    handle.write("\n")
        os.replace(temporary_name, path)
    except BaseException:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
        raise


def write_nested_subsets(
    records: list[FastaRecord],
    output_dir: Path,
    outer: int,
    class_name: str,
    base_seed: int,
) -> tuple[dict[int, Path], dict[int, int]]:
    """Shuffle once, then write every percentage as a prefix."""
    rng = random.Random(class_seed(base_seed, outer, class_name))
    rng.shuffle(records)
    sizes = subset_sizes(len(records))

    paths: dict[int, Path] = {}
    for percentage in PERCENTAGES_DESCENDING:
        output = output_dir / f"train_{percentage}_{class_name}.fa"
        write_records_atomic(output, records[: sizes[percentage]])
        paths[percentage] = output
    return paths, sizes


def iter_fasta_ids(path: Path) -> Iterator[str]:
    """Stream FASTA identifiers without loading sequences."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                yield fasta_identifier(line, path, line_number)


def count_ids(path: Path) -> int:
    return sum(1 for _ in iter_fasta_ids(path))


def assert_id_prefix(smaller: Path, larger: Path) -> None:
    """Verify that smaller IDs exactly match a prefix of larger IDs.

    Exact prefix equality is stronger than set containment and directly checks
    the strategy used to guarantee nested learning-curve subsets.
    """
    larger_ids = iter_fasta_ids(larger)
    for position, smaller_id in enumerate(iter_fasta_ids(smaller), start=1):
        try:
            larger_id = next(larger_ids)
        except StopIteration as error:
            raise AssertionError(
                f"{smaller} has more records than {larger}"
            ) from error
        if smaller_id != larger_id:
            raise AssertionError(
                f"Nested-ID validation failed at record {position}: "
                f"{smaller_id!r} != {larger_id!r}"
            )


def validate_written_subsets(
    paths: dict[int, Path],
    expected_sizes: dict[int, int],
) -> None:
    """Validate counts, monotonicity, and adjacent nested ID containment."""
    observed_sizes = {
        percentage: count_ids(paths[percentage])
        for percentage in PERCENTAGES_ASCENDING
    }

    for percentage, expected in expected_sizes.items():
        observed = observed_sizes[percentage]
        if observed != expected:
            raise AssertionError(
                f"{paths[percentage]}: expected {expected} records, "
                f"observed {observed}"
            )

    counts = [observed_sizes[p] for p in PERCENTAGES_ASCENDING]
    if counts != sorted(counts):
        raise AssertionError(f"Record counts are not monotonic: {observed_sizes}")

    for smaller, larger in zip(
        PERCENTAGES_ASCENDING, PERCENTAGES_ASCENDING[1:]
    ):
        assert_id_prefix(paths[smaller], paths[larger])


def print_class_summary(class_label: str, sizes: dict[int, int]) -> None:
    print(f"{class_label}:")
    for percentage in PERCENTAGES_DESCENDING:
        print(f"{percentage:>3}%: {sizes[percentage]:>10,}")
    print()


def validate_requested_inputs(input_dir: Path, outers: Sequence[int]) -> None:
    """Fail before writing anything if an expected outer training FASTA is absent."""
    missing = []
    for outer in outers:
        for class_name in ("pos", "neg"):
            path = input_dir / f"outer{outer}" / f"train_{class_name}.fa"
            if not path.is_file():
                missing.append(path)
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Missing required training FASTAs:\n{formatted}")


def process_outer(
    outer: int,
    input_dir: Path,
    output_dir: Path,
    base_seed: int,
) -> None:
    input_outer = input_dir / f"outer{outer}"
    output_outer = output_dir / f"outer{outer}"

    print(f"outer{outer}")
    for class_name, label in (("pos", "positive"), ("neg", "negative")):
        input_fasta = input_outer / f"train_{class_name}.fa"
        records = read_fasta_records(input_fasta)
        paths, sizes = write_nested_subsets(
            records=records,
            output_dir=output_outer,
            outer=outer,
            class_name=class_name,
            base_seed=base_seed,
        )
        validate_written_subsets(paths, sizes)
        print_class_summary(label, sizes)

    print("Nested subset validation: PASS")
    print()


def main() -> None:
    args = parse_args()
    outers = tuple(dict.fromkeys(args.outers))
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if input_dir == output_dir:
        raise ValueError("Input and output directories must be different")

    validate_requested_inputs(input_dir, outers)
    for outer in outers:
        process_outer(
            outer=outer,
            input_dir=input_dir,
            output_dir=output_dir,
            base_seed=args.seed,
        )


if __name__ == "__main__":
    main()
