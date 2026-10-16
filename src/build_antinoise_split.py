#!/usr/bin/env python3
"""Build a split-aware antinoise background_genome_mono executable.

The upstream program opens and samples every chromosome for a recognized genome.
This builder creates a temporary patched copy accepting one comma-separated allowed-
chromosome argument. It opens and samples only those chromosomes, allowing a reference
directory to contain just the permitted ``chrN.plain`` files. The upstream checkout is
never modified.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=repo.parent / "antinoise" / "src")
    parser.add_argument("--output", type=Path, default=repo / "bin" / "background_genome_mono_split.exe")
    parser.add_argument("--compiler", default="g++")
    return parser.parse_args()


def replace_once(text: str, old: str, new: str, label: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"Expected exactly one {label} patch anchor; found {text.count(old)}")
    return text.replace(old, new)


def patched_source(source: str) -> str:
    helper = r"""
int ChromAllowed(const char *chromosome, const char *allowed)
{
    const char *cursor = allowed;
    while (*cursor != '\0')
    {
        while (*cursor == ',' || *cursor == ' ' || *cursor == '\t')cursor++;
        const char *start = cursor;
        while (*cursor != '\0' && *cursor != ',')cursor++;
        const char *name = start;
        if (cursor - start >= 3 && strncmp(start, "chr", 3) == 0)name += 3;
        const int length = (int)(cursor - name);
        if ((int)strlen(chromosome) == length && strncmp(name, chromosome, length) == 0)return 1;
        if (*cursor == ',')cursor++;
    }
    return 0;
}
"""
    source = replace_once(source, "int main(int argc, char *argv[])\n{", helper + "\nint main(int argc, char *argv[])\n{", "helper")
    source = replace_once(source, "if (argc != 15)", "if (argc != 16)", "argc")
    source = replace_once(
        source,
        'puts("8double stop_fraction 9int sequence_length (the same OR (>0 && <3000) 10,11,12,13file_out (mono-, di-) x (forground_vs_background, statistics) 14file_out log");',
        'puts("8double stop_fraction 9int sequence_length 10,11,12,13 statistics 14 log 15 allowed chromosomes (comma-separated, without chr prefix)");',
        "usage",
    )
    source = replace_once(
        source,
        "strcpy(file_log, argv[14]);//out_file sta dinucl content\t",
        "strcpy(file_log, argv[14]);//out_file sta dinucl content\n\tconst char *allowed_chromosomes = argv[15];",
        "allowed argument",
    )
    anchor = "\tif (genome_rec == 0)\n\t{\n\t\tprintf(\"Genome %s is not recognized\\n\",genome);\n\t\texit(1);\n\t}\t\n\tint rand_max_win_linux = 32767;"
    replacement = "\tif (genome_rec == 0)\n\t{\n\t\tprintf(\"Genome %s is not recognized\\n\",genome);\n\t\texit(1);\n\t}\n\tint allowed_count = 0;\n\tfor (i = 0; i < n_chr; i++)\n\t{\n\t\tif (ChromAllowed(name_chr[i], allowed_chromosomes))allowed_count++;\n\t\telse sizelo1[i] = 0;\n\t}\n\tif (allowed_count == 0)\n\t{\n\t\tprintf(\"No recognized chromosomes in allowed set %s\\n\", allowed_chromosomes);\n\t\texit(1);\n\t}\n\tint rand_max_win_linux = 32767;"
    source = replace_once(source, anchor, replacement, "chromosome filter")
    open_anchor = """	for (i = 0; i < n_chr; i++)
	{
		strcpy(filechr[i], path_fasta);
		strcat(filechr[i], filesta);
		strcat(filechr[i], name_chr[i]);
		strcat(filechr[i], fileend);
		if ((in_seq[i] = fopen(filechr[i], "rt")) == NULL)
		{
			printf("Input file %s can't be opened!\\n", filechr[i]);
			exit(1);
		}
	}"""
    open_replacement = """	for (i = 0; i < n_chr; i++)
	{
		in_seq[i] = NULL;
		if (!ChromAllowed(name_chr[i], allowed_chromosomes))continue;
		strcpy(filechr[i], path_fasta);
		strcat(filechr[i], filesta);
		strcat(filechr[i], name_chr[i]);
		strcat(filechr[i], fileend);
		if ((in_seq[i] = fopen(filechr[i], "rt")) == NULL)
		{
			printf("Input file %s can't be opened!\\n", filechr[i]);
			exit(1);
		}
	}"""
    source = replace_once(source, open_anchor, open_replacement, "allowed file opening")
    source = replace_once(
        source,
        "\tfor (i = 0; i < n_chr; i++)fclose(in_seq[i]);",
        "\tfor (i = 0; i < n_chr; i++)if (in_seq[i] != NULL)fclose(in_seq[i]);",
        "allowed file closing",
    )
    # Advance the counted prefix before an early continue, otherwise the next
    # foreground can count that same suffix again and get a false A/T match.
    source = replace_once(source, "double fatc = (double)(cat) / sort[i].len;",
                          "len_cur = sort[i].len;\n\t\t\t\t\tdouble fatc = (double)(cat) / sort[i].len;",
                          "A/T counted-prefix tracking")
    return source


def main() -> None:
    args = parse_args()
    source_dir = args.source_dir.resolve()
    cpp = source_dir / "background_genome_mono.cpp"
    headers = (source_dir / "chromosome_name.h", source_dir / "chromosome_len.h")
    for path in (cpp, *headers):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="antinoise_split_build.") as directory:
        build = Path(directory)
        (build / cpp.name).write_text(patched_source(cpp.read_text()))
        for header in headers:
            shutil.copy2(header, build / header.name)
        subprocess.run([args.compiler, "-O3", "-std=c++11", str(build / cpp.name), "-o", str(args.output)], check=True)
    args.output.chmod(0o755)
    print(f"Built {args.output}")


if __name__ == "__main__":
    main()
