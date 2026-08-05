from pathlib import Path


CIF_EXTENSIONS = {".cif", ".mmcif"}
FASTA_EXTENSIONS = {".fasta", ".fa", ".faa"}


def find_cif_files(cif_dir: Path) -> list[Path]:
    """Возвращает все CIF/mmCIF-файлы в указанной папке."""
    return sorted(
        path
        for path in cif_dir.iterdir()
        if path.is_file() and path.suffix.lower() in CIF_EXTENSIONS
    )


def find_fasta_files(fasta_dir: Path) -> list[Path]:
    """Возвращает все FASTA-файлы в указанной папке."""
    return sorted(
        path
        for path in fasta_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )


def build_fasta_index(fasta_dir: Path) -> dict[str, Path]:
    """
    Создаёт индекс FASTA-файлов по имени без расширения.

    Например:
    7abc.cif -> ключ 7abc
    7abc.fasta -> ключ 7abc
    """
    index = {}

    for fasta_path in find_fasta_files(fasta_dir):
        key = fasta_path.stem

        if key in index:
            raise ValueError(
                f"Найдено несколько FASTA с одинаковым именем: {key}"
            )

        index[key] = fasta_path

    return index


def match_cif_to_fasta(
    cif_path: Path,
    fasta_index: dict[str, Path],
) -> Path | None:
    """Возвращает FASTA, связанный с CIF по имени файла."""
    return fasta_index.get(cif_path.stem)
