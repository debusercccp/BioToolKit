"""
main.py — CLI entry point for BioToolkit.

Each chapter is an independent handler function.
Run:  python main.py
"""

import sys
import bio_logic as bio

try:
    import matplotlib.pyplot as plt
    _HAS_MATPLOTLIB = True
except ImportError:
    _HAS_MATPLOTLIB = False

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_dna(prompt: str = "DNA sequence") -> str:
    raw = input(f"{prompt}: ").strip().upper()
    cleaned = "".join(b for b in raw if b in "ACGT")
    if not cleaned:
        raise ValueError("Empty or invalid DNA sequence.")
    return cleaned


def _read_int(prompt: str, default: int) -> int:
    raw = input(f"{prompt} (default {default}): ").strip()
    return int(raw) if raw else default


def _read_dna_list() -> list[str]:
    print("Enter sequences one per line (empty line to stop):")
    seqs: list[str] = []
    while True:
        line = input().strip().upper()
        if not line:
            break
        cleaned = "".join(b for b in line if b in "ACGT")
        if cleaned:
            seqs.append(cleaned)
    if not seqs:
        raise ValueError("No sequences provided.")
    return seqs


# ---------------------------------------------------------------------------
# Chapter handlers
# ---------------------------------------------------------------------------

def run_cap1() -> None:
    """OriC finder: skew minimum + DnaA box search."""
    print("\n--- [CAP 1] ORIGIN FINDER ---")
    genome = _read_dna()
    skew = bio.compute_skew(genome)
    min_val = min(skew)
    min_positions = [i for i, v in enumerate(skew) if v == min_val]
    # Use the last minimum as the OriC candidate
    ori_pos = min_positions[-1]
    print(f"Skew minimum at position: {ori_pos}  (value={min_val})")

    window = genome[max(0, ori_pos - 250) : min(len(genome), ori_pos + 250)]
    k = _read_int("k-mer length", 9)
    d = _read_int("Max mismatches", 1)

    boxes = bio.frequent_words_with_mismatches_and_rc(window, k, d)
    print(f"Candidate DnaA boxes ({len(boxes)} found): {boxes}")

    if _HAS_MATPLOTLIB:
        plt.figure(figsize=(12, 4))
        plt.plot(skew, linewidth=0.8)
        plt.axvline(x=ori_pos, color='r', linestyle='--', label=f"OriC ≈ {ori_pos}")
        plt.title("GC Skew")
        plt.xlabel("Position")
        plt.ylabel("#G − #C")
        plt.legend()
        plt.tight_layout()
        plt.show()
    else:
        print("(matplotlib not installed — skipping plot)")


def run_cap2() -> None:
    """Randomized motif search."""
    print("\n--- [CAP 2] MOTIF FINDER ---")
    dna_list = _read_dna_list()
    k = _read_int("Motif length k", 15)
    iterations = _read_int("Search iterations", 1000)

    best_motifs = bio.randomized_motif_search(dna_list, k, iterations)
    print(f"\nBest motifs found (score={bio.get_score(best_motifs)}):")
    for m in best_motifs:
        print(f"  {m}")


def run_cap3_standard() -> None:
    """Standard genome assembly from k-mers."""
    print("\n--- [CAP 3.1] STANDARD ASSEMBLY ---")
    print("Enter k-mers (empty line to stop):")
    kmers: list[str] = []
    while True:
        line = input().strip().upper()
        if not line:
            break
        kmers.append(line)
    if not kmers:
        print("No k-mers provided.")
        return

    graph = bio.build_de_bruijn_from_kmers(kmers)
    path = bio.find_eulerian_path(graph)
    result = bio.path_to_genome(path)
    print(f"Assembled genome: {result}")


def run_cap3_paired() -> None:
    """Paired-end genome assembly."""
    print("\n--- [CAP 3.2] PAIRED-END ASSEMBLY ---")
    k = _read_int("k-mer length", 3)
    d = _read_int("Gap distance d", 1)

    print("Enter pairs in format  read1|read2  (empty line to stop):")
    pairs: list[tuple[str, str]] = []
    while True:
        line = input().strip().upper()
        if not line:
            break
        if "|" not in line:
            print("  Skipping (no '|' separator).")
            continue
        r1, r2 = line.split("|", 1)
        pairs.append((r1, r2))
    if not pairs:
        print("No pairs provided.")
        return

    graph = bio.build_paired_de_bruijn(pairs)
    path = bio.find_eulerian_path(graph)
    result = bio.paired_path_to_genome(path, k, d)
    print(f"Assembled genome: {result}")


def run_cap4() -> None:
    """Protein / mass-spectrometry tools."""
    print("\n--- [CAP 4] PROTEIN TOOLS ---")
    print("  1) Peptide → theoretical spectrum")
    print("  2) Spectrum → cyclopeptide sequencing")
    sub = input("Choice: ").strip()

    if sub == "1":
        peptide = input("Amino acid sequence (e.g. LEQN): ").strip().upper()
        try:
            masses = [bio.MASS_TABLE[aa] for aa in peptide]
        except KeyError as err:
            print(f"Unknown amino acid: {err}")
            return

        print("  1) Linear spectrum")
        print("  2) Circular spectrum")
        kind = input("Choice: ").strip()
        spec = (
            bio.get_linear_spectrum(masses)
            if kind == "1"
            else bio.get_circular_spectrum(masses)
        )
        print(f"Spectrum: {' '.join(map(str, spec))}")

    elif sub == "2":
        raw = input("Experimental spectrum (space-separated integers): ").strip()
        try:
            exp_spec = sorted(int(m) for m in raw.split())
        except ValueError:
            print("Error: please enter integers only.")
            return

        print("Searching… (may take a moment for large spectra)")
        results = bio.cyclopeptide_sequencing(exp_spec)
        if not results:
            print("No matching peptides found.")
        else:
            unique = sorted(set("-".join(map(str, p)) for p in results))
            print(f"{len(unique)} peptide(s) found:")
            for p in unique:
                print(f"  {p}")
    else:
        print("Invalid choice.")


# ---------------------------------------------------------------------------
# Menu
# ---------------------------------------------------------------------------

MENU = {
    "1": ("Origin Search     (Cap 1)", run_cap1),
    "2": ("Motif Search      (Cap 2)", run_cap2),
    "3": ("Standard Assembly (Cap 3.1)", run_cap3_standard),
    "4": ("Paired-End Asm    (Cap 3.2)", run_cap3_paired),
    "5": ("Protein Tools     (Cap 4)", run_cap4),
}


def main() -> None:
    while True:
        print("\n" + "=" * 34)
        print("        B I O T O O L K I T")
        print("=" * 34)
        for key, (label, _) in MENU.items():
            print(f"  {key}) {label}")
        print("  q) Quit")

        choice = input("\nChoice: ").strip().lower()
        if choice == "q":
            print("Bye!")
            sys.exit(0)

        entry = MENU.get(choice)
        if entry is None:
            print("Invalid option.")
            continue

        _, handler = entry
        try:
            handler()
        except (ValueError, KeyError) as err:
            print(f"[Error] {err}")
        except KeyboardInterrupt:
            print("\n(interrupted)")


if __name__ == "__main__":
    main()
