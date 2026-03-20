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


def run_cap3() -> None:
    """Genome assembly via De Bruijn graphs."""
    print("\n--- [CAP 3] GENOME ASSEMBLY ---")
    print("  1) Standard assembly (k-mers)")
    print("  2) Paired-end assembly")
    sub = input("Choice: ").strip()

    if sub == "1":
        print("Enter k-mers (empty line to stop):")
        kmers: list[str] = []
        while True:
            line = input().strip().upper()
            if not line:
                break
            kmers.append(line)
        if not kmers:
            raise ValueError("No k-mers provided.")
        graph = bio.build_de_bruijn_from_kmers(kmers)
        path = bio.find_eulerian_path(graph)
        print(f"Assembled genome: {bio.path_to_genome(path)}")

    elif sub == "2":
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
            raise ValueError("No pairs provided.")
        graph = bio.build_paired_de_bruijn(pairs)
        path = bio.find_eulerian_path(graph)
        print(f"Assembled genome: {bio.paired_path_to_genome(path, k, d)}")

    else:
        print("Invalid choice.")


def run_cap5() -> None:
    """Sequence alignment: global (Needleman-Wunsch) or local (Smith-Waterman)."""
    print("\n--- [CAP 5] SEQUENCE ALIGNMENT ---")
    s1 = input("Sequence 1: ").strip().upper()
    s2 = input("Sequence 2: ").strip().upper()
    if not s1 or not s2:
        print("Empty sequence.")
        return

    print("  1) Global alignment (Needleman-Wunsch)")
    print("  2) Local alignment  (Smith-Waterman)")
    mode = input("Choice: ").strip()
    if mode not in ("1", "2"):
        print("Invalid choice.")
        return

    print("\n  Substitution matrix:")
    matrix_names = list(bio.SUBSTITUTION_MATRICES.keys())
    for i, name in enumerate(matrix_names, 1):
        print(f"  {i}) {name}")
    print(f"  {len(matrix_names)+1}) None (flat match/mismatch scores)")
    mat_choice = input("Choice: ").strip()

    matrix = None
    matrix_label = "IDENTITY (flat)"
    gap = -2
    try:
        idx = int(mat_choice) - 1
        if 0 <= idx < len(matrix_names):
            matrix_label = matrix_names[idx]
            matrix = bio.SUBSTITUTION_MATRICES[matrix_label]
            gap = -4   # more sensible default for PAM/BLOSUM
    except ValueError:
        pass

    if matrix is not None:
        raw_gap = input(f"Gap penalty (default {gap}): ").strip()
        if raw_gap:
            gap = int(raw_gap)

    fn = bio.global_alignment if mode == "1" else bio.local_alignment
    label = "GLOBAL" if mode == "1" else "LOCAL"
    a1, a2, score = fn(s1, s2, matrix=matrix, gap=gap)

    mid = "".join(
        "|" if a1[k] == a2[k] and a1[k] != "-" else " "
        for k in range(len(a1))
    )
    print(f"\nOptimal {label} alignment  [matrix={matrix_label}  gap={gap}]  score={score}")
    print(f"  S1: {a1}")
    print(f"      {mid}")
    print(f"  S2: {a2}")


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

def run_cap6() -> None:
    """Genome rearrangements: greedy sorting by reversals + breakpoint count."""
    print("\n--- [CAP 6] GENOME REARRANGEMENTS ---")
    print("  1) Greedy sorting by reversals")
    print("  2) Count breakpoints")
    sub = input("Choice: ").strip()
    if sub not in ("1", "2"):
        print("Invalid choice.")
        return

    print("Enter signed permutation (e.g.  +1 -3 +2  or  1 -3 2):")
    raw = input().strip()
    if not raw:
        print("Empty input.")
        return

    # Parse tokens: accept +n, -n, n
    tokens = raw.replace(",", " ").split()
    try:
        perm = [int(t) for t in tokens]
    except ValueError as err:
        raise ValueError(f"Could not parse permutation: {err}") from err

    if sub == "1":
        steps = bio.greedy_sorting(perm)
        n_reversals = len(steps)
        print(f"\nApproximate reversal distance: {n_reversals}")
        if steps:
            show = input("Show intermediate steps? [y/N]: ").strip().lower()
            if show == "y":
                print(f"  {'(':>3}  {_fmt_perm(perm)}  ← initial")
                for k, step in enumerate(steps, 1):
                    print(f"  {k:>3}  {_fmt_perm(step)}")

    elif sub == "2":
        bp = bio.count_breakpoints(perm)
        n = len(perm)
        print(f"\nBreakpoints: {bp}  (out of {n + 1} possible adjacent pairs)")
        print(f"Lower bound on reversal distance ≥ {bp // 2}")


def _fmt_perm(p: list[int]) -> str:
    """Format a signed permutation as  (+1 -3 +2)."""
    return "(" + " ".join(f"+{x}" if x > 0 else str(x) for x in p) + ")"


MENU = {
    "1": ("Origin Search      (Cap 1)",   run_cap1),
    "2": ("Motif Search       (Cap 2)",   run_cap2),
    "3": ("Genome Assembly    (Cap 3)",   run_cap3),
    "4": ("Protein Tools      (Cap 4)",   run_cap4),
    "5": ("Sequence Alignment (Cap 5)",   run_cap5),
    "6": ("Genome Rearrangements (Cap 6)", run_cap6),
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
