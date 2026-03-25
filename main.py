"""
main.py — CLI entry point for BioToolkit.

Each chapter is an independent handler function.
Run:  python main.py
"""

import sys
import math as _math
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


def _read_distance_matrix() -> tuple[list[list[float]], list[str]]:
    """
    Read a distance matrix from stdin.
    First line: space-separated labels.
    Following n lines: one row of floats each.
    """
    print("Enter labels (space-separated, e.g.  A B C D):")
    raw_labels = input().strip().split()
    if not raw_labels:
        raise ValueError("No labels provided.")
    n = len(raw_labels)
    print(f"Enter {n} rows of {n} distances (one row per line):")
    matrix: list[list[float]] = []
    for i in range(n):
        row_raw = input().strip().split()
        if len(row_raw) != n:
            raise ValueError(f"Row {i+1}: expected {n} values, got {len(row_raw)}.")
        matrix.append([float(x) for x in row_raw])
    for i in range(n):
        for j in range(i + 1, n):
            if abs(matrix[i][j] - matrix[j][i]) > 1e-9:
                raise ValueError(
                    f"Matrix not symmetric at ({i},{j}): "
                    f"{matrix[i][j]} vs {matrix[j][i]}."
                )
    return matrix, raw_labels


def _print_edges(edges: list[tuple[int, int, float]], labels: list[str]) -> None:
    """Pretty-print tree edges with label resolution."""
    n = len(labels)
    def _name(node: int) -> str:
        return labels[node] if node < n else f"n{node}"
    for u, v, w in edges:
        print(f"  {_name(u):>6} -- {_name(v):<6}  {w:.4f}")


def run_cap7() -> None:
    """Molecular evolution: UPGMA, Neighbor Joining, Small Parsimony."""
    print("\n--- [CAP 7] MOLECULAR EVOLUTION ---")
    print("  1) UPGMA")
    print("  2) Neighbor Joining")
    print("  3) Small Parsimony (Fitch)")
    sub = input("Choice: ").strip()

    if sub in ("1", "2"):
        matrix, labels = _read_distance_matrix()
        n = len(labels)
        if sub == "1":
            tree, edges = bio.upgma(matrix, labels)
            print(f"\nUPGMA tree  ({n} leaves, {len(tree) - n} internal nodes):")
        else:
            tree, edges = bio.neighbor_joining(matrix, labels)
            print(f"\nNeighbor Joining tree  ({n} leaves):")
        _print_edges(edges, labels)
        print(f"\nNewick: {bio.newick(tree, labels)}")

    elif sub == "3":
        print("Enter leaf labels (space-separated):")
        raw_labels = input().strip().split()
        if not raw_labels:
            raise ValueError("No labels provided.")
        n = len(raw_labels)
        print(f"Enter aligned sequences for each leaf ({n} sequences):")
        leaf_sequences: dict[int, str] = {}
        for i, label in enumerate(raw_labels):
            seq = input(f"  {label}: ").strip().upper()
            if not seq:
                raise ValueError(f"Empty sequence for leaf '{label}'.")
            leaf_sequences[i] = seq
        seqs = list(leaf_sequences.values())
        seq_len = len(seqs[0])
        if any(len(s) != seq_len for s in seqs):
            raise ValueError("All sequences must have the same length.")
        print("\nBuilding tree via Neighbor Joining on Hamming distances...")
        dist = [
            [bio.hamming_distance(seqs[i], seqs[j]) for j in range(n)]
            for i in range(n)
        ]
        tree, _ = bio.neighbor_joining(dist, raw_labels)
        score, node_seqs = bio.parsimony_score(tree, leaf_sequences)
        print(f"\nParsimony score: {score}")
        print("\nReconstructed ancestral sequences:")
        for node, chars in node_seqs.items():
            if node >= n:
                print(f"  n{node}: {''.join(chars)}")

    else:
        print("Invalid choice.")


def run_cap8() -> None:
    """Gene clustering: hard K-Means, soft K-Means, hierarchical clustering."""
    print("\n--- [CAP 8] GENE CLUSTERING ---")
    print("  1) Hard K-Means (Lloyd)")
    print("  2) Soft K-Means")
    print("  3) Hierarchical clustering")
    sub = input("Choice: ").strip()
    if sub not in ("1", "2", "3"):
        print("Invalid choice.")
        return

    print("Enter data points (one per line, space-separated floats; empty line to stop):")
    data: list[list[float]] = []
    while True:
        line = input().strip()
        if not line:
            break
        try:
            data.append([float(x) for x in line.split()])
        except ValueError:
            print("  Skipping invalid line.")
    if not data:
        raise ValueError("No data points provided.")

    if sub == "1":
        k = _read_int("Number of clusters k", 2)
        centers, clusters = bio.lloyd_algorithm(data, k)
        print(f"\nK-Means converged ({k} clusters):")
        for ci, (center, cluster) in enumerate(zip(centers, clusters)):
            c_fmt = "  ".join(f"{v:.3f}" for v in center)
            print(f"  Cluster {ci+1}: center=[{c_fmt}]  n={len(cluster)}")

    elif sub == "2":
        k    = _read_int("Number of clusters k", 2)
        beta_raw = input("Stiffness beta (default 1.0): ").strip()
        beta = float(beta_raw) if beta_raw else 1.0
        centers, resp = bio.soft_kmeans(data, k, beta=beta)
        print(f"\nSoft K-Means converged (k={k}, beta={beta}):")
        for ci, center in enumerate(centers):
            c_fmt = "  ".join(f"{v:.3f}" for v in center)
            print(f"  Center {ci+1}: [{c_fmt}]")
        show = input("Show responsibilities? [y/N]: ").strip().lower()
        if show == "y":
            print("  Point  " + "  ".join(f"C{j+1}" for j in range(k)))
            for i, row in enumerate(resp):
                r_fmt = "  ".join(f"{r:.3f}" for r in row)
                print(f"  {i+1:>5}  {r_fmt}")

    elif sub == "3":
        print("  Linkage:")
        print("  1) Single   2) Complete   3) Average")
        lk_choice = input("  Choice (default 3): ").strip()
        linkage = {"1": "single", "2": "complete", "3": "average"}.get(lk_choice, "average")
        history = bio.hierarchical_clustering(data, linkage=linkage)
        n = len(data)
        print(f"\nHierarchical clustering ({linkage} linkage)  — merge order:")
        def _lbl(idx: int) -> str:
            return f"p{idx+1}" if idx < n else f"c{idx-n+1}"
        for step, (a, b, dist) in enumerate(history, 1):
            print(f"  step {step:>2}: merge {_lbl(a):>4} + {_lbl(b):<4}  dist={dist:.4f}")


def run_cap9() -> None:
    """Read mapping via BWT: transform, inverse, exact/inexact matching."""
    print("\n--- [CAP 9] READ MAPPING (BWT) ---")
    print("  1) BWT transform")
    print("  2) BWT inverse")
    print("  3) Exact pattern matching (backward search)")
    print("  4) Inexact matching (up to d mismatches)")
    sub = input("Choice: ").strip()

    if sub == "1":
        text = input("Text: ").strip().upper()
        if not text:
            raise ValueError("Empty input.")
        bwt = bio.bwt_transform(text)
        sa  = bio.suffix_array(text)
        print(f"BWT : {bwt}")
        print(f"SA  : {sa}")

    elif sub == "2":
        bwt = input("BWT string: ").strip()
        if not bwt:
            raise ValueError("Empty input.")
        original = bio.bwt_inverse(bwt)
        print(f"Original: {original}")

    elif sub == "3":
        text    = input("Reference text: ").strip().upper()
        pattern = input("Pattern        : ").strip().upper()
        if not text or not pattern:
            raise ValueError("Text and pattern must be non-empty.")
        hits = bio.bwt_match(text, pattern)
        if hits:
            print(f"Found {len(hits)} occurrence(s) at positions: {hits}")
        else:
            print("Pattern not found.")

    elif sub == "4":
        text    = input("Reference text : ").strip().upper()
        pattern = input("Pattern        : ").strip().upper()
        if not text or not pattern:
            raise ValueError("Text and pattern must be non-empty.")
        d = _read_int("Max mismatches d", 1)
        hits = bio.bwt_match_with_mismatches(text, pattern, max_mismatches=d)
        if hits:
            print(f"Found {len(hits)} hit(s):")
            for pos, mm in hits:
                snippet = (text + "$")[pos : pos + len(pattern)]
                print(f"  pos={pos:>5}  mismatches={mm}  [{snippet}]")
        else:
            print("No hits found.")

    else:
        print("Invalid choice.")


def _read_hmm_params() -> tuple[list[str], list[str], dict, dict, dict]:
    """Interactive input for HMM parameters."""
    print("Enter state names (space-separated):")
    states = input().strip().split()
    if not states:
        raise ValueError("No states provided.")
    print("Enter emission alphabet (space-separated):")
    alphabet = input().strip().split()
    if not alphabet:
        raise ValueError("No alphabet provided.")
    print(f"Start probabilities for {states}:")
    start_p: dict[str, float] = {}
    for s in states:
        v = float(input(f"  P(start={s}): ").strip())
        start_p[s] = v
    print("Transition probabilities:")
    trans_p: dict[str, dict[str, float]] = {}
    for s in states:
        trans_p[s] = {}
        for t in states:
            v = float(input(f"  P({s}->{t}): ").strip())
            trans_p[s][t] = v
    print("Emission probabilities:")
    emit_p: dict[str, dict[str, float]] = {}
    for s in states:
        emit_p[s] = {}
        for a in alphabet:
            v = float(input(f"  P({s} emits {a}): ").strip())
            emit_p[s][a] = v
    return states, alphabet, start_p, trans_p, emit_p


def run_cap10() -> None:
    """HMM & gene finding: Viterbi, Forward, Baum-Welch, CpG islands."""
    print("\n--- [CAP 10] HMM & GENE FINDING ---")
    print("  1) Viterbi (most probable path)")
    print("  2) Forward (log P of observation sequence)")
    print("  3) Baum-Welch (parameter training)")
    print("  4) CpG island detection")
    sub = input("Choice: ").strip()

    if sub in ("1", "2", "3"):
        print("\n  HMM:")
        print("  1) Casino (fair/loaded die)  2) Custom")
        hmm_choice = input("  Choice: ").strip()

        if hmm_choice == "1":
            states   = ["F", "L"]
            alphabet = list("123456")
            start_p  = {"F": 0.5, "L": 0.5}
            trans_p  = {"F": {"F": 0.95, "L": 0.05}, "L": {"F": 0.10, "L": 0.90}}
            emit_p   = {
                "F": {str(i): 1/6 for i in range(1, 7)},
                "L": {"1":0.1,"2":0.1,"3":0.1,"4":0.1,"5":0.1,"6":0.5},
            }
            print("Enter observation sequence (e.g. 315662663):")
        else:
            states, alphabet, start_p, trans_p, emit_p = _read_hmm_params()
            print("Enter observation sequence (space-separated symbols):")

        raw = input().strip()
        obs = raw.split() if " " in raw else list(raw)
        if not obs:
            raise ValueError("Empty observation sequence.")

        if sub == "1":
            path, lp = bio.viterbi(obs, states, start_p, trans_p, emit_p)
            print(f"\nMost probable path:  {''.join(path)}")
            print(f"Log probability   :  {lp:.4f}")

        elif sub == "2":
            _, log_p = bio.forward(obs, states, start_p, trans_p, emit_p)
            print(f"\nlog P(obs | model) = {log_p:.4f}")
            print(f"P(obs | model)     = {_math.exp(log_p):.4e}")

        elif sub == "3":
            iters = _read_int("EM iterations", 100)
            sp2, tp2, ep2 = bio.baum_welch(
                obs, states, alphabet, start_p, trans_p, emit_p, iterations=iters
            )
            print("\nRe-estimated parameters:")
            print("  Transitions:")
            for s in states:
                for t in states:
                    print(f"    {s}->{t}: {tp2[s][t]:.4f}")
            print("  Emissions:")
            for s in states:
                for a in alphabet:
                    print(f"    {s} emits {a}: {ep2[s][a]:.4f}")

    elif sub == "4":
        print("Enter DNA sequence (A/C/G/T):")
        seq = input().strip().upper()
        if not seq:
            raise ValueError("Empty sequence.")
        min_len = _read_int("Minimum island length", 10)
        islands = bio.find_cpg_islands(seq, min_length=min_len)
        if islands:
            print(f"\nFound {len(islands)} CpG island(s):")
            for start, end in islands:
                snippet = seq[start : min(start + 20, end + 1)]
                dots = "..." if end - start + 1 > 20 else ""
                print(f"  [{start:>5}, {end:>5}]  len={end-start+1:>4}  {snippet}{dots}")
        else:
            print("No CpG islands detected.")

    else:
        print("Invalid choice.")



MENU = {
    "1": ("Origin Search      (Cap 1)",   run_cap1),
    "2": ("Motif Search       (Cap 2)",   run_cap2),
    "3": ("Genome Assembly    (Cap 3)",   run_cap3),
    "4": ("Protein Tools      (Cap 4)",   run_cap4),
    "5": ("Sequence Alignment (Cap 5)",   run_cap5),
    "6": ("Genome Rearrangements (Cap 6)", run_cap6),
    "7": ("Molecular Evolution   (Cap 7)", run_cap7),
    "8": ("Gene Clustering       (Cap 8)", run_cap8),
    "9": ("Read Mapping / BWT    (Cap 9)", run_cap9),
    "10": ("HMM & Gene Finding   (Cap 10)", run_cap10),
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
