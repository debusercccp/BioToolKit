#!/usr/bin/env python3
"""
main.py — CLI entry point for BioToolkit.

Each chapter is an independent handler function.
Run:  python main.py
"""
import re
import sys
import math as _math
import os
import bio_logic as bio

try:
    import matplotlib
    matplotlib.use('Agg')              
    import matplotlib.pyplot as plt
    _HAS_MATPLOTLIB = True
    _HAS_DISPLAY = bool(
        os.environ.get('DISPLAY') or
        os.environ.get('WAYLAND_DISPLAY')
    )
except ImportError:
    _HAS_MATPLOTLIB = False
    _HAS_DISPLAY = False

try:
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    _HAS_BIOPYTHON = True
except ImportError:
    _HAS_BIOPYTHON = False
# ---------------------------------------------------------------------------
# Input helpers — robusti
# ---------------------------------------------------------------------------

_ANSI_RE = re.compile(r'\x1b[\x40-\x5f]|'   # ESC + Fe
                       r'\x1b\[[0-?]*[ -/]*[@-~]')  # CSI sequences


def _clean(raw: str, valid: str = "ACGT") -> str:
    """Rimuove ANSI escape, spazi e caratteri non validi; uppercase."""
    s = _ANSI_RE.sub('', raw)
    return ''.join(c for c in s.upper() if c in valid)


def _read_fasta_first(path: str) -> str:
    """Legge la prima sequenza da un file FASTA."""
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        raise ValueError(f"File non trovato: {path}")
    if _HAS_BIOPYTHON:
        with open(path) as fh:
            _, seq = next(SimpleFastaParser(fh), (None, None))
        if not seq:
            raise ValueError("Nessuna sequenza nel file FASTA.")
        return seq.upper()
    # fallback senza biopython
    lines, in_seq = [], False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if in_seq:
                    break
                in_seq = True
            elif in_seq:
                lines.append(line)
    if not lines:
        raise ValueError("Nessuna sequenza nel file FASTA.")
    return ''.join(lines).upper()


def _read_fasta_all(path: str) -> list[str]:
    """Legge tutte le sequenze da un file FASTA."""
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        raise ValueError(f"File non trovato: {path}")
    seqs = []
    if _HAS_BIOPYTHON:
        with open(path) as fh:
            for _, seq in SimpleFastaParser(fh):
                seqs.append(seq.upper())
    else:
        cur: list[str] = []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith('>'):
                    if cur:
                        seqs.append(''.join(cur).upper())
                    cur = []
                elif line:
                    cur.append(line)
        if cur:
            seqs.append(''.join(cur).upper())
    if not seqs:
        raise ValueError("Nessuna sequenza nel file FASTA.")
    return seqs


def _read_dna(prompt: str = "DNA sequence") -> str:
    """
    Legge una sequenza DNA da stdin.
    Accetta:
      - stringa diretta (anche con ANSI escape dal terminale)
      - percorso a file FASTA (se la riga inizia con / o ~ o ./)
      - sequenza multi-riga: termina con riga vuota
    """
    print(f"{prompt}")
    print("  [sequenza diretta | percorso FASTA | multi-riga → invio vuoto per terminare]")
    first = input("  > ").strip()

    # File FASTA
    if first.startswith(('/','~','./')) or (os.path.sep in first):
        return _read_fasta_first(first)

    # Multi-riga: se la prima riga sembra incompleta, continua a leggere
    parts = [_clean(first)]
    while True:
        line = input("  > " if not parts[0] else "    ").strip()
        if not line:
            break
        parts.append(_clean(line))

    cleaned = ''.join(parts)
    if not cleaned:
        raise ValueError("Sequenza vuota o non valida (caratteri ammessi: A C G T).")
    return cleaned


def _read_seq(prompt: str = "Sequence", valid: str = "ACDEFGHIKLMNPQRSTVWY") -> str:
    """Legge una sequenza generica (DNA o proteina) con cleaning ANSI."""
    print(f"{prompt}")
    print("  [sequenza diretta | percorso FASTA | multi-riga → invio vuoto per terminare]")
    first = input("  > ").strip()

    if first.startswith(('/','~','./')) or (os.path.sep in first):
        return _read_fasta_first(first)

    parts = [_clean(first, valid + "ACGT")]   # accetta anche DNA
    while True:
        line = input("    ").strip()
        if not line:
            break
        parts.append(_clean(line, valid + "ACGT"))

    cleaned = ''.join(parts)
    if not cleaned:
        raise ValueError("Sequenza vuota.")
    return cleaned


def _read_dna_list() -> list[str]:
    """
    Legge una lista di sequenze DNA.
    Accetta anche percorso FASTA (tutte le sequenze del file).
    """
    print("Sequenze DNA — una per riga, riga vuota per terminare")
    print("  [oppure: percorso FASTA per caricare tutte le sequenze del file]")
    first = input("  > ").strip()

    if first.startswith(('/','~','./')) or (os.path.sep in first):
        return _read_fasta_all(first)

    seqs: list[str] = []
    line = first
    while True:
        cleaned = _clean(line)
        if cleaned:
            seqs.append(cleaned)
        line = input("  > ").strip()
        if not line:
            break

    if not seqs:
        raise ValueError("Nessuna sequenza fornita.")
    return seqs


def _read_int(prompt: str, default: int) -> int:
    raw = input(f"{prompt} (default {default}): ").strip()
    return int(raw) if raw else default


def _read_text(prompt: str) -> str:
    """Legge testo generico (BWT, patterns) con ANSI cleaning ma senza filtro ACGT."""
    print(f"{prompt}")
    first = input("  > ").strip()
    if first.startswith(('/','~','./')) or (os.path.sep in first):
        return _read_fasta_first(first)
    parts = [_ANSI_RE.sub('', first).strip()]
    while True:
        line = input("    ").strip()
        if not line:
            break
        parts.append(_ANSI_RE.sub('', line).strip())
    result = ''.join(parts).upper()
    if not result:
        raise ValueError("Input vuoto.")
    return result


# ---------------------------------------------------------------------------
# Chapter handlers
# ---------------------------------------------------------------------------

def run_cap1() -> None:
    """OriC finder: skew minimum + DnaA box search."""
    print("\n--- [CAP 1] ORIGIN FINDER ---")
    genome = _read_dna("Genoma (DNA)")
    skew = bio.compute_skew(genome)
    min_val = min(skew)
    min_positions = [i for i, v in enumerate(skew) if v == min_val]
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

        if _HAS_DISPLAY:
            plt.show()
        else:
            out = "skew_plot.png"
            plt.savefig(out, dpi=150)
            plt.close()
            print(f"(nessun display — grafico salvato in '{out}')")
    else:
        print("(matplotlib non installato — skipping plot)")


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
            kmers.append(_ANSI_RE.sub('', line).strip().upper())
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
            line = _ANSI_RE.sub('', input().strip()).upper()
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


def run_cap4() -> None:
    """Protein / mass-spectrometry tools."""
    print("\n--- [CAP 4] PROTEIN TOOLS ---")
    print("  1) Peptide → theoretical spectrum")
    print("  2) Spectrum → cyclopeptide sequencing")
    sub = input("Choice: ").strip()

    if sub == "1":
        peptide = _read_seq("Amino acid sequence (e.g. LEQN)", valid="ACDEFGHIKLMNPQRSTVWY")
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
        raw = _ANSI_RE.sub('', input("Experimental spectrum (space-separated integers): ").strip())
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


def run_cap5() -> None:
    """Sequence alignment: global (Needleman-Wunsch) or local (Smith-Waterman)."""
    print("\n--- [CAP 5] SEQUENCE ALIGNMENT ---")
    s1 = _read_seq("Sequence 1")
    s2 = _read_seq("Sequence 2")

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
            gap = -4
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
    raw = _ANSI_RE.sub('', input().strip())
    if not raw:
        print("Empty input.")
        return

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
    return "(" + " ".join(f"+{x}" if x > 0 else str(x) for x in p) + ")"


def _read_distance_matrix() -> tuple[list[list[float]], list[str]]:
    print("Enter labels (space-separated, e.g.  A B C D):")
    raw_labels = input().strip().split()
    if not raw_labels:
        raise ValueError("No labels provided.")
    n = len(raw_labels)
    print(f"Enter {n} rows of {n} distances (one row per line):")
    matrix: list[list[float]] = []
    for i in range(n):
        row_raw = _ANSI_RE.sub('', input().strip()).split()
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
            seq = _read_dna(f"  {label}")
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
        line = _ANSI_RE.sub('', input().strip())
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
        k = _read_int("Number of clusters k", 2)
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
        text = _read_text("Testo di riferimento (include '$' se già presente)")
        bwt = bio.bwt_transform(text)
        sa  = bio.suffix_array(text)
        print(f"BWT : {bwt}")
        print(f"SA  : {sa}")

    elif sub == "2":
        print("BWT string:")
        bwt = _ANSI_RE.sub('', input("  > ").strip())
        if not bwt:
            raise ValueError("Empty input.")
        original = bio.bwt_inverse(bwt)
        print(f"Original: {original}")

    elif sub == "3":
        text    = _read_text("Reference text")
        pattern = _read_dna("Pattern")
        hits = bio.bwt_match(text, pattern)
        if hits:
            print(f"Found {len(hits)} occurrence(s) at positions: {hits}")
        else:
            print("Pattern not found.")

    elif sub == "4":
        text    = _read_text("Reference text")
        pattern = _read_dna("Pattern")
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

        raw = _ANSI_RE.sub('', input("  > ").strip())
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
        seq = _read_dna("DNA sequence per CpG island detection")
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


def run_cap11() -> None:
    """Advanced alignment: affine gaps, linear space, progressive MSA, tandem repeats."""
    print("\n--- [CAP 11] ADVANCED ALIGNMENT ---")
    print("  1) Affine gap global alignment")
    print("  2) Affine gap local alignment")
    print("  3) Linear space alignment (Hirschberg)")
    print("  4) Multiple sequence alignment (progressive, ClustalW-style)")
    print("  5) Tandem repeats finder")
    sub = input("Choice: ").strip()

    if sub in ("1", "2", "3"):
        s1 = _read_seq("Sequence 1")
        s2 = _read_seq("Sequence 2")

        print("  Matrix:  1) BLOSUM62  2) PAM250  3) IDENTITY  4) None (flat)")
        mc = input("  Choice (default 1): ").strip()
        mat_map = {"1": "BLOSUM62", "2": "PAM250", "3": "IDENTITY"}
        matrix = bio.SUBSTITUTION_MATRICES.get(mat_map.get(mc, "BLOSUM62"))
        mat_label = mat_map.get(mc, "BLOSUM62") if mc in mat_map else "flat"

        go_raw = input("Gap open   (default -11): ").strip()
        ge_raw = input("Gap extend (default  -1): ").strip()
        gap_open   = int(go_raw) if go_raw else -11
        gap_extend = int(ge_raw) if ge_raw else -1

        if sub == "1":
            a1, a2, score = bio.affine_global_alignment(
                s1, s2, matrix=matrix, gap_open=gap_open, gap_extend=gap_extend)
            label = "AFFINE GLOBAL"
        elif sub == "2":
            a1, a2, score = bio.affine_local_alignment(
                s1, s2, matrix=matrix, gap_open=gap_open, gap_extend=gap_extend)
            label = "AFFINE LOCAL"
        else:
            a1, a2, score = bio.linear_space_alignment(
                s1, s2, matrix=matrix, gap_open=gap_open, gap_extend=gap_extend)
            label = "HIRSCHBERG"

        mid = "".join("|" if a1[k] == a2[k] and a1[k] != "-" else " " for k in range(len(a1)))
        print(f"\n{label}  [matrix={mat_label}  go={gap_open}  ge={gap_extend}]  score={score}")
        print(f"  S1: {a1}")
        print(f"      {mid}")
        print(f"  S2: {a2}")

    elif sub == "4":
        seqs:   list[str] = []
        labels: list[str] = []
        print("Enter sequences as  label:SEQUENCE  (empty line to stop)")
        print("  [oppure percorso FASTA per caricare tutte le sequenze]")
        first = input("  > ").strip()
        if first.startswith(('/','~','./')) or (os.path.sep in first):
            seqs = _read_fasta_all(first)
            labels = [f"seq{i+1}" for i in range(len(seqs))]
        else:
            line = first
            while True:
                line_clean = _ANSI_RE.sub('', line).strip()
                if not line_clean:
                    break
                if ":" in line_clean:
                    lab, seq = line_clean.split(":", 1)
                    labels.append(lab.strip())
                    seqs.append(_clean(seq, "ACDEFGHIKLMNPQRSTVWYACGT"))
                else:
                    seqs.append(_clean(line_clean, "ACDEFGHIKLMNPQRSTVWYACGT"))
                    labels.append(f"seq{len(seqs)}")
                line = input("  > ").strip()

        if len(seqs) < 2:
            raise ValueError("Need at least 2 sequences.")

        print("  Matrix:  1) BLOSUM62  2) PAM250  3) None (flat)")
        mc = input("  Choice (default 1): ").strip()
        mat_map = {"1": "BLOSUM62", "2": "PAM250"}
        matrix = bio.SUBSTITUTION_MATRICES.get(mat_map.get(mc, "BLOSUM62"))

        go_raw = input("Gap open   (default -11): ").strip()
        ge_raw = input("Gap extend (default  -1): ").strip()
        gap_open   = int(go_raw) if go_raw else -11
        gap_extend = int(ge_raw) if ge_raw else -1

        aligned, sp, guide_tree = bio.multiple_sequence_alignment(
            seqs, labels=labels, matrix=matrix,
            gap_open=gap_open, gap_extend=gap_extend,
        )
        max_lbl = max(len(l) for l in labels)
        print(f"\nProgressive MSA  (SP score={sp}):")
        for lab, s in zip(labels, aligned):
            print(f"  {lab:<{max_lbl}}  {s}")
        print(f"\nGuide tree: {guide_tree}")

    elif sub == "5":
        # Cap 11 tandem repeats — il bug originale: ANSI da incolla lungo
        seq = _read_dna("DNA sequence per tandem repeat finder")
        min_p   = _read_int("Min period", 1)
        max_p   = _read_int("Max period", 500)
        min_c_raw = input("Min copies (default 1.5): ").strip()
        min_c   = float(min_c_raw) if min_c_raw else 1.5

        # Validazione esplicita con messaggio utile
        invalid = set(seq) - set("ACGT")
        if invalid:
            raise ValueError(
                f"Sequenza contiene caratteri non validi: {invalid}\n"
                f"  Usa _read_dna() per pulire automaticamente l'input."
            )

        repeats = bio.find_tandem_repeats(seq, min_period=min_p, max_period=max_p, min_copies=min_c)
        if not repeats:
            print("No tandem repeats found.")
        else:
            print(f"\nFound {len(repeats)} tandem repeat(s):")
            print(f"  {'start':>6}  {'end':>6}  {'len':>5}  {'period':>6}  {'copies':>6}  unit")
            for r in repeats:
                print(f"  {r['start']:>6}  {r['end']:>6}  {r['end']-r['start']+1:>5}  "
                      f"{r['period']:>6}  {r['copies']:>6}  {r['unit']}")
    else:
        print("Invalid choice.")


# ---------------------------------------------------------------------------
# Menu
# ---------------------------------------------------------------------------

MENU = {
    "1":  ("Origin Search         (Cap 1)",  run_cap1),
    "2":  ("Motif Search          (Cap 2)",  run_cap2),
    "3":  ("Genome Assembly       (Cap 3)",  run_cap3),
    "4":  ("Protein Tools         (Cap 4)",  run_cap4),
    "5":  ("Sequence Alignment    (Cap 5)",  run_cap5),
    "6":  ("Genome Rearrangements (Cap 6)",  run_cap6),
    "7":  ("Molecular Evolution   (Cap 7)",  run_cap7),
    "8":  ("Gene Clustering       (Cap 8)",  run_cap8),
    "9":  ("Read Mapping / BWT    (Cap 9)",  run_cap9),
    "10": ("HMM & Gene Finding    (Cap 10)", run_cap10),
    "11": ("Advanced Alignment    (Cap 11)", run_cap11),
}


def main() -> None:
    while True:
        print("\n" + "=" * 36)
        print("         B I O T O O L K I T")
        print("=" * 36)
        for key, (label, _) in MENU.items():
            print(f"  {key:>2}) {label}")
        print("   q) Quit")

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
