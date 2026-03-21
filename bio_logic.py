"""
bio_logic.py — Core bioinformatics algorithms.

Chapters:
    1. String primitives & OriC analysis
    2. Profile-based motif search
    3. Genome assembly via De Bruijn graphs
    4. Protein / mass-spectrometry tools
    5. Sequence alignment + substitution matrices
    6. Genome rearrangements
"""

import random
import collections

# ---------------------------------------------------------------------------
# CHAPTER 1 — String primitives & OriC
# ---------------------------------------------------------------------------

def hamming_distance(p: str, q: str) -> int:
    """Number of positions where two equal-length strings differ."""
    if len(p) != len(q):
        raise ValueError(f"Strings must have equal length ({len(p)} vs {len(q)})")
    return sum(a != b for a, b in zip(p, q))


def reverse_complement(pattern: str) -> str:
    """Watson-Crick reverse complement of a DNA string."""
    return pattern.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def compute_skew(genome: str) -> list[int]:
    """
    Running #G - #C skew along *genome*.
    Returns a list of length len(genome)+1 (skew[0] = 0).
    """
    delta = {'G': 1, 'C': -1}
    skew = [0]
    for base in genome:
        skew.append(skew[-1] + delta.get(base, 0))
    return skew


def neighbors(pattern: str, d: int) -> set[str]:
    """All DNA strings within Hamming distance *d* of *pattern*."""
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return set("ACGT")
    neighborhood: set[str] = set()
    for text in neighbors(pattern[1:], d):
        if hamming_distance(pattern[1:], text) < d:
            for base in "ACGT":
                neighborhood.add(base + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood


def frequent_words_with_mismatches_and_rc(text: str, k: int, d: int) -> set[str]:
    """
    Most frequent k-mers in *text* allowing up to *d* mismatches,
    counting each k-mer and its reverse complement together.
    """
    counts: dict[str, int] = {}
    for i in range(len(text) - k + 1):
        kmer = text[i : i + k]
        for seq in (kmer, reverse_complement(kmer)):
            for nb in neighbors(seq, d):
                counts[nb] = counts.get(nb, 0) + 1
    if not counts:
        return set()
    max_count = max(counts.values())
    return {p for p, c in counts.items() if c == max_count}


# ---------------------------------------------------------------------------
# CHAPTER 2 — Profile-based motif search
# ---------------------------------------------------------------------------

def create_profile_with_pseudocounts(motifs: list[str]) -> dict[str, list[float]]:
    """
    Laplace +1 pseudocount profile from a list of equal-length motifs.
    Returns {base: [freq_pos_0, ..., freq_pos_k-1]}.
    """
    k = len(motifs[0])
    profile: dict[str, list[float]] = {b: [1.0] * k for b in "ACGT"}
    for motif in motifs:
        for i, base in enumerate(motif):
            profile[base][i] += 1.0
    denom = len(motifs) + 4.0
    for base in "ACGT":
        profile[base] = [v / denom for v in profile[base]]
    return profile


def get_kmer_probability(kmer: str, profile: dict[str, list[float]]) -> float:
    """Probability of *kmer* under *profile*."""
    prob = 1.0
    for i, base in enumerate(kmer):
        prob *= profile[base][i]
    return prob


def find_most_probable_kmer(text: str, k: int, profile: dict[str, list[float]]) -> str:
    """Profile-most-probable k-mer in *text*."""
    best_prob = -1.0
    best_kmer = text[:k]
    for i in range(len(text) - k + 1):
        kmer = text[i : i + k]
        prob = get_kmer_probability(kmer, profile)
        if prob > best_prob:
            best_prob, best_kmer = prob, kmer
    return best_kmer


def get_score(motifs: list[str]) -> int:
    """
    Motif score = total mismatches from the plurality consensus.
    Lower is better.
    """
    k = len(motifs[0])
    score = 0
    for i in range(k):
        col = [m[i] for m in motifs]
        score += len(motifs) - max(col.count(b) for b in "ACGT")
    return score


def randomized_motif_search(dna_list: list[str], k: int, iterations: int = 1000) -> list[str]:
    """
    Run randomized motif search *iterations* times; return the best result found.
    """
    def _single_run() -> list[str]:
        motifs = [seq[random.randint(0, len(seq) - k):][:k] for seq in dna_list]
        best = motifs[:]
        while True:
            profile = create_profile_with_pseudocounts(motifs)
            motifs = [find_most_probable_kmer(seq, k, profile) for seq in dna_list]
            if get_score(motifs) < get_score(best):
                best = motifs[:]
            else:
                return best

    best = _single_run()
    for _ in range(iterations - 1):
        candidate = _single_run()
        if get_score(candidate) < get_score(best):
            best = candidate
    return best


# ---------------------------------------------------------------------------
# CHAPTER 3 — Genome assembly (De Bruijn + Eulerian path)
# ---------------------------------------------------------------------------

def build_de_bruijn_from_kmers(patterns: list[str]) -> dict[str, list[str]]:
    """De Bruijn graph from a list of k-mers."""
    graph: dict[str, list[str]] = {}
    for p in patterns:
        graph.setdefault(p[:-1], []).append(p[1:])
    return graph


def build_paired_de_bruijn(paired_kmers: list[tuple[str, str]]) -> dict[str, list[str]]:
    """
    De Bruijn graph for paired-end reads.
    Each element of *paired_kmers* is a (read1, read2) tuple.
    """
    graph: dict[str, list[str]] = {}
    for r1, r2 in paired_kmers:
        u = f"{r1[:-1]}|{r2[:-1]}"
        v = f"{r1[1:]}|{r2[1:]}"
        graph.setdefault(u, []).append(v)
    return graph


def find_eulerian_path(graph: dict[str, list[str]]) -> list[str]:
    """
    Hierholzer's algorithm for an Eulerian path in a semi-Eulerian graph.
    Falls back to any node as start if the graph is Eulerian (circuit case).
    """
    in_deg: dict[str, int] = collections.defaultdict(int)
    out_deg: dict[str, int] = collections.defaultdict(int)
    nodes: set[str] = set(graph)

    for u, adj in graph.items():
        out_deg[u] += len(adj)
        for v in adj:
            in_deg[v] += 1
            nodes.add(v)

    start = next(
        (n for n in nodes if out_deg[n] - in_deg[n] == 1),
        next(iter(nodes)),
    )

    tmp: dict[str, list[str]] = {u: list(vs) for u, vs in graph.items()}
    stack, path = [start], []
    while stack:
        u = stack[-1]
        if tmp.get(u):
            stack.append(tmp[u].pop())
        else:
            path.append(stack.pop())
    return path[::-1]


def path_to_genome(path: list[str]) -> str:
    """Reconstruct a genome string from a De Bruijn node path."""
    return path[0] + "".join(node[-1] for node in path[1:])


def paired_path_to_genome(path: list[str], k: int, d: int) -> str:
    """
    Reconstruct a genome from a paired De Bruijn path.
    Raises ValueError if the two string reconstructions are inconsistent.
    """
    first  = path[0].split("|")[0] + "".join(n.split("|")[0][-1] for n in path[1:])
    second = path[0].split("|")[1] + "".join(n.split("|")[1][-1] for n in path[1:])
    overlap = len(first) - (k + d)
    if first[k + d:] != second[:overlap]:
        raise ValueError("Paired-end reconstruction failed: overlapping regions do not match.")
    return first + second[overlap:]


# ---------------------------------------------------------------------------
# CHAPTER 4 — Protein / mass spectrometry
# ---------------------------------------------------------------------------

# Integer monoisotopic masses (Da) for the 20 standard amino acids.
MASS_TABLE: dict[str, int] = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186,
}

# Unique residue masses used by the de-novo sequencing algorithm.
AMINO_ACID_MASSES: list[int] = sorted(set(MASS_TABLE.values()))


def get_linear_spectrum(peptide_masses: list[int]) -> list[int]:
    """Theoretical linear spectrum for a peptide given as a list of residue masses."""
    prefix = [0]
    for m in peptide_masses:
        prefix.append(prefix[-1] + m)
    return sorted(
        prefix[j] - prefix[i]
        for i in range(len(prefix))
        for j in range(i + 1, len(prefix))
    )


def get_circular_spectrum(peptide_masses: list[int]) -> list[int]:
    """Theoretical cyclic spectrum for a peptide given as a list of residue masses."""
    prefix = [0]
    for m in peptide_masses:
        prefix.append(prefix[-1] + m)
    total = prefix[-1]
    spectrum = [0]
    for i in range(len(prefix)):
        for j in range(i + 1, len(prefix)):
            mass = prefix[j] - prefix[i]
            spectrum.append(mass)
            if i > 0 and j < len(prefix) - 1:
                spectrum.append(total - mass)
    return sorted(spectrum)


def _is_consistent(peptide_masses: list[int], experimental_spectrum: list[int]) -> bool:
    """True if the linear spectrum of *peptide_masses* is a sub-multiset of *experimental_spectrum*."""
    remaining = list(experimental_spectrum)
    for mass in get_linear_spectrum(peptide_masses):
        try:
            remaining.remove(mass)
        except ValueError:
            return False
    return True


def cyclopeptide_sequencing(experimental_spectrum: list[int]) -> list[list[int]]:
    """
    Branch-and-bound cyclopeptide sequencing.
    Returns all peptide mass-lists whose cyclic spectrum matches *experimental_spectrum*.
    """
    parent_mass = max(experimental_spectrum)
    candidates: list[list[int]] = [[]]
    results: list[list[int]] = []

    while candidates:
        candidates = [c + [m] for c in candidates for m in AMINO_ACID_MASSES]
        survivors: list[list[int]] = []
        for peptide in candidates:
            current_mass = sum(peptide)
            if current_mass == parent_mass:
                if get_circular_spectrum(peptide) == experimental_spectrum:
                    results.append(peptide)
            elif current_mass < parent_mass:
                if _is_consistent(peptide, experimental_spectrum):
                    survivors.append(peptide)
        candidates = survivors

    return results


# ---------------------------------------------------------------------------
# CHAPTER 5 — Sequence alignment + substitution matrices
# ---------------------------------------------------------------------------

# Type alias: maps residue pairs to integer log-odds scores.
SubMatrix = dict[tuple[str, str], int]

# ---- BLOSUM62 --------------------------------------------------------------
_BLOSUM62_RAW: list[tuple] = [
    #      A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    ("A",  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0),
    ("R", -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3),
    ("N", -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3),
    ("D", -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3),
    ("C",  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1),
    ("Q", -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2),
    ("E", -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2),
    ("G",  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3),
    ("H", -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3),
    ("I", -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3),
    ("L", -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1),
    ("K", -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2),
    ("M", -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1),
    ("F", -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1),
    ("P", -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2),
    ("S",  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2),
    ("T",  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0),
    ("W", -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3),
    ("Y", -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1),
    ("V",  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4),
]

# ---- PAM250 ----------------------------------------------------------------
_PAM250_RAW: list[tuple] = [
    #      A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    ("A",  2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0),
    ("R", -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2),
    ("N",  0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2),
    ("D",  0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2),
    ("C", -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2),
    ("Q",  0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2),
    ("E",  0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2),
    ("G",  1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1),
    ("H", -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2),
    ("I", -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4),
    ("L", -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2),
    ("K", -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2),
    ("M", -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2),
    ("F", -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1),
    ("P",  1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1),
    ("S",  1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1),
    ("T",  1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0),
    ("W", -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6),
    ("Y", -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2),
    ("V",  0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4),
]

_AA_ORDER: list[str] = [row[0] for row in _BLOSUM62_RAW]


def _build_matrix(raw: list[tuple]) -> SubMatrix:
    """Expand a triangular raw matrix into a full symmetric (aa1, aa2) -> score dict."""
    order = [row[0] for row in raw]
    m: SubMatrix = {}
    for i, row in enumerate(raw):
        for j, score in enumerate(row[1:]):
            m[(order[i], order[j])] = score
            m[(order[j], order[i])] = score
    return m


BLOSUM62: SubMatrix = _build_matrix(_BLOSUM62_RAW)
PAM250:   SubMatrix = _build_matrix(_PAM250_RAW)
IDENTITY: SubMatrix = {
    (a, b): (1 if a == b else -1)
    for a in (*"ACGT", *_AA_ORDER)
    for b in (*"ACGT", *_AA_ORDER)
}

SUBSTITUTION_MATRICES: dict[str, SubMatrix] = {
    "BLOSUM62": BLOSUM62,
    "PAM250":   PAM250,
    "IDENTITY": IDENTITY,
}


def _sub(matrix: SubMatrix | None, a: str, b: str, match: int, mismatch: int) -> int:
    """Substitution score for residues *a* and *b* under *matrix* (or flat scores)."""
    if matrix is not None:
        return matrix.get((a, b), mismatch)
    return match if a == b else mismatch


def _traceback_global(
    dp: list[list[int]],
    s1: str, s2: str,
    matrix: SubMatrix | None,
    match: int, mismatch: int, gap: int,
) -> tuple[str, str]:
    """Traceback for Needleman-Wunsch."""
    a1: list[str] = []
    a2: list[str] = []
    i, j = len(s1), len(s2)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + _sub(matrix, s1[i-1], s2[j-1], match, mismatch):
            a1.append(s1[i-1]); a2.append(s2[j-1]); i -= 1; j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
            a1.append(s1[i-1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(s2[j-1]); j -= 1
    return "".join(reversed(a1)), "".join(reversed(a2))


def _traceback_local(
    dp: list[list[int]],
    s1: str, s2: str,
    best_pos: tuple[int, int],
    matrix: SubMatrix | None,
    match: int, mismatch: int, gap: int,
) -> tuple[str, str]:
    """Traceback for Smith-Waterman."""
    a1: list[str] = []
    a2: list[str] = []
    i, j = best_pos
    while i > 0 and j > 0 and dp[i][j] > 0:
        if dp[i][j] == dp[i-1][j-1] + _sub(matrix, s1[i-1], s2[j-1], match, mismatch):
            a1.append(s1[i-1]); a2.append(s2[j-1]); i -= 1; j -= 1
        elif dp[i][j] == dp[i-1][j] + gap:
            a1.append(s1[i-1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(s2[j-1]); j -= 1
    return "".join(reversed(a1)), "".join(reversed(a2))


def global_alignment(
    s1: str,
    s2: str,
    matrix: SubMatrix | None = None,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
) -> tuple[str, str, int]:
    """
    Needleman-Wunsch global alignment.

    Args:
        s1, s2:   sequences to align.
        matrix:   substitution matrix (BLOSUM62, PAM250, IDENTITY, or None).
                  When None, flat *match* / *mismatch* scores are used.
        gap:      linear gap penalty (negative integer).

    Returns:
        (aligned_s1, aligned_s2, score)
    """
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = dp[i-1][j-1] + _sub(matrix, s1[i-1], s2[j-1], match, mismatch)
            dp[i][j] = max(diag, dp[i-1][j] + gap, dp[i][j-1] + gap)
    a1, a2 = _traceback_global(dp, s1, s2, matrix, match, mismatch, gap)
    return a1, a2, dp[m][n]


def local_alignment(
    s1: str,
    s2: str,
    matrix: SubMatrix | None = None,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
) -> tuple[str, str, int]:
    """
    Smith-Waterman local alignment.

    Args:
        s1, s2:   sequences to align.
        matrix:   substitution matrix (BLOSUM62, PAM250, IDENTITY, or None).
        gap:      linear gap penalty (negative integer).

    Returns:
        (aligned_s1, aligned_s2, score)
    """
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    best_score = 0
    best_pos = (0, 0)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = dp[i-1][j-1] + _sub(matrix, s1[i-1], s2[j-1], match, mismatch)
            dp[i][j] = max(0, diag, dp[i-1][j] + gap, dp[i][j-1] + gap)
            if dp[i][j] >= best_score:
                best_score = dp[i][j]
                best_pos = (i, j)
    a1, a2 = _traceback_local(dp, s1, s2, best_pos, matrix, match, mismatch, gap)
    return a1, a2, best_score


# ---------------------------------------------------------------------------
# CHAPTER 6 — Genome rearrangements
# ---------------------------------------------------------------------------

def _validate_signed_permutation(p: list[int]) -> None:
    """
    Raise ValueError if *p* is not a valid signed permutation.
    A valid signed permutation of length n contains each of +/-1 ... +/-n exactly once.
    """
    if not p:
        raise ValueError("Permutation cannot be empty.")
    n = len(p)
    if sorted(abs(x) for x in p) != list(range(1, n + 1)):
        raise ValueError(
            f"Not a valid signed permutation of length {n}: "
            f"expected each of {{+/-1 ... +/-{n}}} exactly once, got {p}."
        )


def greedy_sorting(p: list[int]) -> list[list[int]]:
    """
    Greedy sorting by reversals (approximation algorithm).

    Sorts a signed permutation to the identity (+1, +2, ..., +n) by greedily
    bringing each target element into position via a reversal, then fixing its
    sign with a second reversal if needed.

    Returns the list of intermediate permutations (one per reversal applied),
    not including the initial state.

    Raises ValueError for invalid input.
    """
    _validate_signed_permutation(p)
    p = list(p)
    steps: list[list[int]] = []

    for i in range(len(p)):
        target = i + 1
        try:
            idx = next(j for j in range(i, len(p)) if abs(p[j]) == target)
        except StopIteration:
            raise ValueError(
                f"Element +/-{target} not found at or after position {i} -- malformed permutation."
            )
        if idx != i or p[i] != target:
            p[i : idx + 1] = [-x for x in reversed(p[i : idx + 1])]
            steps.append(list(p))
        if p[i] == -target:
            p[i] = target
            steps.append(list(p))

    return steps


def count_breakpoints(p: list[int]) -> int:
    """
    Count breakpoints in signed permutation *p*.

    A breakpoint is any adjacent pair in (0, p[0], ..., p[n-1], n+1)
    where consecutive elements do not differ by exactly +1.

    Raises ValueError for invalid input.
    """
    _validate_signed_permutation(p)
    extended = [0, *p, len(p) + 1]
    return sum(extended[i + 1] - extended[i] != 1 for i in range(len(extended) - 1))


# ---------------------------------------------------------------------------
# CHAPTER 7 — Molecular evolution (phylogenetics)
# ---------------------------------------------------------------------------
#
# Algorithms covered:
#   - UPGMA          (distance-based, ultrametric tree)
#   - Neighbor Joining (distance-based, unrooted additive tree)
#   - Small Parsimony  (Sankoff / Fitch on a fixed unrooted tree)
#
# Internal tree representation
# ----------------------------
# An unrooted tree is stored as a dict[int, dict[int, float]]:
#   adj[u][v] = branch_length_u_to_v
# Leaves are labelled 0 … n-1; internal nodes start at n.
# ---------------------------------------------------------------------------

import numpy as np

# Type alias used throughout this chapter.
Tree = dict[int, dict[int, float]]


# ---- helpers ---------------------------------------------------------------

def _tree_add_edge(tree: Tree, u: int, v: int, w: float) -> None:
    tree.setdefault(u, {})[v] = w
    tree.setdefault(v, {})[u] = w


def _tree_remove_edge(tree: Tree, u: int, v: int) -> None:
    tree[u].pop(v, None)
    tree[v].pop(u, None)


def _newick(tree: Tree, node: int, parent: int | None, labels: list[str]) -> str:
    """Recursive Newick serialisation (works for rooted & unrooted)."""
    n_leaves = len(labels)
    children = [v for v in tree[node] if v != parent]
    if not children:
        return labels[node] if node < n_leaves else f"n{node}"
    parts = []
    for c in children:
        w = tree[node][c]
        parts.append(f"{_newick(tree, c, node, labels)}:{w:.4f}")
    name = labels[node] if node < n_leaves else f"n{node}"
    return f"({','.join(parts)}){name}"


# ---- UPGMA -----------------------------------------------------------------

def upgma(
    distance_matrix: list[list[float]],
    labels: list[str],
) -> tuple[Tree, list[tuple[int, int, float]]]:
    """
    UPGMA (Unweighted Pair Group Method with Arithmetic mean).

    Builds a rooted ultrametric tree from a symmetric distance matrix.

    Args:
        distance_matrix: n×n symmetric matrix (list of lists or 2-D array).
        labels:          sequence labels, length n.

    Returns:
        (tree, edges) where
          tree  — adjacency dict {node: {neighbour: branch_length}}
          edges — list of (parent, child, branch_length) tuples
    """
    n = len(labels)
    if len(distance_matrix) != n or any(len(r) != n for r in distance_matrix):
        raise ValueError("distance_matrix dimensions must match len(labels).")

    mat = np.array(distance_matrix, dtype=float)
    np.fill_diagonal(mat, np.inf)   # ignore self-distances

    cluster_size = {i: 1 for i in range(n)}
    node_height  = {i: 0.0 for i in range(n)}
    active       = list(range(n))
    tree: Tree   = {i: {} for i in range(n)}
    edges: list[tuple[int, int, float]] = []
    next_node    = n

    # index_of[active[k]] = row/col k in the current submatrix
    # We track active node ids; mat rows/cols correspond to active list positions.

    while len(active) > 1:
        # Find minimum off-diagonal entry
        size = len(active)
        sub = mat[:size, :size]
        flat_idx = int(np.argmin(sub))
        ri, ci = divmod(flat_idx, size)
        if ri > ci:
            ri, ci = ci, ri

        u, v = active[ri], active[ci]
        min_dist = mat[ri, ci]

        # New internal node
        new = next_node
        next_node += 1
        node_height[new] = min_dist / 2.0
        tree[new] = {}

        branch_u = node_height[new] - node_height[u]
        branch_v = node_height[new] - node_height[v]
        _tree_add_edge(tree, new, u, branch_u)
        _tree_add_edge(tree, new, v, branch_v)
        edges.append((new, u, branch_u))
        edges.append((new, v, branch_v))

        # Compute new distances (weighted average)
        su, sv = cluster_size[u], cluster_size[v]
        new_dists = np.array([
            (mat[ri, k] * su + mat[ci, k] * sv) / (su + sv)
            for k in range(size)
        ])

        # Shrink matrix: delete rows/cols for ri and ci, append new cluster
        keep = [k for k in range(size) if k not in (ri, ci)]
        mat = mat[np.ix_(keep, keep)]
        new_row = new_dists[keep]
        s = mat.shape[0]
        new_mat = np.empty((s + 1, s + 1))
        new_mat[:s, :s] = mat
        new_mat[s, :s]  = new_row
        new_mat[:s, s]  = new_row
        new_mat[s, s]   = np.inf
        mat = new_mat

        cluster_size[new] = su + sv
        node_height[new]  = min_dist / 2.0

        # Update active list
        active = [active[k] for k in keep]
        active.append(new)

    return tree, edges


# ---- Neighbor Joining ------------------------------------------------------

def neighbor_joining(
    distance_matrix: list[list[float]],
    labels: list[str],
) -> tuple[Tree, list[tuple[int, int, float]]]:
    """
    Neighbor Joining (Saitou & Nei, 1987).

    Builds an unrooted additive tree from a symmetric distance matrix.

    Args:
        distance_matrix: n×n symmetric distance matrix.
        labels:          sequence labels, length n.

    Returns:
        (tree, edges) where
          tree  — adjacency dict {node: {neighbour: branch_length}}
          edges — list of (u, v, branch_length) tuples (undirected)
    """
    n = len(labels)
    if len(distance_matrix) != n or any(len(r) != n for r in distance_matrix):
        raise ValueError("distance_matrix dimensions must match len(labels).")

    mat = np.array(distance_matrix, dtype=float)
    active = list(range(n))
    tree: Tree = {i: {} for i in range(n)}
    edges: list[tuple[int, int, float]] = []
    next_node = n

    while len(active) > 2:
        size = len(active)

        # Compute the Q-matrix
        row_sums = mat[:size, :size].sum(axis=1)
        Q = np.full((size, size), np.inf)
        for i in range(size):
            for j in range(i + 1, size):
                Q[i, j] = (size - 2) * mat[i, j] - row_sums[i] - row_sums[j]
                Q[j, i] = Q[i, j]

        # Find minimum Q entry
        flat_idx = int(np.argmin(Q))
        ri, ci = divmod(flat_idx, size)
        if ri > ci:
            ri, ci = ci, ri

        u, v = active[ri], active[ci]
        d_uv = mat[ri, ci]

        # Branch lengths from new node to u and v
        if size > 2:
            delta = (row_sums[ri] - row_sums[ci]) / (size - 2)
        else:
            delta = 0.0
        branch_u = (d_uv + delta) / 2.0
        branch_v = d_uv - branch_u

        # New internal node
        new = next_node
        next_node += 1
        tree[new] = {}
        _tree_add_edge(tree, new, u, branch_u)
        _tree_add_edge(tree, new, v, branch_v)
        edges.append((new, u, branch_u))
        edges.append((new, v, branch_v))

        # Compute distances from new node to remaining nodes
        new_dists = np.array([
            (mat[ri, k] + mat[ci, k] - d_uv) / 2.0
            for k in range(size)
        ])

        # Shrink matrix
        keep = [k for k in range(size) if k not in (ri, ci)]
        mat = mat[np.ix_(keep, keep)]
        new_row = new_dists[keep]
        s = mat.shape[0]
        new_mat = np.zeros((s + 1, s + 1))
        new_mat[:s, :s] = mat
        new_mat[s, :s]  = new_row
        new_mat[:s, s]  = new_row
        mat = new_mat

        active = [active[k] for k in keep]
        active.append(new)

    # Connect the last two nodes
    if len(active) == 2:
        u, v = active[0], active[1]
        w = float(mat[0, 1])
        _tree_add_edge(tree, u, v, w)
        edges.append((u, v, w))

    return tree, edges


# ---- Small Parsimony (Fitch / Sankoff) ------------------------------------

def small_parsimony(
    tree: Tree,
    leaf_chars: dict[int, str],
    alphabet: str = "ACGT",
) -> tuple[int, dict[int, str]]:
    """
    Small Parsimony on a fixed unrooted binary tree (Fitch algorithm).

    Works by temporarily rooting at an arbitrary internal node, running the
    two-pass Fitch algorithm, then un-rooting.

    Args:
        tree:       adjacency dict from upgma() or neighbor_joining().
        leaf_chars: {leaf_node: character} for each leaf (single site).
        alphabet:   characters to consider (default "ACGT").

    Returns:
        (parsimony_score, node_assignments) where
          parsimony_score  — minimum number of substitutions
          node_assignments — {node: character} for every node
    """
    if not tree:
        raise ValueError("Tree is empty.")

    leaves = set(leaf_chars)
    internal = [v for v in tree if v not in leaves]
    if not internal:
        raise ValueError("Tree has no internal nodes.")

    # Root at the first internal node found
    root = internal[0]

    # Build a parent map via BFS from root (treats tree as directed)
    parent: dict[int, int | None] = {root: None}
    order: list[int] = []
    queue = collections.deque([root])
    while queue:
        node = queue.popleft()
        order.append(node)
        for nb in tree[node]:
            if nb not in parent:
                parent[nb] = node
                queue.append(nb)

    # ----- Downward pass (leaves → root): compute Fitch sets -----
    fitch: dict[int, set[str]] = {}

    score = 0
    for node in reversed(order):
        if node in leaves:
            fitch[node] = {leaf_chars[node]}
            continue
        children = [nb for nb in tree[node] if nb != parent[node]]
        if not children:
            fitch[node] = set(alphabet)
            continue
        child_sets = [fitch[c] for c in children]
        intersection = child_sets[0].intersection(*child_sets[1:])
        if intersection:
            fitch[node] = intersection
        else:
            fitch[node] = child_sets[0].union(*child_sets[1:])
            score += 1

    # ----- Upward pass (root → leaves): assign characters -----
    assignment: dict[int, str] = {}

    for node in order:
        if node in leaves:
            assignment[node] = leaf_chars[node]
            continue
        par = parent[node]
        candidates = fitch[node]
        if par is None:
            # Root: pick any candidate
            assignment[node] = min(candidates)
        else:
            par_char = assignment[par]
            assignment[node] = par_char if par_char in candidates else min(candidates)

    return score, assignment


def parsimony_score(
    tree: Tree,
    leaf_sequences: dict[int, str],
    alphabet: str = "ACGT",
) -> tuple[int, dict[int, list[str]]]:
    """
    Total parsimony score over all sites of aligned sequences.

    Args:
        tree:            adjacency dict.
        leaf_sequences:  {leaf_node: sequence_string} (equal-length strings).
        alphabet:        characters to consider.

    Returns:
        (total_score, node_sequences) where node_sequences maps every node to
        its reconstructed sequence as a list of characters.
    """
    lengths = {len(s) for s in leaf_sequences.values()}
    if len(lengths) != 1:
        raise ValueError("All leaf sequences must have the same length.")
    seq_len = lengths.pop()

    total = 0
    node_seqs: dict[int, list[str]] = {n: [] for n in tree}

    for site in range(seq_len):
        leaf_chars = {node: seq[site] for node, seq in leaf_sequences.items()}
        site_score, assignment = small_parsimony(tree, leaf_chars, alphabet)
        total += site_score
        for node, char in assignment.items():
            node_seqs[node].append(char)

    return total, node_seqs


def newick(tree: Tree, labels: list[str]) -> str:
    """
    Serialise *tree* in Newick format.
    Roots at the node with the highest index (last internal node created).
    """
    root = max(tree)
    return _newick(tree, root, None, labels) + ";"
