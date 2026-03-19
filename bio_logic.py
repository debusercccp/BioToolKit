"""
bio_logic.py — Core bioinformatics algorithms.

Chapters:
    1. String / OriC analysis
    2. Profile-based motif search
    3. Genome assembly via De Bruijn graphs
    4. Protein / mass-spectrometry tools
"""

import random
import collections
from typing import Optional

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
    complement = str.maketrans("ACGT", "TGCA")
    return pattern.translate(complement)[::-1]


def compute_skew(genome: str) -> list[int]:
    """
    Running #G - #C skew along *genome*.
    Returns a list of length len(genome)+1 (index 0 = 0).
    """
    skew = [0]
    delta = {'G': 1, 'C': -1}
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
        diff = hamming_distance(pattern[1:], text)
        if diff < d:
            for base in "ACGT":
                neighborhood.add(base + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood


def frequent_words_with_mismatches_and_rc(
    text: str, k: int, d: int
) -> set[str]:
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
    t = len(motifs)
    k = len(motifs[0])
    profile: dict[str, list[float]] = {b: [1.0] * k for b in "ACGT"}
    for motif in motifs:
        for i, base in enumerate(motif):
            profile[base][i] += 1.0
    denom = t + 4.0
    for base in "ACGT":
        profile[base] = [v / denom for v in profile[base]]
    return profile


def get_kmer_probability(kmer: str, profile: dict[str, list[float]]) -> float:
    """Probability of *kmer* under *profile*."""
    prob = 1.0
    for i, base in enumerate(kmer):
        prob *= profile[base][i]
    return prob


def find_most_probable_kmer(
    text: str, k: int, profile: dict[str, list[float]]
) -> str:
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
    Score = total mismatches from the plurality consensus.
    Lower is better.
    """
    k = len(motifs[0])
    score = 0
    for i in range(k):
        col = [m[i] for m in motifs]
        score += len(motifs) - max(col.count(b) for b in "ACGT")
    return score


def randomized_motif_search(
    dna_list: list[str], k: int, iterations: int = 1000
) -> list[str]:
    """
    Run randomized motif search *iterations* times and return the best result.
    """
    t = len(dna_list)

    def _single_run() -> list[str]:
        motifs = [seq[random.randint(0, len(seq) - k) : ][:k] for seq in dna_list]
        best = motifs[:]
        while True:
            profile = create_profile_with_pseudocounts(motifs)
            motifs = [find_most_probable_kmer(seq, k, profile) for seq in dna_list]
            if get_score(motifs) < get_score(best):
                best = motifs[:]
            else:
                return best

    best_motifs = _single_run()
    for _ in range(iterations - 1):
        candidate = _single_run()
        if get_score(candidate) < get_score(best_motifs):
            best_motifs = candidate
    return best_motifs


# ---------------------------------------------------------------------------
# CHAPTER 3 — Genome assembly (De Bruijn + Eulerian path)
# ---------------------------------------------------------------------------

def build_de_bruijn_from_kmers(patterns: list[str]) -> dict[str, list[str]]:
    """De Bruijn graph from a list of k-mers."""
    graph: dict[str, list[str]] = {}
    for p in patterns:
        prefix, suffix = p[:-1], p[1:]
        graph.setdefault(prefix, []).append(suffix)
    return graph


def build_paired_de_bruijn(
    paired_kmers: list[tuple[str, str]]
) -> dict[str, list[str]]:
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
    Hierholzer's algorithm for an Eulerian *path* (semi-Eulerian graph).
    Raises ValueError if no valid start node is found.
    """
    in_deg: dict[str, int] = collections.defaultdict(int)
    out_deg: dict[str, int] = collections.defaultdict(int)
    nodes: set[str] = set(graph)

    for u, neighbors_list in graph.items():
        out_deg[u] += len(neighbors_list)
        for v in neighbors_list:
            in_deg[v] += 1
            nodes.add(v)

    # Start from the node with out - in == 1 (Eulerian path start)
    start = next(
        (n for n in nodes if out_deg[n] - in_deg[n] == 1),
        next(iter(nodes)),  # fallback: any node (Eulerian circuit)
    )

    tmp: dict[str, list[str]] = {u: list(vs) for u, vs in graph.items()}
    stack = [start]
    path: list[str] = []

    while stack:
        u = stack[-1]
        if tmp.get(u):
            stack.append(tmp[u].pop())
        else:
            path.append(stack.pop())

    return path[::-1]


def path_to_genome(path: list[str]) -> str:
    """Reconstruct genome string from a node path in the De Bruijn graph."""
    return path[0] + "".join(node[-1] for node in path[1:])


def paired_path_to_genome(path: list[str], k: int, d: int) -> str:
    """
    Reconstruct a genome from a paired De Bruijn path.
    Raises ValueError if the two string reconstructions are inconsistent.
    """
    first  = path[0].split("|")[0] + "".join(n.split("|")[0][-1] for n in path[1:])
    second = path[0].split("|")[1] + "".join(n.split("|")[1][-1] for n in path[1:])

    overlap_start = k + d
    if first[overlap_start:] != second[: len(first) - overlap_start]:
        raise ValueError(
            "Paired-end reconstruction failed: overlapping regions do not match."
        )
    return first + second[len(first) - overlap_start :]


# ---------------------------------------------------------------------------
# CHAPTER 4 — Protein / mass spectrometry
# ---------------------------------------------------------------------------

# Integer monoisotopic masses (Da) for the 20 standard amino acids
MASS_TABLE: dict[str, int] = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186,
}

# Unique mass values used during de-novo sequencing
AMINO_ACID_MASSES: list[int] = sorted(set(MASS_TABLE.values()))


def get_linear_spectrum(peptide_masses: list[int]) -> list[int]:
    """Theoretical *linear* spectrum for a peptide given as a list of residue masses."""
    prefix = [0]
    for m in peptide_masses:
        prefix.append(prefix[-1] + m)
    spectrum = [0]
    for i in range(len(prefix)):
        for j in range(i + 1, len(prefix)):
            spectrum.append(prefix[j] - prefix[i])
    return sorted(spectrum)


def get_circular_spectrum(peptide_masses: list[int]) -> list[int]:
    """Theoretical *circular* (cyclic) spectrum for a peptide."""
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
    linear = get_linear_spectrum(peptide_masses)
    remaining = list(experimental_spectrum)
    for mass in linear:
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
        # Branch: extend each candidate by one amino acid
        candidates = [c + [m] for c in candidates for m in AMINO_ACID_MASSES]

        survivors: list[list[int]] = []
        for peptide in candidates:
            current_mass = sum(peptide)
            if current_mass == parent_mass:
                if get_circular_spectrum(peptide) == experimental_spectrum:
                    results.append(peptide)
                # Do NOT keep: already at parent mass, no further extension useful
            elif current_mass < parent_mass:
                if _is_consistent(peptide, experimental_spectrum):
                    survivors.append(peptide)
            # current_mass > parent_mass → pruned (not added to survivors)

        candidates = survivors

    return results
