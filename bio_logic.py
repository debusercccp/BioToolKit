import random

# ==========================================
# --- CAPITOLO 1: STRINGHE E SKEW (OriC) ---
# ==========================================

def hamming_distance(p: str, q: str) -> int:
    """Calcola la distanza di Hamming tra due stringhe."""
    return sum(1 for a, b in zip(p, q) if a != b)

def reverse_complement(pattern: str) -> str:
    """Genera il filamento complementare inverso."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(pattern))

def compute_skew(genome: str) -> list:
    """Calcola il G-C Skew per identificare l'origine di replicazione."""
    skew = [0]
    for base in genome:
        if base == 'G':
            skew.append(skew[-1] + 1)
        elif base == 'C':
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew

def neighbors(pattern, d):
    """Genera ricorsivamente tutti i k-mer con distanza di Hamming <= d."""
    if d == 0: return {pattern}
    if len(pattern) == 1: return {'A', 'C', 'G', 'T'}
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        if hamming_distance(pattern[1:], text) < d:
            for base in ['A', 'C', 'G', 'T']:
                neighborhood.add(base + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood

def frequent_words_with_mismatches_and_rc(text: str, k: int, d: int) -> set:
    """Trova i k-mer più frequenti considerando mismatch e RC."""
    counts = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        for sequence in [kmer, reverse_complement(kmer)]:
            for neighbor in neighbors(sequence, d):
                counts[neighbor] = counts.get(neighbor, 0) + 1
    if not counts: return set()
    max_count = max(counts.values())
    return {p for p, c in counts.items() if c == max_count}

# ==========================================
# --- CAPITOLO 2: PROFILI E MOTIF (Clock) --
# ==========================================

def create_profile_with_pseudocounts(motifs: list) -> dict:
    """Crea matrice di probabilità con Regola di Laplace (+1)."""
    t, k = len(motifs), len(motifs[0])
    profile = {base: [1.0] * k for base in "ACGT"}
    for motif in motifs:
        for i, base in enumerate(motif):
            profile[base][i] += 1
    for base in "ACGT":
        for i in range(k):
            profile[base][i] /= (t + 4)
    return profile

def get_kmer_probability(kmer: str, profile: dict) -> float:
    """Calcola probabilità di un k-mer dato un profilo."""
    prob = 1.0
    for i, base in enumerate(kmer):
        prob *= profile[base][i]
    return prob

def find_most_probable_kmer(text: str, k: int, profile: dict) -> str:
    """Trova il k-mer che meglio si adatta al profilo in una sequenza."""
    max_prob, best_kmer = -1.0, text[0:k]
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = get_kmer_probability(kmer, profile)
        if prob > max_prob:
            max_prob, best_kmer = prob, kmer
    return best_kmer

def get_score(motifs: list) -> int:
    """Calcola lo Score (mismatch totali rispetto al consenso)."""
    k = len(motifs[0])
    score = 0
    for i in range(k):
        column = [m[i] for m in motifs]
        max_freq = max(column.count(b) for b in "ACGT")
        score += (len(motifs) - max_freq)
    return score

def randomized_motif_search_instance(dna_list: list, k: int, t: int) -> list:
    """Singola corsa dell'algoritmo Randomized Motif Search."""
    motifs = []
    for seq in dna_list:
        start = random.randint(0, len(seq) - k)
        motifs.append(seq[start:start+k])
    best_motifs = motifs
    while True:
        profile = create_profile_with_pseudocounts(motifs)
        motifs = [find_most_probable_kmer(seq, k, profile) for seq in dna_list]
        if get_score(motifs) < get_score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs
