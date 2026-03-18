import random
import collections

# ==========================================
# --- CAPITOLO 1: STRINGHE E SKEW (OriC) ---
# ==========================================

def hamming_distance(p: str, q: str) -> int:
    return sum(1 for a, b in zip(p, q) if a != b)

def reverse_complement(pattern: str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(pattern))

def compute_skew(genome: str) -> list:
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
    prob = 1.0
    for i, base in enumerate(kmer):
        prob *= profile[base][i]
    return prob

def find_most_probable_kmer(text: str, k: int, profile: dict) -> str:
    max_prob, best_kmer = -1.0, text[0:k]
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = get_kmer_probability(kmer, profile)
        if prob > max_prob:
            max_prob, best_kmer = prob, kmer
    return best_kmer

def get_score(motifs: list) -> int:
    k = len(motifs[0])
    score = 0
    for i in range(k):
        column = [m[i] for m in motifs]
        max_freq = max(column.count(b) for b in "ACGT")
        score += (len(motifs) - max_freq)
    return score

def randomized_motif_search_instance(dna_list: list, k: int, t: int) -> list:
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

# ==========================================
# --- CAPITOLO 3: GENOME ASSEMBLY (Grafi) --
# ==========================================

def build_de_bruijn_from_kmers(patterns: list) -> dict:
    graph = {}
    for pattern in patterns:
        prefix, suffix = pattern[:-1], pattern[1:]
        if prefix not in graph: graph[prefix] = []
        graph[prefix].append(suffix)
    return graph

def build_paired_de_bruijn(paired_kmers: list) -> dict:
    graph = {}
    for pair in paired_kmers:
        u = f"{pair[0][:-1]}|{pair[1][:-1]}"
        v = f"{pair[0][1:]}|{pair[1][1:]}"
        if u not in graph: graph[u] = []
        graph[u].append(v)
    return graph

def find_eulerian_path(graph: dict) -> list:
    in_degree = collections.defaultdict(int)
    out_degree = collections.defaultdict(int)
    nodes = set(graph.keys())
    for u in graph:
        out_degree[u] = len(graph[u])
        for v in graph[u]:
            in_degree[v] += 1
            nodes.add(v)
    start_node = next(iter(nodes))
    for node in nodes:
        if out_degree[node] > in_degree[node]:
            start_node = node
            break
    stack, path = [start_node], []
    temp_graph = {u: list(v) for u, v in graph.items()}
    while stack:
        u = stack[-1]
        if u in temp_graph and temp_graph[u]:
            stack.append(temp_graph[u].pop())
        else:
            path.append(stack.pop())
    return path[::-1]

def path_to_genome(path: list) -> str:
    genome = path[0]
    for i in range(1, len(path)):
        genome += path[i][-1]
    return genome

def paired_path_to_genome(path: list, k: int, d: int) -> str:
    first_string = path[0].split('|')[0]
    second_string = path[0].split('|')[1]
    for i in range(1, len(path)):
        first_string += path[i].split('|')[0][-1]
        second_string += path[i].split('|')[1][-1]
    overlap_len = len(first_string) - (k + d)
    if first_string[k+d:] == second_string[:overlap_len]:
        return first_string + second_string[overlap_len:]
    return "Error: String reconstruction failed (non-matching overlaps)"
