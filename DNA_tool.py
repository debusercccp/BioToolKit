import matplotlib.pyplot as plt

# --- 1. FUNZIONI DI BASE & MAPPATURA ---

def hammingDistance(p: str, q: str) -> int:
    """Calcola quante basi differiscono tra due stringhe."""
    return sum(1 for a, b in zip(p, q) if a != b)

def reverseComplement(pattern: str) -> str:
    """Genera il filamento complementare inverso (fondamentale per DnaA boxes)."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(pattern))

# --- 2. ANALISI DELLO SKEW (Ricerca macro-regione OriC) ---

def compute_skew(genome: str) -> list:
    """Calcola il G-C Skew per identificare l'inversione dei filamenti."""
    skew = [0]
    for base in genome:
        if base == 'G':
            skew.append(skew[-1] + 1)
        elif base == 'C':
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew

# --- 3. GENERAZIONE VICINATO (L'anima della ricerca con mismatch) ---

def neighbors(pattern, d):
    """Genera ricorsivamente tutti i k-mer con distanza di Hamming <= d."""
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        if hammingDistance(pattern[1:], text) < d:
            for base in ['A', 'C', 'G', 'T']:
                neighborhood.add(base + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood

# --- 4. RICERCA DNAA BOXES (Versione Ottimizzata con Dizionario) ---

def frequentWordsWithMismatchesAndRC(text: str, k: int, d: int) -> set:
    """Trova i k-mer più frequenti considerando mismatch e filamento opposto."""
    counts = {} 
    
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        # Cerchiamo sia il k-mer che il suo Reverse Complement
        for sequence in [kmer, reverseComplement(kmer)]:
            vicini = neighbors(sequence, d)
            for neighbor in vicini:
                counts[neighbor] = counts.get(neighbor, 0) + 1
            
    if not counts:
        return set()

    max_count = max(counts.values())
    return {pattern for pattern, count in counts.items() if count == max_count}

# --- 5. MAIN INTERACTOR ---

if __name__ == "__main__":
    print("\n" + "="*45)
    print("      DNA TOOLKIT: ORIGIN FINDER PRO        ")
    print("="*45)
    
    # Input e pulizia
    raw_input = input("DNA (sequenza o genoma): ").strip().upper()
    raw_genome = "".join(base for base in raw_input if base in "ACGT")
    
    if not raw_genome:
        print("ERRORE: Sequenza non valida!")
        exit()

    genome_id = input("Nome della sequenza: ") or "Genoma Ignoto"
    
    # [1/3] SKEW ANALYSIS
    print("\n[1/3] Analisi Skew Diagram...")
    skew_data = compute_skew(raw_genome)
    min_val = min(skew_data)
    
    # Trova TUTTI i minimi per non limitarsi all'indice 0
    all_min_pos = [i for i, val in enumerate(skew_data) if val == min_val]
    
    # Spesso l'OriC è verso la fine della zona di calo, prendiamo l'ultimo dei minimi
    pos = all_min_pos[-1] 
    print(f"-> Minimo globale ({min_val}) trovato in {len(all_min_pos)} punti.")
    print(f"-> Candidato OriC scelto: posizione {pos}")
    
    # [2/3] WINDOW DEFINITION
    print("\n[2/3] Estrazione finestra intorno al minimo...")
    window_size = 500
    start = max(0, pos - (window_size // 2))
    end = min(len(raw_genome), pos + (window_size // 2))
    window = raw_genome[start:end]
    print(f"-> Analisi finestra [{start}:{end}] (lunghezza: {len(window)} bp)")
    
    # [3/3] DNAA BOX SEARCH
    k_len = int(input("\nLunghezza k-mer (es. 9): ") or 9)
    max_d = int(input("Mismatch massimi (es. 1): ") or 1)
    
    print(f"\n[3/3] Ricerca dei {k_len}-mer più frequenti con d={max_d}...")
    candidati = frequentWordsWithMismatchesAndRC(window, k_len, max_d)
    
    print("\n" + "—"*40)
    print(f"RISULTATI FINALI PER: {genome_id}")
    print(f"—"*40)
    print(f"Punto di inversione Skew: {pos}")
    print(f"Potenziali DnaA Boxes: {candidati}")
    print("—"*40)

    # Plotting
    plt.figure(figsize=(10, 5))
    plt.plot(skew_data, label='Skew Value', color='blue', linewidth=1.5)
    plt.scatter(pos, min_val, color='red', zorder=5, label='OriC Candidate')
    plt.axvline(x=pos, color='red', linestyle='--', alpha=0.5)
    plt.title(f"G-C Skew Diagram - {genome_id}")
    plt.xlabel("Posizione nel genoma")
    plt.ylabel("Skew (G-C)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
