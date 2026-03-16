import bio_logic as bio
import matplotlib.pyplot as plt

def run_cap1():
    print("\n--- [CAP 1] ORIGIN FINDER PRO ---")
    genome = input("DNA (sequenza): ").strip().upper()
    genome = "".join(b for b in genome if b in "ACGT")
    
    skew_data = bio.compute_skew(genome)
    min_pos = [i for i, v in enumerate(skew_data) if v == min(skew_data)][-1]
    
    print(f"-> Analisi Skew completata. Minimo trovato a: {min_pos}")
    
    # Estrazione finestra (500bp)
    start = max(0, min_pos - 250)
    end = min(len(genome), min_pos + 250)
    window = genome[start:end]
    
    k = int(input("Lunghezza k-mer (es. 9): ") or 9)
    d = int(input("Mismatch (es. 1): ") or 1)
    
    boxes = bio.frequent_words_with_mismatches_and_rc(window, k, d)
    print(f"\n DnaA Boxes candidate: {boxes}")

    plt.plot(skew_data)
    plt.axvline(x=min_pos, color='r', linestyle='--')
    plt.title("G-C Skew Analysis")
    plt.show()

def run_cap2():
    print("\n--- [CAP 2] MOTIF FINDER (RANDOMIZED) ---")
    dna_list = []
    print("Inserisci sequenze (Invio vuoto per finire):")
    while True:
        line = input().strip().upper()
        if not line: break
        dna_list.append("".join(b for b in line if b in "ACGT"))
    
    if not dna_list: return
    
    k = int(input("Lunghezza motivo (k): ") or 15)
    best_motifs = bio.randomized_motif_search_instance(dna_list, k, len(dna_list))
    
    # Esecuzione iterata per precisione
    for _ in range(100):
        current = bio.randomized_motif_search_instance(dna_list, k, len(dna_list))
        if bio.get_score(current) < bio.get_score(best_motifs):
            best_motifs = current
            
    print("\n Migliori motivi trovati:")
    for m in best_motifs: print(m)
    print(f"Score: {bio.get_score(best_motifs)}")

def main():
    while True:
        print("\n" + "="*30 + "\n  BIOTOOLKIT \n" + "="*30)
        print("1) Ricerca OriC (Cap 1)")
        print("2) Ricerca Motivi (Cap 2)")
        print("q) Esci")
        choice = input("\nScelta: ").lower()
        if choice == '1': run_cap1()
        elif choice == '2': run_cap2()
        elif choice == 'q': break

if __name__ == "__main__":
    main()
