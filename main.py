import bio_logic as bio
import matplotlib.pyplot as plt

def run_cap1():
    print("\n--- [CAP 1] ORIGIN FINDER PRO ---")
    genome = input("DNA sequence: ").strip().upper()
    genome = "".join(b for b in genome if b in "ACGT")
    skew_data = bio.compute_skew(genome)
    min_pos = [i for i, v in enumerate(skew_data) if v == min(skew_data)][-1]
    print(f"Skew analysis complete. Minimum at: {min_pos}")
    window = genome[max(0, min_pos - 250):min(len(genome), min_pos + 250)]
    k = int(input("k-mer length (default 9): ") or 9)
    d = int(input("Max mismatches (default 1): ") or 1)
    boxes = bio.frequent_words_with_mismatches_and_rc(window, k, d)
    print(f"Candidate DnaA Boxes: {boxes}")
    plt.plot(skew_data)
    plt.axvline(x=min_pos, color='r', linestyle='--')
    plt.show()

def run_cap2():
    print("\n--- [CAP 2] MOTIF FINDER ---")
    dna_list = []
    print("Enter sequences (Empty line to stop):")
    while True:
        line = input().strip().upper()
        if not line: break
        dna_list.append("".join(b for b in line if b in "ACGT"))
    if not dna_list: return
    k = int(input("Motif length (k): ") or 15)
    best_motifs = bio.randomized_motif_search_instance(dna_list, k, len(dna_list))
    for _ in range(100):
        current = bio.randomized_motif_search_instance(dna_list, k, len(dna_list))
        if bio.get_score(current) < bio.get_score(best_motifs):
            best_motifs = current
    print("\nBest motifs found:")
    for m in best_motifs: print(m)
    print(f"Score: {bio.get_score(best_motifs)}")

def run_cap3():
    print("\n--- [CAP 3.1] STANDARD ASSEMBLY ---")
    kmers = []
    print("Enter k-mers (Empty line to stop):")
    while True:
        line = input().strip().upper()
        if not line: break
        kmers.append(line)
    if not kmers: return
    try:
        graph = bio.build_de_bruijn_from_kmers(kmers)
        path = bio.find_eulerian_path(graph)
        result = bio.path_to_genome(path)
        print(f"Assembled genome: {result}")
    except Exception as e:
        print(f"Assembly error: {e}")

def run_cap3_paired():
    print("\n--- [CAP 3.2] PAIRED-END ASSEMBLY ---")
    k = int(input("k-mer length: ") or 3)
    d = int(input("Gap distance (d): ") or 1)
    pairs = []
    print("Enter pairs (format k1|k2):")
    while True:
        line = input().strip().upper()
        if not line: break
        if '|' in line: pairs.append(line.split('|'))
    if not pairs: return
    try:
        graph = bio.build_paired_de_bruijn(pairs)
        path = bio.find_eulerian_path(graph)
        result = bio.paired_path_to_genome(path, k, d)
        print(f"Assembled genome: {result}")
    except Exception as e:
        print(f"Assembly error: {e}")

def main():
    while True:
        print("\n" + "="*30 + "\n  BIOTOOLKIT MENU\n" + "="*30)
        print("1) Origin Search (Cap 1)")
        print("2) Motif Search (Cap 2)")
        print("3) Standard Assembly (Cap 3.1)")
        print("4) Paired-End Assembly (Cap 3.2)")
        print("q) Quit")
        choice = input("\nChoice: ").lower()
        if choice == '1': run_cap1()
        elif choice == '2': run_cap2()
        elif choice == '3': run_cap3()
        elif choice == '4': run_cap3_paired()
        elif choice == 'q': break

if __name__ == "__main__":
    main()
