# BioToolkit

Implementazione didattica degli algoritmi del libro
**Bioinformatics Algorithms** di Compeau & Pevzner (capitoli 1–11).

```
python3 main.py
```

Dipendenze: `numpy` (cap 7–8).  Tutto il resto usa solo la standard library.

---

## Struttura

```
bio_logic.py   — algoritmi puri, nessun I/O
main.py        — menu interattivo, un handler per capitolo
shell.nix      — Configurazione dell'ambiente di sviluppo riproducibile.
```

---

## Utilizzo

Utilizzo Standard (Linux/macOS/Windows)

Assicurati di avere numpy installato:

```
pip install numpy
python3 main.py
```

Utilizzo con Nix (NixOS o sistemi con Nix)

Il repository include uno shell.nix per caricare un ambiente isolato con tutte le dipendenze (Python 3.12, NumPy, Matplotlib):

```
nix-shell
python main.py
```

## Capitolo 1 — Origin Search (`run_cap1`)

Individua l'origine di replicazione (OriC) di un genoma batterico.

| Funzione | Descrizione |
|---|---|
| `hamming_distance(p, q)` | Numero di posizioni che differiscono tra due stringhe di uguale lunghezza |
| `reverse_complement(pattern)` | Complemento inverso Watson-Crick |
| `compute_skew(genome)` | Curva #G − #C lungo il genoma; il minimo indica l'OriC |
| `neighbors(pattern, d)` | Tutti i k-mer a distanza di Hamming ≤ d da `pattern` |
| `frequent_words_with_mismatches_and_rc(text, k, d)` | k-mer più frequenti con mismatch e reverse complement |

**Flusso nel menu**: inserisci la sequenza → calcola skew → cerca DnaA box nella finestra attorno al minimo → mostra il grafico (richiede matplotlib).

---

## Capitolo 2 — Motif Search (`run_cap2`)

Trova motivi conservati in un insieme di sequenze DNA (es. siti di legame per fattori di trascrizione).

| Funzione | Descrizione |
|---|---|
| `create_profile_with_pseudocounts(motifs)` | Profilo di frequenza con pseudoconti di Laplace (+1) |
| `get_kmer_probability(kmer, profile)` | Probabilità di un k-mer sotto un profilo |
| `find_most_probable_kmer(text, k, profile)` | k-mer più probabile in una sequenza data il profilo |
| `get_score(motifs)` | Score = somma dei mismatch dalla consensus; minore è meglio |
| `randomized_motif_search(dna_list, k, iterations)` | Ricerca randomizzata ripetuta `iterations` volte, restituisce il miglior risultato |

---

## Capitolo 3 — Genome Assembly (`run_cap3`)

Ricostruzione del genoma da reads tramite grafi di De Bruijn e cammini Euleriani.

| Funzione | Descrizione |
|---|---|
| `build_de_bruijn_from_kmers(patterns)` | Grafo di De Bruijn da una lista di k-mer |
| `build_paired_de_bruijn(paired_kmers)` | Grafo di De Bruijn per paired-end reads |
| `find_eulerian_path(graph)` | Algoritmo di Hierholzer per cammino Euleriano |
| `path_to_genome(path)` | Ricostruisce il genoma da un cammino nel grafo |
| `paired_path_to_genome(path, k, d)` | Ricostruzione da cammino paired-end con gap `d` |

**Sotto-menu**: 1) assembly standard da k-mer, 2) paired-end assembly.

---

## Capitolo 4 — Protein Tools (`run_cap4`)

Strumenti per la spettrometria di massa e il sequenziamento di peptidi.

| Funzione | Descrizione |
|---|---|
| `MASS_TABLE` | Masse monoisotopiche intere dei 20 amminoacidi standard |
| `get_linear_spectrum(peptide_masses)` | Spettro teorico lineare di un peptide |
| `get_circular_spectrum(peptide_masses)` | Spettro teorico ciclico (per antibiotici) |
| `cyclopeptide_sequencing(experimental_spectrum)` | Branch-and-bound: trova tutti i peptidi ciclici che corrispondono allo spettro |

**Sotto-menu**: 1) peptide → spettro (lineare o ciclico), 2) spettro → peptide (sequenziamento).

---

## Capitolo 5 — Sequence Alignment (`run_cap5`)

Allineamento a coppie di sequenze DNA o proteiche.

| Funzione | Descrizione |
|---|---|
| `global_alignment(s1, s2, matrix, gap)` | Needleman-Wunsch con matrice di sostituzione opzionale |
| `local_alignment(s1, s2, matrix, gap)` | Smith-Waterman con matrice di sostituzione opzionale |
| `BLOSUM62` | Matrice log-odds per sequenze con > 30% identità |
| `PAM250` | Matrice evolutiva per sequenze distanti (< 25% identità) |
| `IDENTITY` | Matrice piatta +1/−1 per DNA o uso generico |
| `SUBSTITUTION_MATRICES` | Dict `{"BLOSUM62": ..., "PAM250": ..., "IDENTITY": ...}` |

Il menu chiede modalità (global/local), matrice e gap penalty prima di allineare.

---

## Capitolo 6 — Genome Rearrangements (`run_cap6`)

Distanza evolutiva tra genomi tramite inversioni.

| Funzione | Descrizione |
|---|---|
| `greedy_sorting(p)` | Ordina una permutazione con segni tramite inversioni greedy; restituisce i passi intermedi |
| `count_breakpoints(p)` | Conta i breakpoint; `breakpoints / 2` è il lower bound sulla reversal distance reale |

**Input**: permutazione con segni come `+1 -3 +2` (il `+` è opzionale, si accettano anche virgole).

---

## Capitolo 7 — Molecular Evolution (`run_cap7`)

Filogenesi e ricostruzione degli antenati.

| Funzione | Descrizione |
|---|---|
| `upgma(distance_matrix, labels)` | UPGMA: albero radicato ultrametrico da matrice di distanze |
| `neighbor_joining(distance_matrix, labels)` | Neighbor Joining: albero non radicato additivo (Saitou & Nei 1987) |
| `small_parsimony(tree, leaf_chars)` | Algoritmo di Fitch su albero fissato (singolo sito) |
| `parsimony_score(tree, leaf_sequences)` | Parsimonia su sequenze allineate, tutti i siti |
| `newick(tree, labels)` | Serializzazione in formato Newick (compatibile FigTree, iTOL) |

**Sotto-menu**: 1) UPGMA, 2) NJ (entrambi richiedono matrice di distanze), 3) small parsimony (inserisci sequenze, albero costruito via NJ automaticamente).

---

## Capitolo 8 — Gene Clustering (`run_cap8`)

Clustering di profili di espressione genica.

| Funzione | Descrizione |
|---|---|
| `euclidean_distance(v, w)` | Distanza euclidea tra due vettori |
| `lloyd_algorithm(data, k, iterations)` | Hard K-Means; cluster vuoti mantengono il centro precedente |
| `soft_kmeans(data, k, beta, iterations)` | Soft K-Means EM con parametro di rigidità β; β alto → comportamento hard |
| `hierarchical_clustering(data, linkage)` | Clustering agglomerativo; `linkage` ∈ `{"single","complete","average"}` |

**Input**: punti come righe di float separati da spazio (es. `2.1 0.5`), una per riga, riga vuota per terminare.

---

## Capitolo 9 — Read Mapping / BWT (`run_cap9`)

Indicizzazione e mapping di reads su un genoma di riferimento.

| Funzione | Descrizione |
|---|---|
| `suffix_array(text)` | Suffix array via prefix doubling O(n log² n) |
| `bwt_transform(text)` | BWT via suffix array (O(n log² n) tempo, O(n) spazio — nessuna rotazione esplicita) |
| `bwt_inverse(bwt)` | Inversione BWT via LF-mapping con rank table precomputata O(n) |
| `bwt_match(text, pattern)` | Backward search esatta O(p · \|Σ\|) con occurrence table |
| `bwt_match_with_mismatches(text, pattern, max_mismatches)` | Ricerca inesatta via pigeonhole + verifica Hamming |

**Differenza dalla versione naïve**: `bwt_transform` non costruisce le n rotazioni esplicite (O(n²) memoria); `bwt_inverse` usa LF-mapping O(n) invece di `list.index()` O(n²).

---

## Capitolo 10 — HMM & Gene Finding (`run_cap10`)

Hidden Markov Models per gene prediction e analisi di sequenze.

| Funzione | Descrizione |
|---|---|
| `viterbi(obs, states, start_p, trans_p, emit_p)` | Cammino di stati più probabile (log-space, backpointer efficiente) |
| `forward(obs, states, start_p, trans_p, emit_p)` | log P(osservazioni \| modello) via forward algorithm |
| `backward(obs, states, trans_p, emit_p)` | Probabilità backward β[t][s] in log-space |
| `baum_welch(obs, states, alphabet, ...)` | Training EM dei parametri HMM; converge su `\|ΔlogL\| < tol` |
| `cpg_island_hmm()` | HMM pre-costruito 2-stati (CpG / non-CpG) da Durbin et al. 1998 |
| `find_cpg_islands(sequence, min_length)` | Predice isole CpG con Viterbi, filtra per lunghezza minima |

**HMM pre-costruiti nel menu**: dado carico/onesto (Viterbi, Forward, Baum-Welch); input custom disponibile per tutti gli algoritmi.

---

## Capitolo 11 — Advanced Alignment (`run_cap11`)

Allineamento con penalità di gap affini e allineamento multiplo.

| Funzione | Descrizione |
|---|---|
| `affine_global_alignment(s1, s2, matrix, gap_open, gap_extend)` | Needleman-Wunsch con penalità affine (3 matrici DP: M, X, Y) |
| `affine_local_alignment(s1, s2, matrix, gap_open, gap_extend)` | Smith-Waterman con penalità affine |
| `linear_space_alignment(s1, s2, ...)` | Hirschberg divide-and-conquer: stessa qualità di NW affine ma O(n) memoria |
| `multiple_sequence_alignment(sequences, ...)` | MSA via centre-star heuristic (approssimazione 2−2/k dell'ottimo) |

**Penalità affine**: `gap_open + gap_extend × L` per un gap di lunghezza L.  Default: `gap_open=-11, gap_extend=-1` (valori tipici BLAST con BLOSUM62).

**Differenza cap 5 vs cap 11**: il cap 5 usa penalità lineare (ogni gap costa la stessa cosa); il cap 11 usa penalità affine (aprire un gap costa molto, estenderlo costa poco) — molto più biologicamente realistico.

---

## Note tecniche

- Tutti gli algoritmi sono in `bio_logic.py`, completamente separati dall'I/O.
- Gli handler in `main.py` usano `ValueError` / `KeyError` per la gestione degli errori; il loop principale li cattura e stampa il messaggio senza terminare il programma.
- Le matrici di sostituzione (BLOSUM62, PAM250) sono hardcoded come `dict[tuple[str,str], int]`; non richiedono file esterni.
- Il formato Newick in output (cap 7) è compatibile con FigTree, iTOL, e qualsiasi parser standard.
