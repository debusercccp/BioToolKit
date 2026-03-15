use std::collections::{HashMap, HashSet};
use std::io::{self, Write};

// --- 1. MAPPATURA E BASI ---

fn pattern_to_num(pattern: &str) -> u64 {
    let base_map = |c| match c {
        'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3,
        _ => 0,
    };
    let mut number = 0;
    let k = pattern.len() as u32;
    for (i, c) in pattern.chars().enumerate() {
        number += base_map(c) * 4u64.pow(k - 1 - i as u32);
    }
    number
}

fn number_to_pattern(mut index: u64, k: usize) -> String {
    let num_to_base = ['A', 'C', 'G', 'T'];
    let mut pattern = String::with_capacity(k);
    for _ in 0..k {
        pattern.push(num_to_base[(index % 4) as usize]);
        index /= 4;
    }
    pattern.chars().rev().collect() // Rust costruisce al contrario, quindi invertiamo
}

fn hamming_distance(p: &str, q: &str) -> usize {
    p.chars().zip(q.chars()).filter(|(a, b)| a != b).count()
}

fn reverse_complement(pattern: &str) -> String {
    pattern.chars().rev().map(|base| {
        match base {
            'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A',
            _ => base,
        }
    }).collect()
}

// --- 2. ANALISI DELLO SKEW ---

fn compute_skew(genome: &str) -> Vec<i32> {
    let mut skew = Vec::with_capacity(genome.len() + 1);
    skew.push(0);
    let mut current = 0;
    for base in genome.chars() {
        match base {
            'G' => current += 1,
            'C' => current -= 1,
            _ => (),
        }
        skew.push(current);
    }
    skew
}

// --- 3. GENERAZIONE VICINATO ---

fn neighbors(pattern: &str, d: usize) -> HashSet<String> {
    if d == 0 {
        let mut hs = HashSet::new();
        hs.insert(pattern.to_string());
        return hs;
    }
    if pattern.len() == 1 {
        return vec!["A", "C", "G", "T"].into_iter().map(|s| s.to_string()).collect();
    }

    let mut neighborhood = HashSet::new();
    let suffix = &pattern[1..];
    let first_char = pattern.chars().next().unwrap();
    let suffix_neighbors = neighbors(suffix, d);

    for text in suffix_neighbors {
        if hamming_distance(&pattern[1..], &text) < d {
            for base in ["A", "C", "G", "T"] {
                neighborhood.insert(format!("{}{}", base, text));
            }
        } else {
            neighborhood.insert(format!("{}{}", first_char, text));
        }
    }
    neighborhood
}

// --- 4. RICERCA DNAA BOXES ---

fn frequent_words_with_mismatches_and_rc(text: &str, k: usize, d: usize) -> HashSet<String> {
    let mut counts: HashMap<String, u32> = HashMap::new();

    for i in 0..=(text.len() - k) {
        let kmer = &text[i..i+k];
        
        // Unione dei vicini del kmer corrente e del suo RC
        let mut combined_vicini = neighbors(kmer, d);
        combined_vicini.extend(neighbors(&reverse_complement(kmer), d));

        for neighbor in combined_vicini {
            *counts.entry(neighbor).or_insert(0) += 1;
        }
    }

    if counts.is_empty() { return HashSet::new(); }
    let max_count = *counts.values().max().unwrap();
    counts.into_iter()
        .filter(|(_, count)| *count == max_count)
        .map(|(pattern, _)| pattern)
        .collect()
}

// --- 5. MAIN INTERACTOR ---

fn main() {
    println!("\n=============================================");
    println!("      DNA TOOLKIT: ORIGIN FINDER (RUST PRO)  ");
    println!("=============================================");

    print!("DNA (sequenza): ");
    io::stdout().flush().unwrap();
    let mut input_dna = String::new();
    io::stdin().read_line(&mut input_dna).unwrap();
    
    let raw_genome: String = input_dna.trim().to_uppercase()
        .chars().filter(|&c| "ACGT".contains(c)).collect();

    if raw_genome.is_empty() {
        println!("ERRORE: Sequenza non valida!");
        return;
    }

    // [1/3] SKEW
    let skew_data = compute_skew(&raw_genome);
    let min_val = *skew_data.iter().min().unwrap();
    let all_min_pos: Vec<usize> = skew_data.iter().enumerate()
        .filter(|&(_, &val)| val == min_val)
        .map(|(i, _)| i).collect();
    
    let pos = *all_min_pos.last().unwrap();
    println!("\n[1/3] Minimo Skew: {} trovato a posizione {}", min_val, pos);

    // [2/3] FINESTRA
    let window_size = 500;
    let start = if pos > window_size / 2 { pos - window_size / 2 } else { 0 };
    let end = std::cmp::min(raw_genome.len(), pos + window_size / 2);
    let window = &raw_genome[start..end];
    println!("[2/3] Analisi finestra {}:{}", start, end);

    // [3/3] PARAMETRI K-MER
    let k_len = 9; // Puoi aggiungere un input per renderlo dinamico
    let max_d = 1;

    println!("[3/3] Ricerca {}-mer con d={}...", k_len, max_d);
    let candidati = frequent_words_with_mismatches_and_rc(window, k_len, max_d);

    println!("\n---------------------------------------------");
    println!("RISULTATI FINALI:");
    println!("Punto di inversione: {}", pos);
    println!("DnaA Boxes probabili: {:?}", candidati);
    println!("---------------------------------------------");
}
