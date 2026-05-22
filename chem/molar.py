#!/usr/bin/env python3
"""
Calcolatore di Molarità, Moli e Massa
  - Calcolo Moli     : n = m / PM
  - Calcolo Molarità : M = n / V
  - Calcolo inverso  : PM da (g, n) | Massa da (PM, n) | Volume da (M, n)
  - Database         : tavola periodica completa per ricerca PM per simbolo
  - Parser formula   : calcola PM da formula chimica (es. H2SO4, NaCl, C6H12O6)
"""

import argparse
import re
import sys

# ── Tavola periodica (simbolo → peso atomico IUPAC 2021) ──────────────────────

ELEMENTI: dict[str, float] = {
    "H":  1.008,   "He": 4.0026,  "Li": 6.941,   "Be": 9.0122,
    "B":  10.811,  "C":  12.011,  "N":  14.007,  "O":  15.999,
    "F":  18.998,  "Ne": 20.180,  "Na": 22.990,  "Mg": 24.305,
    "Al": 26.982,  "Si": 28.086,  "P":  30.974,  "S":  32.06,
    "Cl": 35.45,   "Ar": 39.948,  "K":  39.098,  "Ca": 40.078,
    "Sc": 44.956,  "Ti": 47.867,  "V":  50.942,  "Cr": 51.996,
    "Mn": 54.938,  "Fe": 55.845,  "Co": 58.933,  "Ni": 58.693,
    "Cu": 63.546,  "Zn": 65.38,   "Ga": 69.723,  "Ge": 72.630,
    "As": 74.922,  "Se": 78.971,  "Br": 79.904,  "Kr": 83.798,
    "Rb": 85.468,  "Sr": 87.62,   "Y":  88.906,  "Zr": 91.224,
    "Nb": 92.906,  "Mo": 95.95,   "Tc": 98.0,    "Ru": 101.07,
    "Rh": 102.91,  "Pd": 106.42,  "Ag": 107.87,  "Cd": 112.41,
    "In": 114.82,  "Sn": 118.71,  "Sb": 121.76,  "Te": 127.60,
    "I":  126.90,  "Xe": 131.29,  "Cs": 132.91,  "Ba": 137.33,
    "La": 138.91,  "Ce": 140.12,  "Pr": 140.91,  "Nd": 144.24,
    "Pm": 145.0,   "Sm": 150.36,  "Eu": 151.96,  "Gd": 157.25,
    "Tb": 158.93,  "Dy": 162.50,  "Ho": 164.93,  "Er": 167.26,
    "Tm": 168.93,  "Yb": 173.04,  "Lu": 174.97,  "Hf": 178.49,
    "Ta": 180.95,  "W":  183.84,  "Re": 186.21,  "Os": 190.23,
    "Ir": 192.22,  "Pt": 195.08,  "Au": 196.97,  "Hg": 200.59,
    "Tl": 204.38,  "Pb": 207.2,   "Bi": 208.98,  "Po": 209.0,
    "At": 210.0,   "Rn": 222.0,   "Fr": 223.0,   "Ra": 226.0,
    "Ac": 227.0,   "Th": 232.04,  "Pa": 231.04,  "U":  238.03,
    "Np": 237.0,   "Pu": 244.0,   "Am": 243.0,   "Cm": 247.0,
    "Bk": 247.0,   "Cf": 251.0,   "Es": 252.0,   "Fm": 257.0,
    "Md": 258.0,   "No": 259.0,   "Lr": 262.0,   "Rf": 267.0,
    "Db": 270.0,   "Sg": 271.0,   "Bh": 270.0,   "Hs": 277.0,
    "Mt": 276.0,   "Ds": 281.0,   "Rg": 282.0,   "Cn": 285.0,
    "Nh": 286.0,   "Fl": 289.0,   "Mc": 290.0,   "Lv": 293.0,
    "Ts": 294.0,   "Og": 294.0,
}

# Composti comuni con nome per la ricerca rapida con -cn
COMPOSTI_COMUNI: dict[str, str] = {
    "acqua":                 "H2O",
    "cloruro_di_sodio":      "NaCl",
    "glucosio":              "C6H12O6",
    "saccarosio":            "C12H22O11",
    "acido_solforico":       "H2SO4",
    "acido_cloridrico":      "HCl",
    "acido_nitrico":         "HNO3",
    "acido_acetico":         "C2H4O2",
    "idrossido_di_sodio":    "NaOH",
    "idrossido_di_potassio": "KOH",
    "carbonato_di_sodio":    "Na2CO3",
    "bicarbonato_di_sodio":  "NaHCO3",
    "solfato_di_rame":       "CuSO4",
    "cloruro_di_calcio":     "CaCl2",
    "acido_fosforico":       "H3PO4",
    "ammoniaca":             "NH3",
    "etanolo":               "C2H6O",
    "metanolo":              "CH4O",
    "acetone":               "C3H6O",
    "urea":                  "CH4N2O",
    "glicina":               "C2H5NO2",
    "alanina":               "C3H7NO2",
    "acido_citrico":         "C6H8O7",
    "acido_lattico":         "C3H6O3",
    "EDTA":                  "C10H16N2O8",
    "tris":                  "C4H11NO3",
}


# ── Parser formula chimica ─────────────────────────────────────────────────────

def _parse_formula(formula: str, pos: int = 0) -> tuple[dict, int]:
    """
    Parser ricorsivo per formule chimiche, supporta:
      - Elementi singoli e multi-carattere : NaCl, Fe2O3
      - Coefficienti numerici              : H2SO4, C6H12O6
      - Gruppi con parentesi               : Ca(OH)2, Al2(SO4)3
      - Parentesi quadre                   : [Fe(CN)6]3-  (ignora la carica)
    Restituisce un dict {simbolo: conteggio} e la posizione finale.
    """
    counts: dict[str, float] = {}

    def add(d: dict, sym: str, n: float):
        d[sym] = d.get(sym, 0) + n

    while pos < len(formula):
        c = formula[pos]

        if c in ("(", "["):
            # Gruppo: ricorsione fino alla parentesi chiusa corrispondente
            close = ")" if c == "(" else "]"
            sub, pos = _parse_formula(formula, pos + 1)
            if pos >= len(formula) or formula[pos] != close:
                raise ValueError(f"Parentesi non chiusa nella formula '{formula}'.")
            pos += 1  # salta la parentesi chiusa
            # Coefficiente dopo la parentesi
            num_start = pos
            while pos < len(formula) and formula[pos].isdigit():
                pos += 1
            multiplier = int(formula[num_start:pos]) if pos > num_start else 1
            for sym, cnt in sub.items():
                add(counts, sym, cnt * multiplier)

        elif c.isupper():
            # Simbolo elemento: una maiuscola + eventuali minuscole
            sym_start = pos
            pos += 1
            while pos < len(formula) and formula[pos].islower():
                pos += 1
            sym = formula[sym_start:pos]
            if sym not in ELEMENTI:
                raise ValueError(f"Elemento sconosciuto '{sym}' nella formula.")
            # Coefficiente dopo il simbolo
            num_start = pos
            while pos < len(formula) and formula[pos].isdigit():
                pos += 1
            n = int(formula[num_start:pos]) if pos > num_start else 1
            add(counts, sym, n)

        elif c in (")", "]"):
            # Fine del gruppo corrente (gestito dal livello superiore)
            break

        elif c in ("+", "-") or c.isdigit():
            # Carica ionica o numero isolato: ignorati a livello di radice
            pos += 1

        elif c == ".":
            # Idrati: es. CuSO4·5H2O → continua a parsare
            pos += 1

        else:
            raise ValueError(f"Carattere inatteso '{c}' nella formula '{formula}'.")

    return counts, pos


def pm_da_formula(formula: str) -> tuple[float, dict]:
    """
    Calcola il Peso Molecolare da una formula chimica.
    Restituisce (PM, {simbolo: (conteggio, peso_atomico)}).
    """
    counts, _ = _parse_formula(formula)
    dettaglio: dict[str, tuple[int, float]] = {}
    pm = 0.0
    for sym, n in counts.items():
        pa = ELEMENTI[sym]
        pm += n * pa
        dettaglio[sym] = (int(n), pa)
    return pm, dettaglio


# ── Core math ──────────────────────────────────────────────────────────────────

def calcola_moli(massa: float, pm: float) -> float:
    """n = m / PM"""
    if pm <= 0:
        raise ValueError("Il Peso Molecolare deve essere > 0.")
    if massa < 0:
        raise ValueError("La massa non puo' essere negativa.")
    return massa / pm


def calcola_molarita(moli: float, volume: float) -> float:
    """M = n / V"""
    if volume <= 0:
        raise ValueError("Il volume deve essere > 0.")
    if moli < 0:
        raise ValueError("Le moli non possono essere negative.")
    return moli / volume


# ── Output helpers ─────────────────────────────────────────────────────────────

SEP = "  " + "-" * 52


def stampa_composizione(formula: str, dettaglio: dict, pm: float):
    print(f"\n  Formula        : {formula}")
    print(f"  {'Elemento':<10} {'Atomi':>6}  {'Peso atomico':>14}  {'Contributo':>12}")
    print(SEP)
    for sym, (n, pa) in dettaglio.items():
        contributo = n * pa
        print(f"  {sym:<10} {n:>6}  {pa:>12.4f} u  {contributo:>10.4f} u")
    print(SEP)
    print(f"  PM calcolato   : {pm:.4f} g/mol")


def stampa_lista_composti():
    print(f"\n  {'Nome':<28} {'Formula':<16}  PM (g/mol)")
    print("  " + "-" * 58)
    for nome, formula in COMPOSTI_COMUNI.items():
        try:
            pm, _ = pm_da_formula(formula)
            print(f"  {nome:<28} {formula:<16}  {pm:.4f}")
        except ValueError:
            print(f"  {nome:<28} {formula:<16}  (errore)")
    print()


def stampa_lista_elementi():
    print(f"\n  {'Simbolo':<6} {'Peso atomico (u)':>18}")
    print("  " + "-" * 28)
    for sym, pa in ELEMENTI.items():
        print(f"  {sym:<6} {pa:>18.4f}")
    print()


# ── Argparse ───────────────────────────────────────────────────────────────────

def getArgs():
    parser = argparse.ArgumentParser(
        description="Calcolatore chimico moli/molarita' con database tavola periodica.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Esempi:
  # PM da formula (parser automatico)
  %(prog)s -f H2SO4
  %(prog)s -f Ca(OH)2
  %(prog)s -f C6H12O6

  # Composto per nome (database interno)
  %(prog)s -cn glucosio

  # Moli da grammi e formula
  %(prog)s -g 5.0 -f NaCl

  # Molarita' completa: grammi + formula + volume
  %(prog)s -g 10.5 -f C6H12O6 -v 0.5

  # PM manuale (senza formula)
  %(prog)s -g 5.0 -pm 58.44 -v 0.25

  # Inverso: PM da grammi e moli
  %(prog)s -g 3.65 -n 0.0625

  # Inverso: massa da PM e moli
  %(prog)s -pm 180.16 -n 0.25

  # Inverso: volume da molarita' e moli
  %(prog)s -n 0.5 -M 2.0

  # Liste di riferimento
  %(prog)s --lista-composti
  %(prog)s --lista-elementi
""",
    )

    # Sorgente PM
    pm_group = parser.add_mutually_exclusive_group()
    pm_group.add_argument("-f", "--formula", metavar="FORMULA",
                          help="Formula chimica per calcolo PM automatico (es. NaCl, Ca(OH)2)")
    pm_group.add_argument("-cn", "--composto-nome", metavar="NOME",
                          help="Nome dal database composti comuni (usa --lista-composti)")
    pm_group.add_argument("-pm", "--pm", type=float, metavar="FLOAT",
                          help="Peso Molecolare manuale (g/mol)")

    # Grandezze chimiche
    parser.add_argument("-g", "--grammi", type=float, metavar="FLOAT",
                        help="Massa della sostanza (g)")
    parser.add_argument("-n", "--moli", type=float, metavar="FLOAT",
                        help="Quantita' di sostanza (mol)")
    parser.add_argument("-v", "--volume", type=float, metavar="FLOAT",
                        help="Volume della soluzione (L)")
    parser.add_argument("-M", "--molarita", type=float, metavar="FLOAT",
                        help="Molarita' nota (mol/L) per calcolo inverso del volume")

    # Liste
    parser.add_argument("--lista-composti", action="store_true",
                        help="Mostra il database dei composti comuni con i loro PM")
    parser.add_argument("--lista-elementi", action="store_true",
                        help="Mostra tutti gli elementi con i pesi atomici")

    return parser.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = getArgs()

    if args.lista_composti:
        stampa_lista_composti()
        return

    if args.lista_elementi:
        stampa_lista_elementi()
        return

    # Nessun argomento
    if not any([args.formula, args.composto_nome, args.pm,
                args.grammi, args.moli, args.volume, args.molarita]):
        print("Errore: specifica almeno un parametro. Usa -h per gli esempi.")
        return

    # ── Risoluzione PM ────────────────────────────────────────────────────────
    pm_val    = args.pm
    dettaglio = None
    formula   = None

    if args.formula or args.composto_nome:
        if args.composto_nome:
            key = args.composto_nome.lower()
            if key not in COMPOSTI_COMUNI:
                print(f"Errore: composto '{key}' non trovato. "
                      f"Usa --lista-composti per vedere le opzioni.")
                return
            formula = COMPOSTI_COMUNI[key]
            print(f"\n  Composto riconosciuto: {args.composto_nome} -> {formula}")
        else:
            formula = args.formula

        try:
            pm_val, dettaglio = pm_da_formula(formula)
        except ValueError as e:
            print(f"Errore nel parsing della formula: {e}")
            return

    # ── Stampa composizione se disponibile ───────────────────────────────────
    if dettaglio is not None:
        stampa_composizione(formula, dettaglio, pm_val)

    # ── Calcoli ───────────────────────────────────────────────────────────────
    print(f"\n  [ Calcoli ]")
    print(SEP)

    moli_eff = args.moli  # moli disponibili per il calcolo successivo

    # Caso 1: g + PM  ->  moli
    if args.grammi is not None and pm_val is not None and args.moli is None:
        try:
            moli_eff = calcola_moli(args.grammi, pm_val)
            print(f"  Massa             : {args.grammi:.4f} g")
            print(f"  PM                : {pm_val:.4f} g/mol")
            print(f"  -> Moli           : {moli_eff:.6f} mol")
        except ValueError as e:
            print(f"  Errore: {e}")
            return

    # Caso 2: g + moli  ->  PM inverso
    elif args.grammi is not None and args.moli is not None and pm_val is None:
        if args.moli <= 0:
            print("  Errore: le moli devono essere > 0 per calcolare il PM.")
            return
        pm_calc = args.grammi / args.moli
        print(f"  Massa             : {args.grammi:.4f} g")
        print(f"  Moli              : {args.moli:.6f} mol")
        print(f"  -> PM calcolato   : {pm_calc:.4f} g/mol")

    # Caso 3: PM + moli  ->  massa inversa
    elif pm_val is not None and args.moli is not None and args.grammi is None:
        g_calc = args.moli * pm_val
        print(f"  PM                : {pm_val:.4f} g/mol")
        print(f"  Moli              : {args.moli:.6f} mol")
        print(f"  -> Massa calcolata: {g_calc:.4f} g")

    # Caso 4: solo moli fornite a mano
    elif args.moli is not None and args.grammi is None and pm_val is None:
        print(f"  Moli              : {args.moli:.6f} mol")

    # Caso 5: solo PM (senza grammi o moli) - utile se si vuole solo il PM
    elif pm_val is not None and args.grammi is None and args.moli is None:
        if dettaglio is None:
            # PM fornito a mano, gia' stampato sopra come parte della composizione
            print(f"  PM                : {pm_val:.4f} g/mol")

    # ── Molarita' e calcoli da volume ─────────────────────────────────────────

    # Caso A: n + V  ->  Molarita'
    if args.volume is not None and moli_eff is not None:
        try:
            molarita = calcola_molarita(moli_eff, args.volume)
            print(f"  Volume            : {args.volume:.4f} L")
            print(SEP)
            print(f"  MOLARITA'         : {molarita:.6f} mol/L")
        except ValueError as e:
            print(f"  Errore: {e}")

    # Caso B: n + M  ->  Volume inverso
    elif args.molarita is not None and moli_eff is not None:
        if args.molarita <= 0:
            print("  Errore: la molarita' deve essere > 0.")
        else:
            vol_calc = moli_eff / args.molarita
            print(f"  Molarita' nota    : {args.molarita:.6f} mol/L")
            print(SEP)
            print(f"  -> Volume         : {vol_calc:.6f} L  ({vol_calc * 1000:.4f} mL)")

    # Caso C: V + M  ->  moli inverse
    elif args.volume is not None and args.molarita is not None and moli_eff is None:
        moli_calc = args.molarita * args.volume
        print(f"  Volume            : {args.volume:.4f} L")
        print(f"  Molarita'         : {args.molarita:.6f} mol/L")
        print(SEP)
        print(f"  -> Moli           : {moli_calc:.6f} mol")
        if pm_val is not None:
            g_calc = moli_calc * pm_val
            print(f"  -> Massa          : {g_calc:.4f} g")

    # Info: moli disponibili ma niente volume
    elif moli_eff is not None and args.volume is None and args.molarita is None:
        if args.grammi is not None or pm_val is not None:
            print(SEP)
            print("  Info: aggiungi -v (volume in L) o -M (molarita') per ulteriori calcoli.")

    print(SEP + "\n")


if __name__ == "__main__":
    main()
