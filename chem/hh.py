#!/usr/bin/env python3
"""
Henderson-Hasselbalch esteso
  - Database pKa (monprotici e poliprotici)
  - Calcolo forward  : pH da concentrazioni, con scelta dello step
  - Calcolo inverso  : rapporto o concentrazione mancante da pH noto
  - Speciazione      : frazioni α_i per ogni specie a pH dato
"""

import argparse
import math

# ── Database ───────────────────────────────────────────────────────────────────

PKA_DATABASE = {
    # Monprotici
    "acido_acetico":    {"pKa": [4.76],              "formula": "CH₃COOH",
                         "species": ["CH₃COOH",  "CH₃COO⁻"]},
    "acido_formico":    {"pKa": [3.75],              "formula": "HCOOH",
                         "species": ["HCOOH",    "HCOO⁻"]},
    "acido_lattico":    {"pKa": [3.86],              "formula": "C₃H₆O₃",
                         "species": ["C₃H₅O₃H", "C₃H₅O₃⁻"]},
    "ammoniaca":        {"pKa": [9.25],              "formula": "NH₄⁺",
                         "species": ["NH₄⁺",     "NH₃"]},
    # Diprotici
    "acido_carbonico":  {"pKa": [6.35, 10.33],       "formula": "H₂CO₃",
                         "species": ["H₂CO₃",   "HCO₃⁻",   "CO₃²⁻"]},
    "acido_ossalico":   {"pKa": [1.25, 4.27],        "formula": "C₂H₂O₄",
                         "species": ["H₂Ox",     "HOx⁻",    "Ox²⁻"]},
    "acido_maleico":    {"pKa": [1.92, 6.23],        "formula": "C₄H₄O₄",
                         "species": ["H₂Mal",    "HMal⁻",   "Mal²⁻"]},
    "acido_succinico":  {"pKa": [4.21, 5.64],        "formula": "C₄H₆O₄",
                         "species": ["H₂Suc",    "HSuc⁻",   "Suc²⁻"]},
    "glicina":          {"pKa": [2.34, 9.60],        "formula": "Gly",
                         "species": ["H₂Gly⁺",  "HGly",    "Gly⁻"]},
    "alanina":          {"pKa": [2.34, 9.69],        "formula": "Ala",
                         "species": ["H₂Ala⁺",  "HAla",    "Ala⁻"]},
    # Triprotici
    "acido_citrico":    {"pKa": [3.13, 4.76, 6.40],  "formula": "C₆H₈O₇",
                         "species": ["H₃Cit",    "H₂Cit⁻",  "HCit²⁻",  "Cit³⁻"]},
    "acido_fosforico":  {"pKa": [2.15, 7.20, 12.35], "formula": "H₃PO₄",
                         "species": ["H₃PO₄",   "H₂PO₄⁻",  "HPO₄²⁻",  "PO₄³⁻"]},
    "istidina":         {"pKa": [1.82, 6.00, 9.17],  "formula": "His",
                         "species": ["H₃His²⁺", "H₂His⁺",  "HHis",    "His⁻"]},
    "lisina":           {"pKa": [2.18, 8.95, 10.53], "formula": "Lys",
                         "species": ["H₃Lys²⁺", "H₂Lys⁺",  "HLys",    "Lys⁻"]},
    "acido_aspartico":  {"pKa": [1.88, 3.65, 9.60],  "formula": "Asp",
                         "species": ["H₃Asp⁺",  "H₂Asp",   "HAsp⁻",   "Asp²⁻"]},
    "acido_glutammico": {"pKa": [2.19, 4.25, 9.67],  "formula": "Glu",
                         "species": ["H₃Glu⁺",  "H₂Glu",   "HGlu⁻",   "Glu²⁻"]},
}


# ── Core math ──────────────────────────────────────────────────────────────────

def hh_forward(pKa: float, Bconc: float, Aconc: float) -> float:
    """pH = pKa + log10([A⁻] / [HA])"""
    if Aconc <= 0:
        raise ValueError("La concentrazione dell'acido deve essere > 0.")
    if Bconc <= 0:
        raise ValueError("La concentrazione della base deve essere > 0.")
    return pKa + math.log10(Bconc / Aconc)


def hh_inverse_ratio(pH: float, pKa: float) -> float:
    """Dato pH e pKa → rapporto [A⁻]/[HA]"""
    return 10 ** (pH - pKa)


def hh_inverse_conc(pH: float, pKa: float,
                    known_conc: float, known_is_base: bool) -> float:
    """
    Dato pH, pKa e una concentrazione nota → l'altra concentrazione.
      known_is_base=True  : conosco [A⁻], voglio [HA]
      known_is_base=False : conosco [HA], voglio [A⁻]
    """
    ratio = hh_inverse_ratio(pH, pKa)   # [A⁻]/[HA]
    if known_is_base:
        return known_conc / ratio        # [HA] = [A⁻] / ratio
    else:
        return known_conc * ratio        # [A⁻] = [HA] * ratio


def speciation(pH: float, pKa_list: list) -> list:
    """
    Frazioni di speciazione α_i per un acido n-protico.

    Derivazione:
      Per H_n A  →  n+1 specie: H_n A, H_{n-1}A⁻, …, A^{n-}
      α_i = (h^{n-i} · ∏_{j<i} Ka_j) / D
      D   = Σ_{i=0}^{n} (h^{n-i} · ∏_{j<i} Ka_j)
    """
    h  = 10.0 ** (-pH)
    n  = len(pKa_list)
    Ka = [10.0 ** (-pk) for pk in pKa_list]

    terms, prod = [], 1.0
    for i in range(n + 1):
        terms.append((h ** (n - i)) * prod)
        if i < n:
            prod *= Ka[i]

    D = sum(terms)
    return [t / D for t in terms]


def suggest_step(pH: float, pKa_list: list) -> int:
    """Restituisce l'indice (0-based) del pKa più vicino al pH dato."""
    return min(range(len(pKa_list)), key=lambda i: abs(pKa_list[i] - pH))


# ── Output helpers ─────────────────────────────────────────────────────────────

SEP = "  " + "─" * 50

def print_speciation(pH: float, compound: dict, name: str):
    fracs   = speciation(pH, compound["pKa"])
    species = compound["species"]
    dom_idx = fracs.index(max(fracs))

    print(f"\n  Speciazione di {name}  (formula: {compound['formula']})  a pH {pH:.2f}")
    print(SEP)
    print(f"  {'Specie':<16} {'α':>10}  {'%':>8}")
    print(SEP)
    for sp, fr in zip(species, fracs):
        marker = " ◀ dominante" if sp == species[dom_idx] else ""
        print(f"  {sp:<16} {fr:>10.6f}  {fr*100:>7.2f}%{marker}")
    print(SEP)


def print_forward_result(pH: float, pKa: float, step_idx: int,
                         pair: list, compound: dict | None):
    print(f"\n  Step {step_idx + 1}: pKa = {pKa}")
    print(f"  Equilibrio : {pair[0]} ⇌ {pair[1]} + H⁺")
    print(f"  pH calcolato: {pH:.4f}")
    if compound and len(compound["pKa"]) > 1:
        print_speciation(pH, compound, "")


def list_compounds():
    print(f"\n  {'Nome':<22} {'pKa (ogni step)':^34}  Formula")
    print("  " + "─" * 70)
    for name, data in PKA_DATABASE.items():
        pka_str = " / ".join(f"{p:.2f}" for p in data["pKa"])
        print(f"  {name:<22} {pka_str:^34}  {data['formula']}")
    print()


# ── Argparse ───────────────────────────────────────────────────────────────────

def getArgs():
    parser = argparse.ArgumentParser(
        description="Henderson-Hasselbalch esteso: poliprotici, inverso, speciazione",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Esempi:
  # Lista composti
  %(prog)s -l

  # Forward monprotico (manuale)
  %(prog)s -pKa 4.76 -bc 0.1 -ac 0.05

  # Forward poliprotico step 2 (fosfato)
  %(prog)s -c acido_fosforico -s 2 -bc 0.08 -ac 0.12

  # Speciazione a pH 7.4
  %(prog)s -c acido_fosforico --pH 7.4 --speciate

  # Inverso: rapporto [A-]/[HA] dato pH
  %(prog)s -c acido_acetico --pH 5.0 --inv-ratio

  # Inverso: calcola [HA] noto [A-]=0.1 M e pH=7.4 (step 2 fosfato)
  %(prog)s -c acido_fosforico -s 2 --pH 7.4 --inv-base 0.1
""",
    )

    parser.add_argument("-l", "--list", action="store_true",
                        help="Elenca i composti nel database")

    # Sorgente pKa
    src = parser.add_mutually_exclusive_group()
    src.add_argument("-c", "--compound", metavar="NOME",
                     help="Composto dal database")
    src.add_argument("-pKa", "--pka", type=float, metavar="FLOAT",
                     help="pKa manuale (monprotico)")

    parser.add_argument("-s", "--step", type=int, default=1, metavar="N",
                        help="Step di dissociazione per poliprotici (1-based, default=1)")

    # Forward
    fwd = parser.add_argument_group("Calcolo forward (pH da concentrazioni)")
    fwd.add_argument("-bc", "--baseC", type=float, metavar="FLOAT",
                     help="Concentrazione base coniugata [A⁻]  (M)")
    fwd.add_argument("-ac", "--acidC", type=float, metavar="FLOAT",
                     help="Concentrazione acido [HA]           (M)")

    # Inverso
    inv = parser.add_argument_group("Calcolo inverso")
    inv.add_argument("--pH", type=float, metavar="FLOAT",
                     help="pH target")
    inv.add_argument("--inv-ratio", action="store_true",
                     help="Calcola [A⁻]/[HA] dato pH e pKa")
    inv.add_argument("--inv-base", type=float, metavar="FLOAT",
                     help="Noto [A⁻] (M): calcola [HA] mancante")
    inv.add_argument("--inv-acid", type=float, metavar="FLOAT",
                     help="Noto [HA] (M): calcola [A⁻] mancante")

    # Speciazione
    parser.add_argument("--speciate", action="store_true",
                        help="Calcola speciazione a --pH dato (richiede -c e --pH)")

    return parser.parse_args()


def resolve_pKa(args):
    """Restituisce (pKa_list, compound_dict | None)."""
    if args.compound:
        key = args.compound.lower()
        if key not in PKA_DATABASE:
            raise ValueError(f"Composto '{key}' non trovato. Usa -l per la lista.")
        return PKA_DATABASE[key]["pKa"], PKA_DATABASE[key]
    if args.pka is not None:
        return [args.pka], None
    return None, None


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = getArgs()

    if args.list:
        list_compounds()
        return

    try:
        pKa_list, compound = resolve_pKa(args)
    except ValueError as e:
        print(f"Errore: {e}")
        return

    # ── Speciazione pura ─────────────────────────────────────────────────────
    if args.speciate:
        if compound is None:
            print("Errore: --speciate richiede -c (composto dal database).")
            return
        if args.pH is None:
            print("Errore: --speciate richiede --pH.")
            return
        print_speciation(args.pH, compound, args.compound)
        return

    # ── Calcolo inverso ──────────────────────────────────────────────────────
    inverse_requested = args.inv_ratio or \
                        args.inv_base is not None or \
                        args.inv_acid is not None

    if inverse_requested:
        if args.pH is None or pKa_list is None:
            print("Errore: il calcolo inverso richiede --pH e (-c o -pKa).")
            return

        step_idx = args.step - 1
        if step_idx >= len(pKa_list):
            print(f"Errore: step {args.step} non esiste "
                  f"(questo composto ha {len(pKa_list)} step).")
            return

        pKa = pKa_list[step_idx]
        sp  = compound["species"] if compound else ["HA", "A⁻"]

        print(f"\n  Composto : {args.compound or 'manuale'}")
        print(f"  Step     : {args.step}   pKa = {pKa}")
        print(f"  pH target: {args.pH}")
        print(SEP)

        if args.inv_ratio:
            ratio = hh_inverse_ratio(args.pH, pKa)
            print(f"  [{sp[step_idx+1]}] / [{sp[step_idx]}] = {ratio:.6f}  ({ratio:.4e})")
            pct_base = ratio / (1 + ratio) * 100
            print(f"  → {pct_base:.1f}% in forma deprotonata ({sp[step_idx+1]})")

        elif args.inv_base is not None:
            acid = hh_inverse_conc(args.pH, pKa, args.inv_base, known_is_base=True)
            print(f"  [{sp[step_idx+1]}] noto  = {args.inv_base:.6f} M")
            print(f"  [{sp[step_idx]}]  calcolato = {acid:.6f} M")

        elif args.inv_acid is not None:
            base = hh_inverse_conc(args.pH, pKa, args.inv_acid, known_is_base=False)
            print(f"  [{sp[step_idx]}]  noto  = {args.inv_acid:.6f} M")
            print(f"  [{sp[step_idx+1]}] calcolato = {base:.6f} M")

        print(SEP)

        # Speciazione bonus se poliprotico
        if compound and len(pKa_list) > 1:
            print_speciation(args.pH, compound, args.compound)
        return

    # ── Calcolo forward ──────────────────────────────────────────────────────
    if args.baseC is None or args.acidC is None or pKa_list is None:
        print("Errore: specifica -bc e -ac (e -c o -pKa) per il calcolo forward.")
        return

    step_idx = args.step - 1
    if step_idx >= len(pKa_list):
        print(f"Errore: step {args.step} non esiste "
              f"(questo composto ha {len(pKa_list)} step).")
        return

    pKa = pKa_list[step_idx]
    sp  = compound["species"] if compound else ["HA", "A⁻"]

    try:
        pH = hh_forward(pKa, args.baseC, args.acidC)
    except ValueError as e:
        print(f"Errore: {e}")
        return

    print_forward_result(pH, pKa, step_idx, [sp[step_idx], sp[step_idx + 1]], compound)


if __name__ == "__main__":
    main()
