{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  name = "biotoolkit-env";

  # Pacchetti necessari per far girare il toolkit
  buildInputs = [
    (pkgs.python3.withPackages (ps: with ps; [
      numpy        # Necessario per Cap 7-8 (Matrici di distanza e Clustering)
      matplotlib   # Consigliato per visualizzare lo Skew del DNA (Cap 1)
    ]))
  ];

  # Variabili d'ambiente e comandi all'avvio
  shellHook = ''
    echo "---   BioToolKit Development Environment ---"
    echo "Python version: $(python --version)"
    echo "Dependencies: numpy, matplotlib loaded."
    echo "Digitare 'python main.py' per avviare il toolkit."
    echo "---------------------------------------------"
  '';
}
