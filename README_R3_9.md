# Molecular descriptor generation for Zn adsorption screening

`generate_additive_descriptors.py` converts SMILES into the original 35 descriptor columns and exports the final 13 separately in the order listed in Table S4. The script performs the descriptor-generation step of the workflow. Exact adsorption-energy predictions additionally require the original fitting-set imputer and fitted ExtraTrees forest described in Section 4.5.

## Setup

Use Python 3.12 in a virtual environment and install the RDKit release used in the study.

```sh
python -m pip install rdkit==2026.9.1
```

For Excel input, also run `python -m pip install openpyxl==3.1.5`.

## Run the seven-compound example

```sh
python generate_additive_descriptors.py --example --out-prefix seven
```

The example uses the exact input SMILES in Table S5. In the verified run, the final 13 columns reproduced the 81 available values and 10 unavailable entries used in the seven-compound assessment.

## Run a new list

Give the file one `SMILES` column and an optional `Name` column. Use a header row and one structure per row. Preserve the input SMILES for an exact repeat of a conformer calculation. The declared stereochemistry, charge, protonation and counterions define the represented species.

```sh
python generate_additive_descriptors.py --input additives.csv --out-prefix additives
python generate_additive_descriptors.py --input additives.xlsx --sheet Additives --out-prefix additives_xlsx
```

Use `--sheet` to select a named worksheet. The first worksheet is the default. Choose a fresh output prefix for each run.

| Output | Contents |
| --- | --- |
| `<prefix>_descriptors35.csv` | Input identities, the 35 descriptor columns in their original order, and structure, charge, conformer and missing-value flags. |
| `<prefix>_model13.csv` | Input identities followed by the final 13 raw descriptor columns in model order. |

## Descriptor conventions

- RDKit graph descriptors use the parsed SMILES. Charge descriptors and conformer generation use explicit hydrogens. Canonical SMILES are recorded for identity checks.
- The 3D procedure uses one ETKDGv3 conformer with seed 42, random starting coordinates and one thread. It uses MMFF94 when fully parameterized and UFF when MMFF94 parameters are unavailable and UFF is supported. The optimization must converge within 500 iterations. Volume uses a 0.2 Å grid and 2.0 Å margin. Inertia and asphericity use atomic masses.
- The procedure skips 3D calculations for disconnected species, species without heavy atoms, structures with more than 80 heavy atoms, and inputs containing Li, Na, K, Mg, Ca, Zn, Ag, Al, Fe, Cu, Mo, W or Te. Invalid SMILES, unsupported charge or conformer calculations and undefined ratios are flagged. Blank cells represent unavailable values.
- The stored chalcogen counts equal O + S in the development set. The additional-compound recipe also counts Se, Te, Po and Lv. These heavier chalcogens were absent from the fitting set.
- This procedure reproduces the final 13 inputs used for the additional-compound assessment. The historical development analysis used the stored descriptor matrix.

If applying the original model, use the 13 numerical columns with the original fitting-set median imputer and fitted forest. The original ExtraTrees pipeline used no scaling. Preserve the input flags when interpreting predictions for molecules with unavailable descriptors.
