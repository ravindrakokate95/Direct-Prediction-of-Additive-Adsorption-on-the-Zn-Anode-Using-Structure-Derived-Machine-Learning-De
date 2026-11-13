#!/usr/bin/env python3
"""Generate structure-derived descriptors for Zn adsorption screening.

Use Python 3.12 and install rdkit==2026.9.1. For Excel input, also install
openpyxl==3.1.5. Input needs a SMILES column. Name is optional.

Examples:
  python generate_additive_descriptors.py --example --out-prefix seven
  python generate_additive_descriptors.py --input additives.csv --out-prefix additives
  python generate_additive_descriptors.py --input additives.xlsx --sheet Additives --out-prefix additives

Each run writes <prefix>_descriptors35.csv and <prefix>_model13.csv.
The 13 model columns follow Table S4. Blank cells are unavailable values.
Conformer settings follow Section 4.2 and the seven-compound assessment.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


RDKIT_VERSION = "2026.09.1"
FEATURES35 = (
    "MinPartialCharge", "NumHAcceptors", "PartialCharge_range", "NHOHCount",
    "TPSA_per_heavy_atom", "coordination_atom_fraction", "Halogen_fraction",
    "MolVolume_3D", "MaxPartialCharge", "MaxEStateIndex", "n_chalcogen_total",
    "MolLogP", "LogP_to_TPSA_ratio", "NumRotatableBonds", "BertzCT", "PMI1_3D",
    "PMI3_3D", "fr_Al_OH", "MeanAbsHeavyAtomPartialCharge", "n_halogen_total",
    "RotatableBond_fraction", "fr_NH2", "n_N", "MinAbsEStateIndex", "SPS", "n_Na",
    "fr_C_O", "NumHeteroatoms", "n_Cl", "InertialShapeFactor_3D", "MR_per_volume",
    "fr_COO", "Asphericity_3D", "n_F", "BalabanJ",
)
FEATURES13 = (
    "NHOHCount", "n_N", "fr_NH2", "MR_per_volume", "fr_Al_OH", "n_chalcogen_total",
    "PMI1_3D", "MaxPartialCharge", "SPS", "InertialShapeFactor_3D",
    "MinAbsEStateIndex", "fr_C_O", "n_F",
)
GRAPH_DESCRIPTORS = (
    "NumHAcceptors", "NHOHCount", "MaxEStateIndex", "MolLogP", "NumRotatableBonds",
    "BertzCT", "fr_Al_OH", "fr_NH2", "MinAbsEStateIndex", "SPS", "fr_C_O",
    "NumHeteroatoms", "fr_COO", "BalabanJ",
)
MODEL_VOLUME_AND_INERTIA = ("MolVolume_3D", "PMI1_3D", "InertialShapeFactor_3D", "MR_per_volume")
UNSUPPORTED_3D_ELEMENTS = frozenset({
    "Li", "Na", "K", "Mg", "Ca", "Zn", "Ag", "Al", "Fe", "Cu", "Mo", "W", "Te",
})
DIAGNOSTIC_COLUMNS = (
    "Canonical_SMILES", "Structure_status", "Charge_status", "Conformer_status",
    "Force_field", "Missing_descriptors35", "Missing_final_descriptors", "RDKit_version",
)
EXAMPLE = (
    ("4-hydroxybenzoic acid", "O=C(O)c1ccc(O)cc1"),
    ("cetyltrimethylammonium bromide", "CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]"),
    ("diethylene glycol dimethyl ether", "COCCOCCOC"),
    ("dimethyl malonate", "COC(=O)CC(=O)OC"),
    ("methyl(1-phenylethyl)tellane", "C[Te]C(C)c1ccccc1"),
    ("N-vinyl-2-pyrrolidone", "C=CN1CCCC1=O"),
    ("sodium saccharinate", "[Na+].O=C1[N-]S(=O)(=O)c2ccccc12"),
)


def descriptor_value(function, molecule, **kwargs):
    try:
        value = float(function(molecule, **kwargs))
    except (RuntimeError, ValueError, ArithmeticError):
        return None
    return value if math.isfinite(value) else None


def ratio(numerator, denominator):
    if numerator is None or denominator is None or denominator == 0:
        return None
    return numerator / denominator


def charge_descriptors(molecule_with_h, result):
    # Use the same RDKit extrema as the fitted-model assessment.
    high = descriptor_value(Descriptors.MaxPartialCharge, molecule_with_h)
    low = descriptor_value(Descriptors.MinPartialCharge, molecule_with_h)
    charges = [
        float(atom.GetProp("_GasteigerCharge")) if atom.HasProp("_GasteigerCharge") else math.nan
        for atom in molecule_with_h.GetAtoms()
    ]
    if not any(map(math.isfinite, charges)):
        result["Charge_status"] = "unavailable"
        return

    result["MaxPartialCharge"] = high
    result["MinPartialCharge"] = low
    if high is not None and low is not None:
        result["PartialCharge_range"] = high - low

    heavy_charges = [
        charge for atom, charge in zip(molecule_with_h.GetAtoms(), charges)
        if atom.GetAtomicNum() > 1
    ]
    if heavy_charges and all(map(math.isfinite, heavy_charges)):
        result["MeanAbsHeavyAtomPartialCharge"] = sum(map(abs, heavy_charges)) / len(heavy_charges)
    result["Charge_status"] = "ok" if all(map(math.isfinite, charges)) else "incomplete"


def conformer_descriptors(molecule, molecule_with_h, result):
    if result["Structure_status"] == "disconnected_species":
        result["Conformer_status"] = "disconnected_species"
        return
    heavy_atoms = molecule.GetNumHeavyAtoms()
    if heavy_atoms == 0:
        result["Conformer_status"] = "no_heavy_atoms"
        return
    if heavy_atoms > 80:
        result["Conformer_status"] = "over_80_heavy_atoms"
        return
    if any(atom.GetSymbol() in UNSUPPORTED_3D_ELEMENTS for atom in molecule.GetAtoms()):
        result["Conformer_status"] = "excluded_element"
        return

    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 1
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(molecule_with_h, params) != 0:
            result["Conformer_status"] = "embedding_failed"
            return

        if AllChem.MMFFHasAllMoleculeParams(molecule_with_h):
            result["Force_field"] = "MMFF94"
            code = AllChem.MMFFOptimizeMolecule(molecule_with_h, mmffVariant="MMFF94", maxIters=500)
        elif AllChem.UFFHasAllMoleculeParams(molecule_with_h):
            result["Force_field"] = "UFF"
            code = AllChem.UFFOptimizeMolecule(molecule_with_h, maxIters=500)
        else:
            result["Conformer_status"] = "force_field_unavailable"
            return
        if code != 0:
            result["Conformer_status"] = "optimization_not_converged"
            return

        volume = descriptor_value(AllChem.ComputeMolVolume, molecule_with_h, gridSpacing=0.2, boxMargin=2.0)
        if volume is None or volume <= 0:
            result["Conformer_status"] = "volume_unavailable"
            return
        result["MolVolume_3D"] = volume
        result["PMI1_3D"] = descriptor_value(rdMolDescriptors.CalcPMI1, molecule_with_h, useAtomicMasses=True)
        result["InertialShapeFactor_3D"] = descriptor_value(
            rdMolDescriptors.CalcInertialShapeFactor, molecule_with_h, useAtomicMasses=True
        )
        result["MR_per_volume"] = ratio(descriptor_value(Descriptors.MolMR, molecule), volume)
        if any(result[key] is None for key in MODEL_VOLUME_AND_INERTIA):
            for key in MODEL_VOLUME_AND_INERTIA:
                result[key] = None
            result["Conformer_status"] = "volume_or_inertia_unavailable"
            return

        result["PMI3_3D"] = descriptor_value(rdMolDescriptors.CalcPMI3, molecule_with_h, useAtomicMasses=True)
        result["Asphericity_3D"] = descriptor_value(
            rdMolDescriptors.CalcAsphericity, molecule_with_h, useAtomicMasses=True
        )
        result["Conformer_status"] = (
            "ok" if result["PMI3_3D"] is not None and result["Asphericity_3D"] is not None
            else "additional_3D_descriptor_unavailable"
        )
    except (RuntimeError, ValueError, ArithmeticError):
        result["Conformer_status"] = "conformer_calculation_failed"


def calculate_descriptors(name, smiles):
    result = dict.fromkeys(FEATURES35)
    result.update(
        Name=name, SMILES=smiles, Canonical_SMILES="", Structure_status="ok",
        Charge_status="not_calculated", Conformer_status="not_generated", Force_field="",
        RDKit_version=rdBase.rdkitVersion,
    )
    params = Chem.SmilesParserParams()
    params.parseName = False
    molecule = Chem.MolFromSmiles(smiles, params) if smiles else None
    if molecule is None or not molecule.GetNumAtoms():
        result["Structure_status"] = "invalid_SMILES" if smiles else "missing_SMILES"
        return result

    # Canonical SMILES are metadata. Calculations use the original atom order.
    result["Canonical_SMILES"] = Chem.MolToSmiles(molecule, isomericSmiles=True)
    if len(Chem.GetMolFrags(molecule)) != 1:
        result["Structure_status"] = "disconnected_species"

    for key in GRAPH_DESCRIPTORS:
        result[key] = descriptor_value(getattr(Descriptors, key), molecule)

    atom_numbers = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    heavy_atoms = molecule.GetNumHeavyAtoms()
    tpsa = descriptor_value(Descriptors.TPSA, molecule)
    result.update(
        n_N=atom_numbers.count(7), n_F=atom_numbers.count(9),
        n_Cl=atom_numbers.count(17), n_Na=atom_numbers.count(11),
        n_halogen_total=sum(n in (9, 17, 35, 53) for n in atom_numbers),
        n_chalcogen_total=sum(n in (8, 16, 34, 52, 84, 116) for n in atom_numbers),
        TPSA_per_heavy_atom=ratio(tpsa, heavy_atoms),
        LogP_to_TPSA_ratio=ratio(result["MolLogP"], tpsa),
        coordination_atom_fraction=ratio(sum(n in (7, 8, 15, 16) for n in atom_numbers), heavy_atoms),
        RotatableBond_fraction=ratio(result["NumRotatableBonds"], heavy_atoms),
    )
    result["Halogen_fraction"] = ratio(result["n_halogen_total"], heavy_atoms)

    try:
        molecule_with_h = Chem.AddHs(molecule)
    except (RuntimeError, ValueError):
        result["Charge_status"] = "hydrogen_addition_failed"
        result["Conformer_status"] = "hydrogen_addition_failed"
        return result
    charge_descriptors(molecule_with_h, result)
    conformer_descriptors(molecule, molecule_with_h, result)
    return result


def parse_rows(rows):
    header = [str(value).strip() if value is not None else "" for value in next(rows, ())]
    if header.count("SMILES") != 1 or header.count("Name") > 1:
        raise ValueError("Use one SMILES column and, optionally, one Name column in the first row.")
    smiles_index = header.index("SMILES")
    name_index = header.index("Name") if "Name" in header else None

    def cell(row, index):
        if index is None or index >= len(row) or row[index] is None:
            return ""
        return str(row[index]).strip()

    for line, row in enumerate(rows, start=2):
        if not any(value is not None and str(value).strip() for value in row):
            continue
        yield cell(row, name_index) or f"molecule_{line - 1}", cell(row, smiles_index)


def read_additives(path, sheet):
    if path.suffix.lower() == ".csv":
        if sheet is not None:
            raise ValueError("--sheet is for .xlsx input.")
        with path.open(encoding="utf-8-sig", newline="") as stream:
            return list(parse_rows(csv.reader(stream)))
    if path.suffix.lower() == ".xlsx":
        try:
            from openpyxl import load_workbook
        except ImportError as exc:
            raise ValueError("For Excel input, install openpyxl==3.1.5.") from exc
        workbook = load_workbook(path, read_only=True, data_only=True)
        try:
            worksheet = workbook[sheet] if sheet is not None else workbook.worksheets[0]
            return list(parse_rows(worksheet.values))
        finally:
            workbook.close()
    raise ValueError("Use a .csv or .xlsx input file.")


def write_outputs(prefix, rows):
    columns35 = ("Name", "SMILES", *FEATURES35, *DIAGNOSTIC_COLUMNS)
    columns13 = ("Name", "SMILES", *FEATURES13)
    paths = (Path(f"{prefix}_descriptors35.csv"), Path(f"{prefix}_model13.csv"))
    for path in paths:
        if path.exists():
            raise FileExistsError(f"{path} already exists. Choose another --out-prefix.")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    created = []
    try:
        for path, columns in zip(paths, (columns35, columns13)):
            with path.open("x", encoding="utf-8", newline="") as stream:
                created.append(path)
                writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
    except (OSError, KeyboardInterrupt):
        for path in created:
            path.unlink()
        raise
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input", type=Path, help="CSV or Excel file with a SMILES column.")
    inputs.add_argument("--example", action="store_true", help="Use the seven compounds in Table S5.")
    parser.add_argument("--sheet", help="Excel worksheet name. Defaults to the first worksheet.")
    parser.add_argument("--out-prefix", required=True, type=Path, help="Path prefix for the two CSVs.")
    args = parser.parse_args()

    if rdBase.rdkitVersion != RDKIT_VERSION:
        parser.error(f"Use RDKit {RDKIT_VERSION} (pip install rdkit==2026.9.1). Found {rdBase.rdkitVersion}.")
    if args.example and args.sheet is not None:
        parser.error("--sheet is for .xlsx input.")
    try:
        additives = EXAMPLE if args.example else read_additives(args.input, args.sheet)
        if not additives:
            raise ValueError("The input contains no molecule rows.")
        rows = []
        for name, smiles in additives:
            result = calculate_descriptors(name, smiles)
            result["Missing_descriptors35"] = ", ".join(key for key in FEATURES35 if result[key] is None)
            result["Missing_final_descriptors"] = ", ".join(key for key in FEATURES13 if result[key] is None)
            rows.append(result)
        paths = write_outputs(args.out_prefix, rows)
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        parser.error(str(exc))
    incomplete = sum(bool(row["Missing_final_descriptors"]) for row in rows)
    print(f"{len(rows)} input rows processed. {incomplete} have unavailable final descriptors.")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
