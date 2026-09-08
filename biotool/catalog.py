"""Original educational descriptions of compatible, publicly available PDB entries."""

from dataclasses import dataclass


@dataclass(frozen=True)
class ProteinExample:
    pdb_id: str
    name: str
    organism: str
    description: str


EXAMPLES = (
    ProteinExample(
        "1CRN", "Crambin", "Crambe hispanica subsp. abyssinica",
        "A small plant protein: explore its compact shape and composition "
        "through 46 observed residues.",
    ),
    ProteinExample(
        "1UBQ", "Human ubiquitin", "Homo sapiens",
        "Helps tag proteins for degradation. "
        "A compact example with 76 observed residues.",
    ),
    ProteinExample(
        "4INS", "Porcine insulin", "Sus scrofa",
        "Compare insulin chains: this file contains four chains "
        "and 102 observed residues. Zinc is not shown.",
    ),
    ProteinExample(
        "1LYZ", "Lysozyme", "Gallus gallus",
        "A protective enzyme in egg white that acts on bacterial cell walls. "
        "Contains 129 observed residues.",
    ),
    ProteinExample(
        "1MBN", "Myoglobin", "Physeter macrocephalus",
        "Explore the composition of a globular protein with 153 observed residues. "
        "The heme group is not shown in this view.",
    ),
    ProteinExample(
        "1GFL", "Green fluorescent protein (GFP)", "Aequorea victoria",
        "A protein with a beta-barrel fold. The file contains two chains "
        "and 460 observed residues; the view does not represent the full chromophore chemistry.",
    ),
    ProteinExample(
        "2HHB", "Human hemoglobin", "Homo sapiens",
        "Compare four chains and 574 observed residues with myoglobin. "
        "Heme groups are not shown.",
    ),
)
