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
        "1CRN", "Crambina", "Crambe hispanica subsp. abyssinica",
        "Una proteína vegetal pequeña: 46 residuos observados para explorar "
        "su forma compacta y composición.",
    ),
    ProteinExample(
        "1UBQ", "Ubiquitina humana", "Homo sapiens",
        "Participa en el etiquetado de proteínas para su degradación. "
        "Un ejemplo compacto de 76 residuos observados.",
    ),
    ProteinExample(
        "4INS", "Insulina porcina", "Sus scrofa",
        "Compara las cadenas de la insulina: este archivo contiene cuatro cadenas "
        "y 102 residuos observados. El zinc no se muestra.",
    ),
    ProteinExample(
        "1LYZ", "Lisozima", "Gallus gallus",
        "Una enzima de defensa de la clara de huevo que actúa sobre la pared "
        "bacteriana. Contiene 129 residuos observados.",
    ),
    ProteinExample(
        "1MBN", "Mioglobina", "Physeter macrocephalus",
        "Explora la composición de una proteína globular de 153 residuos observados. "
        "El grupo hemo no se muestra en esta vista.",
    ),
    ProteinExample(
        "1GFL", "Proteína fluorescente GFP", "Aequorea victoria",
        "Una proteína con plegamiento de barril beta. El archivo aporta dos cadenas "
        "y 460 residuos observados; la vista no representa toda la química del cromóforo.",
    ),
    ProteinExample(
        "2HHB", "Hemoglobina humana", "Homo sapiens",
        "Compara cuatro cadenas y 574 residuos observados con la mioglobina. "
        "Los grupos hemo no se muestran.",
    ),
)
