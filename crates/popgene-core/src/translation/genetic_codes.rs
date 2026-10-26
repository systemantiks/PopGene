//! NCBI gc.prt version 4.6, retrieved 2026-10-25.
//! Source: https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
//! TCAG order; Starts preserves separate initiation and termination flags.
pub(super) struct CodeData {
    pub id: u16,
    pub name: &'static str,
    pub amino_acids: &'static [u8; 64],
    pub starts: &'static [u8; 64],
}
pub(super) const CODES: &[CodeData] = &[
    CodeData {
        id: 1,
        name: "Standard",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M---------------M----------------------------",
    },
    CodeData {
        id: 2,
        name: "Vertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        starts: b"----------**--------------------MMMM----------**---M------------",
    },
    CodeData {
        id: 3,
        name: "Yeast Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**----------------------MM---------------M------------",
    },
    CodeData {
        id: 4,
        name: "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate\n Mitochondrial; Mycoplasma; Spiroplasma",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--MM------**-------M------------MMMM---------------M------------",
    },
    CodeData {
        id: 5,
        name: "Invertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        starts: b"---M------**--------------------MMMM---------------M------------",
    },
    CodeData {
        id: 6,
        name: "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
        amino_acids: b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    CodeData {
        id: 9,
        name: "Echinoderm Mitochondrial; Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    CodeData {
        id: 10,
        name: "Euplotid Nuclear",
        amino_acids: b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    CodeData {
        id: 11,
        name: "Bacterial, Archaeal and Plant Plastid",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M------------MMMM---------------M------------",
    },
    CodeData {
        id: 12,
        name: "Alternative Yeast Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    CodeData {
        id: 13,
        name: "Ascidian Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        starts: b"---M------**----------------------MM---------------M------------",
    },
    CodeData {
        id: 14,
        name: "Alternative Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"-----------*-----------------------M----------------------------",
    },
    CodeData {
        id: 15,
        name: "Blepharisma Macronuclear",
        amino_acids: b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    CodeData {
        id: 16,
        name: "Chlorophycean Mitochondrial",
        amino_acids: b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    CodeData {
        id: 21,
        name: "Trematode Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    CodeData {
        id: 22,
        name: "Scenedesmus obliquus Mitochondrial",
        amino_acids: b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"------*---*---*--------------------M----------------------------",
    },
    CodeData {
        id: 23,
        name: "Thraustochytrium Mitochondrial",
        amino_acids: b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--*-------**--*-----------------M--M---------------M------------",
    },
    CodeData {
        id: 24,
        name: "Rhabdopleuridae Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M------**-------M---------------M---------------M------------",
    },
    CodeData {
        id: 25,
        name: "Candidate Division SR1 and Gracilibacteria",
        amino_acids: b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**-----------------------M---------------M------------",
    },
    CodeData {
        id: 26,
        name: "Pachysolen tannophilus Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    CodeData {
        id: 27,
        name: "Karyorelict Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    CodeData {
        id: 28,
        name: "Condylostoma Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*--------------------M----------------------------",
    },
    CodeData {
        id: 29,
        name: "Mesodinium Nuclear",
        amino_acids: b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    CodeData {
        id: 30,
        name: "Peritrich Nuclear",
        amino_acids: b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    CodeData {
        id: 31,
        name: "Blastocrithidia Nuclear",
        amino_acids: b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    CodeData {
        id: 32,
        name: "Balanophoraceae Plastid",
        amino_acids: b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------*---*----M------------MMMM---------------M------------",
    },
    CodeData {
        id: 33,
        name: "Cephalodiscidae Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M-------*-------M---------------M---------------M------------",
    },
];
