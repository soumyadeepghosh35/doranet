import dataclasses
import typing


@dataclasses.dataclass(frozen=True, slots=True)
class OperatorSmarts:
    name: str
    smarts: str
    reactants_stoi: typing.Optional[tuple[int, ...]] = None
    products_stoi: typing.Optional[tuple[int, ...]] = None
    enthalpy_correction: typing.Optional[float] = None
    ring_issue: typing.Optional[bool] = False
    kekulize_flag: typing.Optional[bool] = False
    Retro_Not_Aromatic: typing.Optional[bool] = False
    number_of_steps: typing.Optional[int] = 1


op_smarts = (
    # Alkenes ##############################################################
    # Hydrogenation of Alkene
    OperatorSmarts(
        "Hydrogenation of Alkene", "[C+0:1]=[C+0:2].[H][H]>>[*:1][*:2]",
        (1, 1),
        (1,)
    ),
    # Oxidative Cleavage of Alkenes
    OperatorSmarts(
        "Oxidative Cleavage of Alkenes",
        "[C+0:1]=!@[C+0:2].[O+0:4]=[O+0:5]>>[*:1]=[*:4].[*:2]=[*:5]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Oxidative Cleavage of Alkenes, Intramolecular",
        "[C+0:1]=@[C+0:2].[O+0:4]=[O+0:5]>>([*:1]=[*:4].[*:2]=[*:5])",
        (1, 1),
        (1,),
    ),
    # Olefin Metathesis
    # Cross Metathesis (CM)
    OperatorSmarts(
        "Olefin Cross Metathesis (CM)",
        "[C+0:1]=!@[C+0:2].[C+0:3]=!@[C+0:4]>>[*:1]=[*:3].[*:2]=[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Ring-Opening Metathesis (ROM)
    # Ring-Opening Metathesis 1:
    OperatorSmarts(
        "Olefin Ring-Opening Metathesis (ROM) 1",
        "[C+0:1]=@[C+0:2].[C+0:3]=!@[C+0:4]>>([*:1]=[*:3].[*:2]=[*:4])",
        (1, 1),
        (1,),
    ),
    # Ring-Opening Metathesis 2:
    OperatorSmarts(
        "Olefin Ring-Opening Metathesis (ROM) 2",
        "[C+0:1]=@[C+0:2].[C+0:3]=!@[C+0:4]>>([*:1]=[*:3].[*:2]=[*:3]).([*:1]=[*:4].[*:2]=[*:4])",
        (1, 1),
        (0.5, 0.5),
    ),
    # Ring-Closing Metathesis (RCM)
    OperatorSmarts(
        "Olefin Ring-Closing Metathesis (RCM)",
        "([C+0:1]=!@[C+0:3].[C+0:2]=!@[C+0:4])>>[*:1]=[*:2].[*:3]=[*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Reductive Cross-Coupling of Alkenes https://www.nature.com/articles/nature14006 https://www.nature.com/articles/s41467-022-30286-8
    OperatorSmarts(
        "Reductive Cross-Coupling of Alkenes",
        "[C+0:1]=[C+0:2].[C+0:3]=[C+0:4].[H][H]>>[*:1][*:2][*:3][*:4]",
        (1, 1, 1),
        (1,),
    ),
    # Epoxidation of Alkene
    OperatorSmarts(
        "Epoxidation of Alkene",
        "[C+0:1]=[C+0:2].[O+0:3]=[O+0]>>[*:1]1[*:2][*:3]1",
        (1, 0.5),
        (1,),
    ),
    # Hydration of Alkene, Ethers from Addition of Alcohols to Alkenes
    OperatorSmarts(
        "Hydration of Alkene, Ethers from Addition of Alcohols to Alkenes",
        "[C+0:1]=[C+0:2].[O+0!H0;!$(*-C=O):3]>>[*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Ethers from Addition of Alcohols to Alkenes, Intramolecular",
        "([C+0:1]=[C+0:2].[O+0H:3][C!$(*=O)+0:4])>>[*:1][*:2][*:3][*:4]",
        (1,),
        (1,),
    ),
    # Hydrohalogenation of Alkenes
    OperatorSmarts(
        "Hydrohalogenation of Alkenes",
        "[C+0:1]=[C+0:2].[F,Cl,Br,I;+0H:3]>>[*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    # Halogenation of Alkenes
    OperatorSmarts(
        "Halogenation of Alkenes",
        "[C+0:1]=[C+0:2].[F,Cl,Br,I;+0:3][F,Cl,Br,I;+0:4]>>[*:1]([*:3])[*:2][*:4]",
        (1, 1),
        (1,),
    ),
    # Diels-Alder Reaction
    OperatorSmarts(
        "Diels-Alder Reaction with Alkenes",
        "[C+0:1]=[C+0:2][C+0:3]=[C+0:4].[C+0:5]=[C+0:6]>>[*:1]1[*:2]=[*:3][*:4][*:6][*:5]1",
        (1, 1),
        (1,),
        kekulize_flag=True,
    ),
    OperatorSmarts(
        "Diels-Alder Reaction with Alkenes, Intramolecular",
        "([C+0:1]=[C+0:2][C+0:3]=[C+0:4].[C+0:5]=[C+0:6])>>[*:1]1[*:2]=[*:3][*:4][*:6][*:5]1",
        (1,),
        (1,),
        kekulize_flag=True,
        ring_issue=True,
    ),
    OperatorSmarts(
        "Diels-Alder Reaction with Alkynes",
        "[C+0:1]=[C+0:2][C+0:3]=[C+0:4].[C+0:5]#[C+0:6]>>[*:1]1[*:2]=[*:3][*:4][*:6]=[*:5]1",
        (1, 1),
        (1,),
        kekulize_flag=True,
    ),
    OperatorSmarts(
        "Diels-Alder Reaction with Alkynes, Intramolecular",
        "([C+0:1]=[C+0:2][C+0:3]=[C+0:4].[C+0:5]#[C+0:6])>>[*:1]1[*:2]=[*:3][*:4][*:6]=[*:5]1",
        (1,),
        (1,),
        kekulize_flag=True,
        ring_issue=True,
    ),
    # Halohydrin Formation
    OperatorSmarts(
        "Halohydrin Formation",
        "[C+0:1]=[C+0:2].[F,Cl,Br,I;+0:3][F,Cl,Br,I;+0:4].[O+0H2:5]>>[*:1]([*:5])[*:2][*:3].[*:4]",
        (1, 1, 1),
        (1, 1),
    ),
    # Diol Formation by Oxidation
    OperatorSmarts(
        "Diol Formation by Oxidation",
        "[C+0:1]=[C+0:2].[O+0:3]=[O+0].[O+0H2:4]>>[*:1]([*:3])[*:2][*:4]",
        (1, 0.5, 1),
        (1,),
    ),
    # Conjugated Dienes 1,4 Addition
    OperatorSmarts(
        "Conjugated Dienes 1,4 Addition with Hydrogen Halide",
        "[C+0:1]=[C+0:2][C+0:3]=[C+0:4].[F,Cl,Br,I;+0H:5]>>[*:5][*:1][*:2]=[*:3][*:4]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Conjugated Dienes 1,4 Addition with Halogens",
        "[C+0:1]=[C+0:2][C+0:3]=[C+0:4].[F,Cl,Br,I;+0:5][F,Cl,Br,I;+0:6]>>[*:5][*:1][*:2]=[*:3][*:4][*:6]",
        (1, 1),
        (1,),
    ),
    # Claisen Rearrangement, Cope Rearrangements
    OperatorSmarts(
        "Claisen Rearrangement, Cope Rearrangements",
        "[C+0:1]=[C+0:2][O,C;+0:3][C+0:4][C+0:5]=[C+0:6]>>[*:3]=[*:2][*:1][*:6][*:5]=[*:4]",
        (1,),
        (1,),
        ring_issue=True,
    ),
    # Aromatic Claisen Rearrangement
    OperatorSmarts(
        "Aromatic Claisen Rearrangement 1",
        "[c+0H:1]1:[c+0:2]:[c+0:3]:[c+0:4]:[c+0:5]:[c+0:6]:1-[O+0:7][C+0:8][C+0:9]=[C+0:10]>>[*:1]1([*:10][*:9]=[*:8]):[*:2]:[*:3]:[*:4]:[*:5]:[*:6]:1-[*:7]",
        (1,),
        (1,),
        ring_issue=True,
    ),
    OperatorSmarts(
        "Aromatic Claisen Rearrangement with Ortho Substituents",
        "[c+0H0:1]1:[c+0:2]:[c+0:3]:[c+0:4]:[c+0H0:5]:[c+0:6]:1-[O+0:7][C+0:8][C+0:9]=[C+0:10]>>[*:1]1:[*:2]:[*:3]([*:8][*:9]=[*:10]):[*:4]:[*:5]:[*:6]:1-[*:7]",
        (1,),
        (1,),
        ring_issue=True,
    ),
    # Tsuji–Trost Reaction
    # Nu = Active Methylenes
    OperatorSmarts(
        "Tsuji–Trost Reaction, Nu = Active Methylenes",
        "[C+0:1]=[C+0:2][C+0;!$(*[OH]);!$(*=O):3]-!@[F,Cl,Br,I,O$(*C=O);+0:4].[C+0:5](=[O+0:6])[C+0!H0:7][C+0:8]=[O+0:9]>>[*:1]=[*:2][*:3][*:7]([*:5]=[*:6])[*:8]=[*:9].[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Nu = Enolates
    OperatorSmarts(
        "Tsuji–Trost Reaction, Nu = Enolates",
        "[C+0:1]=[C+0:2][C+0;!$(*[OH]);!$(*=O):3]-!@[F,Cl,Br,I,O$(*C=O);+0:4].[C+0:5][C+0!H0:7][C+0;!$(*CC=O):8]=[O+0:9]>>[*:1]=[*:2][*:3][*:7]([*:5])[*:8]=[*:9].[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Nu = Amines
    OperatorSmarts(
        "Tsuji–Trost Reaction, Nu = Amines",
        "[C+0:1]=[C+0:2][C+0;!$(*[OH]);!$(*=O):3]-!@[F,Cl,Br,I,O$(*C=O);+0:4].[N+0!H0:7][C,c;+0:8]>>[*:1]=[*:2][*:3][*:7][*:8].[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Nu = Phenols
    OperatorSmarts(
        "Tsuji–Trost Reaction, Nu = Phenols",
        "[C+0:1]=[C+0:2][C+0;!$(*[OH]);!$(*=O):3]-!@[F,Cl,Br,I,O$(*C=O);+0:4].[O+0H:7][c+0:8]>>[*:1]=[*:2][*:3][*:7][*:8].[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Intramolecular Tsuji-Trost Reaction, Nu = Active Methylenes
    OperatorSmarts(
        "Intramolecular Tsuji-Trost Reaction, Nu = Active Methylenes",
        "([C+0:1]=[C+0:2][C+0;!$(*[OH]);!$(*=O):3]-!@[F,Cl,Br,I,O$(*C=O);+0:4].[C+0:5](=[O+0:6])[C+0!H0:7][C+0:8]=[O+0:9])>>[*:1]=[*:2][*:3][*:7]([*:5]=[*:6])[*:8]=[*:9].[*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Intramolecular Tsuji-Trost Reaction, Nu = Enolates
    OperatorSmarts(
        "Intramolecular Tsuji-Trost Reaction, Nu = Enolates",
        "([C+0:1]=[C+0:2][C+0;!$(*[OH]);!$(*=O):3]-!@[F,Cl,Br,I,O$(*C=O);+0:4].[C+0:5][C+0!H0:7][C+0;!$(*CC=O):8]=[O+0:9])>>[*:1]=[*:2][*:3][*:7]([*:5])[*:8]=[*:9].[*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Intramolecular Tsuji-Trost Reaction, Nu = Amines
    OperatorSmarts(
        "Intramolecular Tsuji-Trost Reaction, Nu = Amines",
        "([C+0:1]=[C+0:2][C+0;!$(*[OH]);!$(*=O):3]-!@[F,Cl,Br,I,O$(*C=O);+0:4].[N+0!H0:7][C,c;+0:8])>>[*:1]=[*:2][*:3][*:7][*:8].[*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Intramolecular Tsuji-Trost Reaction, Nu = Phenols
    OperatorSmarts(
        "Intramolecular Tsuji-Trost Reaction, Nu = Phenols",
        "([C+0:1]=[C+0:2][C+0;!$(*[OH]);!$(*=O):3]-!@[F,Cl,Br,I,O$(*C=O);+0:4].[O+0H:7][c+0:8])>>[*:1]=[*:2][*:3][*:7][*:8].[*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Alkynes ##############################################################
    # Halogen Addition of Alkynes
    OperatorSmarts(
        "Halogen Addition of Alkynes",
        "[C+0:1]#[C+0:2].[F,Cl,Br,I;+0:3][F,Cl,Br,I;+0:4]>>[*:1]([*:3])=[*:2][*:4]",
        (1, 1),
        (1,),
    ),
    # Hydrogen Halides Addition of Alkynes
    OperatorSmarts(
        "Hydrogen Halides Addition of Alkynes",
        "[C+0:1]#[C+0:2].[F,Cl,Br,I;+0H:3]>>[*:1]([*:3])=[*:2]",
        (1, 1),
        (1,),
    ),
    # Hydration of Alkynes
    OperatorSmarts(
        "Hydration of Alkynes",
        "[C+0:1]#[C+0:2].[O+0H2:3]>>[*:1]([*:3])=[*:2]",
        (1, 1),
        (1,),
    ),
    # Hydrogenation of Alkynes
    OperatorSmarts(
        "Hydrogenation of Alkynes", "[C+0:1]#[C+0:2].[H][H]>>[*:1]=[*:2]", (1, 1), (1,)
    ),
    OperatorSmarts(
        "Hydrogenation of Alkynes to Alkanes",
        "[C+0:1]#[C+0:2].[H][H]>>[*:1][*:2]",
        (1, 2),
        (1,),
    ),
    # Gentle Alkyne Oxidation
    OperatorSmarts(
        "Gentle Alkyne Oxidation",
        "[C+0:1]#[C+0:2].[O+0:4]=[O+0:5]>>[*:1](=[*:4])[*:2]=[*:5]",
        (1, 1),
        (1,),
    ),
    # Oxidative Cleavage of Alkynes
    OperatorSmarts(
        "Oxidative Cleavage of Alkynes",
        "[*:6][C+0:1]#!@[C+0:2][*:7].[O+0H2:3].[O+0:4]=[O+0:5]>>[*:6][*:1](=[*:3])[*:4].[*:7][*:2](=[*:5])[O]",
        (1, 1, 1.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Oxidative Cleavage of Alkynes Intramolecular",
        "[*:6][C+0:1]#@[C+0:2][*:7].[O+0H2:3].[O+0:4]=[O+0:5]>>([*:6][*:1](=[*:3])[*:4].[*:7][*:2](=[*:5])[O])",
        (1, 1, 1.5),
        (1,),
    ),
    OperatorSmarts(
        "Oxidative Cleavage of Terminal Alkynes",
        "[*:6][C+0:1]#[C+0H:2].[O+0:4]=[O+0:5]>>[*:6][*:1](=[*:4])[*:5].[O]=[*:2]=[O]",
        (1, 2),
        (1, 1),
    ),
    # Alkynes by Dehydrohalogenation
    OperatorSmarts(
        "Alkynes by Dehydrohalogenation of Alkanes",
        "[C+0!H0:1]([F,Cl,Br,I;+0:2])[C+0!H0:3]([F,Cl,Br,I;+0])>>[*:1]#[*:3].[*:2]",
        (1,),
        (1, 2),
    ),
    OperatorSmarts(
        "Alkynes by Dehydrohalogenation of Alkenes",
        "[C+0:1]([F,Cl,Br,I;+0:2])=[C+0!H0:3]>>[*:1]#[*:3].[*:2]",
        (1,),
        (1, 1),
    ),
    # Alkylation of Acetylide Anions with Methyl and Primary Haloalkanes
    # enthalpy correction from DFT: NaNH2 -> NH3 + NaBr - HBr  dH = -54.50 kcal/mol
    OperatorSmarts(
        "Alkylation of Acetylide Anions with Methyl and Primary Haloalkanes",
        "[C+0:1]#[C+0H:2].[F,Cl,Br,I;+0:3][C;H2,H3;+0:4]>>[*:1]#[*:2][*:4].[*:3]",
        (1, 1),
        (1, 1),
        enthalpy_correction=-54.50,
    ),
    # Sonogashira Coupling
    OperatorSmarts(
        "Sonogashira Coupling",
        "[C+0:1]#[C+0H:2].[F,Cl,Br,I;+0:3][c,C$(*=C),C$(*=O);+0:4]>>[*:1]#[*:2][*:4].[*:3]",
        (1, 1),
        (1, 1),
    ),
    # Addition of Alkynes to Aldehydes and Activated Ketones   doi.org/10.1021/ol026282k  doi.org/10.1021/jo701643h
    OperatorSmarts(
        "Addition of Alkynes to Aldehydes",
        "[C+0:1]#[C+0H:2].[CX3+0H:3]=[O+0:4]>>[*:1]#[*:2][*:3][*:4]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Addition of Alkynes to Activated Ketones",
        "[C+0:1]#[C+0H:2].[CX3$(*=O)+0:5][CX3!$(*O)+0:3]=[O+0:4]>>[*:1]#[*:2][*:3]([*:4])[*:5]",
        (1, 1),
        (1,),
    ),
    # Haloalkanes, Halogenation ###########################################
    # Alcohols React with Hydrogen Halides to form Haloalkanes, Carboxylic Acids Conversion to Acid Chlorides
    OperatorSmarts(
        "Alcohols React with Hydrogen Halides to form Haloalkanes, Carboxylic Acids Conversion to Acid Chlorides",
        "[C+0:1][O+0H:2].[F,Cl,Br,I;+0H:3]>>[*:1][*:3].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Haloalkane Hydrolysis, Acid Chlorides Hydrolysis, Haloarenes Hydrolysis
    OperatorSmarts(
        "Haloalkane Hydrolysis, Acid Chlorides Hydrolysis, Haloarenes Hydrolysis",
        "[C,c;+0:1][F,Cl,Br,I;+0:2].[O+0H2:3]>>[*:1][*:3].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Wurtz Reaction, Coupling of Halides with Gilman Reagent
    OperatorSmarts(
        "Wurtz Reaction, Coupling of Halides with Gilman Reagent",
        "[C,c;+0:1][F,Cl,Br,I;+0:2].[C,c;+0:3][F,Cl,Br,I;+0:4]>>[*:1][*:3].[*:2][*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Wurtz Reaction, Coupling of Halides with Gilman Reagent, Intramolecular",
        "([C,c;+0:1][F,Cl,Br,I;+0:2].[C,c;+0:3][F,Cl,Br,I;+0:4])>>[*:1][*:3].[*:2][*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Dehydrohalogenation
    OperatorSmarts(
        "Dehydrohalogenation",
        "[C+0:1]([F,Cl,Br,I;+0:2])[C+0!H0:3]>>[*:1]=[*:3].[*:2]",
        (1,),
        (1, 1),
    ),
    # Aromatic Halogenation
    OperatorSmarts(
        "Aromatic Halogenation",
        "[c+0H:1].[F,Cl,Br,I;+0:2][F,Cl,Br,I;+0:3]>>[*:1][*:2].[*:3]",
        (1, 1),
        (1, 1),
    ),
    # Friedel–Crafts Reaction
    OperatorSmarts(
        "Friedel–Crafts Reaction",
        "[c+0H:1].[CX4,C$(*=O);+0:2][F,Cl,Br,I;+0:3]>>[*:1][*:2].[*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Friedel–Crafts Reaction, Intramolecular",
        "([c+0H:1].[CX4,C$(*=O);+0:2][F,Cl,Br,I;+0:3])>>[*:1][*:2].[*:3]",
        (1,),
        (1, 1),
    ),
    OperatorSmarts(
        "Friedel–Crafts Acylation with Carboxylic Acids",  # doi.org/10.1021/ol800752v
        "[c+0H:1].[C$(*=O)+0:2][O+0H:3]>>[*:1][*:2].[*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Friedel–Crafts Acylation with Carboxylic Acids, Intramolecular",
        "([c+0H:1].[C$(*=O)+0:2][O+0H:3])>>[*:1][*:2].[*:3]",
        (1,),
        (1, 1),
    ),
    OperatorSmarts(
        "Friedel–Crafts Acylation with Acid Anhydrides",
        "[c+0H:1].[C$(*=O)+0:2][O+0:3]!@[C$(*=O)+0:4]>>[*:1][*:2].[*:3][*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Friedel–Crafts Reaction with Alkenes",
        "[c+0H:1].[C+0:2]=[C+0:3]>>[*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Friedel–Crafts Reaction with Alkenes, Intramolecular",
        "([c+0H:1].[C+0:2]=[C+0:3])>>[*:1][*:2][*:3]",
        (1,),
        (1,),
    ),
    # Heck Reaction     possible with alkyl halides: doi.org/10.1039/C4CC00297K
    OperatorSmarts(
        "Heck Reaction",
        "[C+0:1]=[C+0!H0:2].[C,c;+0:3][F,Cl,Br,I;+0:4]>>[*:1]=[*:2][*:3].[*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Heck Reaction, Intramolecular",
        "([C+0:1]=[C+0!H0:2].[C,c;+0:3][F,Cl,Br,I;+0:4])>>[*:1]=[*:2][*:3].[*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Suzuki Coupling (Combined with Grignard Reagent and Borate Ester Formation, or Hydroboration of Alkenes or Alkynes)
    # NOT BALANCED, needs to add a correction term for enthalpy calculation
    # Boranes from R-X
    # enthalpy correction from DFT: Mg + B(OCH3)3 -> BrMgOCH3 + B(OCH3)2Br  dH = -42.13 kcal/mol
    OperatorSmarts(
        "Suzuki Coupling with Boranes from Alkyl, Alkenyl, or Aryl Halides",
        "[c,C!$(*~O);!$(*#*);+0:1][F,Cl,Br,I;+0].[c,C$(*=,#C);+0:2][F,Cl,Br,I;+0]>>[*:1][*:2]",
        (1, 1),
        (1,),
        enthalpy_correction=-42.13,
    ),
    # Boranes from alkenes or alkynes
    # enthalpy correction from DFT: (Sia)2BH -> (Sia)2BBr   dH = -26.11 kcal/mol
    OperatorSmarts(
        "Suzuki Coupling with Boranes from Alkenes",
        "[C+0:1]=[C+0:2].[c,C$(*=,#C);+0:3][F,Cl,Br,I;+0]>>[*:1]-[*:2][*:3]",
        (1, 1),
        (1,),
        enthalpy_correction=-26.11,
    ),
    OperatorSmarts(
        "Suzuki Coupling with Boranes from Alkynes",
        "[C+0:1]#[C+0:2].[c,C$(*=,#C);+0:3][F,Cl,Br,I;+0]>>[*:1]=[*:2][*:3]",
        (1, 1),
        (1,),
        enthalpy_correction=-26.11,
    ),
    # Reduction of Haloalkanes           # aryl: # doi.org/10.1021/ja01191a035
    OperatorSmarts(
        "Reduction of Haloalkanes",
        "[CX4,c;+0:1][F,Cl,Br,I;+0:2].[H][H]>>[*:1].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Grignard Reagent with Oxiranes
    # NOT BALANCED, needs to add a correction term for enthalpy calculation
    # enthalpy correction from DFT: Mg + HBr -> MgBr2  dH = -61.15 kcal/mol
    OperatorSmarts(
        "Grignard Reagent with Oxiranes",
        "[c,C!$(*~O);!$(*#*);+0:1][F,Cl,Br,I;+0].[C+0:2]1[C+0:3][O+0:4]1>>[*:1][*:2][*:3][*:4]",
        (1, 1),
        (1,),
        enthalpy_correction=-61.15,
    ),
    # Dichlorocarbene with Alkenes
    OperatorSmarts(
        "Dichlorocarbene with Alkenes",
        "[C+0H:1]([F,Cl,Br,I;+0:2])([F,Cl,Br,I;+0:3])[F,Cl,Br,I;+0:4].[C+0:5]=[C+0:6]>>[*:1]1([*:2])([*:3])[*:5][*:6]1.[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Simmons-Smith Reaction
    OperatorSmarts(
        "Simmons-Smith Reaction",
        "[C+0H2:1]([F,Cl,Br,I;+0:2])[F,Cl,Br,I;+0:3].[C+0:5]=[C+0:6]>>[*:1]1[*:5][*:6]1.[*:2].[*:3]",
        (1, 1),
        (1, 1, 1),
    ),
    # Esterification with Alkyl and Aryl Halides        aryl:doi.org/10.1021/acs.joc.1c00204
    OperatorSmarts(
        "Esterification with Alkyl and Aryl Halides",
        "[C+0:2](=[O+0:3])[O+0H:4].[CX4,c;+0:5][F,Cl,Br,I;+0:6]>>[*:2](=[*:3])[*:4][*:5].[*:6]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Esterification with Alkyl and Aryl Halides, Intramolecular",
        "([C+0:2](=[O+0:3])[O+0H:4].[CX4,c;+0:5][F,Cl,Br,I;+0:6])>>[*:2](=[*:3])[*:4][*:5].[*:6]",
        (1,),
        (1, 1),
    ),
    # Alcohols #############################################################
    # Glycol Cleavage by Oxidation
    OperatorSmarts(
        "Glycol Cleavage by Oxidation",
        "[C+0:1]([O+0H:5])-!@[C+0:2][O+0H:6].[O+0:3]=[O+0]>>[*:1]=[*:5].[*:2]=[*:6].[*:3]",
        (1, 0.5),
        (1, 1, 1),
    ),
    OperatorSmarts(
        "Glycol Cleavage by Oxidation Intramolecular",
        "[C+0:1]([O+0H:5])-@[C+0:2][O+0H:6].[O+0:3]=[O+0]>>([*:1]=[*:5].[*:2]=[*:6]).[*:3]",
        (1, 0.5),
        (1, 1),
    ),
    # Hydrodeoxygenation of Alcohol, Classic Synthesis of Aldehydes from Carboxylic Acids
    OperatorSmarts(
        "Hydrodeoxygenation of Alcohol, Classic Synthesis of Aldehydes from Carboxylic Acids",
        "[C,c;+0:1][O+0H:2].[H][H]>>[*:1].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Dehydration of Alcohols
    OperatorSmarts(
        "Dehydration of Alcohol",
        "[C+0!H0:1][C+0:2][O+0H:3]>>[*:1]=[*:2].[*:3]",
        (1,),
        (1, 1),
    ),
    # Dehydration of Alcohols, 2-step
    OperatorSmarts(
        "Dehydration of Alcohol, 2-step",
        "([C+0!H0:1][C+0:2][O+0H:3].[C+0!H0:4][C+0:5][O+0H])>>([*:1]=[*:2].[*:4]=[*:5]).[*:3]",
        (1,),
        (1, 2),
        number_of_steps=2,
    ),
    # Selective Oxidation of Alcohols
    OperatorSmarts(
        "Selective Oxidation of Alcohols",
        "[C+0!H0:1][O+0H:2].[O+0:3]=[O+0]>>[*:1]=[*:2].[*:3]",
        (1, 0.5),
        (1, 1),
    ),
    # Diol Carboxylation
    OperatorSmarts(
        "Diol Carboxylation",
        "[O+0H:1][C+0:2][C+0:3][O+0H:4].[O+0:5]=[C+0:6]=[O+0:7]>>[*:1].[*:2]1[*:3][*:4][*:6](=[*:7])[*:5]1",
        (1, 1),
        (1, 1),
    ),
    # add n=1, doi.org/10.1021/cs500301d
    OperatorSmarts(
        "Diol Carboxylation 2",
        "[O+0H:1][C+0:2][C+0:8][C+0:3][O+0H:4].[O+0:5]=[C+0:6]=[O+0:7]>>[*:1].[*:2]1[*:8][*:3][*:4][*:6](=[*:7])[*:5]1",
        (1, 1),
        (1, 1),
    ),
    # Hydrogenolysis of Primary Alcohol    doi:10.1021/ja01330a506
    OperatorSmarts(
        "Hydrogenolysis of Primary Alcohol",
        "[C+0:1][C+0H2:3][O+0H:2].[H][H]>>[*:1].[*:3].[*:2]",
        (1, 2),
        (1, 1, 1),
    ),
    # Pinacol Rearrangement
    # no H on C2
    OperatorSmarts(
        "Pinacol Rearrangement, no H on Carbon#2",
        "[C+0:1][C+0:2]([C+0:3])([O+0H:4])[C+0X4:5][O+0H:8]>>[*:1][*:2](=[*:4])[*:5]([*:3]).[*:8]",
        (1,),
        (1, 1),
    ),
    # 1 H on C2
    OperatorSmarts(
        "Pinacol Rearrangement, 1 H on Carbon#2",
        "[C+0H:2]([C+0:3])([O+0H:4])[C+0X4:5][O+0H:8]>>[*:2](=[*:4])([*:3])[*:5].[*:8]",
        (1,),
        (1, 1),
    ),
    # 2 Hs on C2
    OperatorSmarts(
        "Pinacol Rearrangement, 2 Hs on Carbon#2",
        "[C+0H2:2]([O+0H:4])[C+0X4:5][O+0H:8]>>[*:2](=[*:4])[*:5].[*:8]",
        (1,),
        (1, 1),
    ),
    # Dehydration of Geminal Diol
    OperatorSmarts(
        "Dehydration of Geminal Diol",
        "[C+0:1]([O+0H:2])[O+0H:3]>>[*:1]=[*:2].[*:3]",
        (1,),
        (1, 1),
    ),
    # Formation of Acetals from Hemiacetals
    OperatorSmarts(
        "Formation of Acetals from Hemiacetals",
        "[C+0X4:1]([O+0H0X2:2])[O+0H:3].[O+0H:4][C+0X4;!$(*OC):5]>>[*:1]([*:2])[*:4][*:5].[*:3]",
        (1, 1),
        (1, 1),
    ),
    # Oxidation Of Primary Alcohols to Carboxylic Acids
    OperatorSmarts(
        "Oxidation Of Primary Alcohols to Carboxylic Acids",
        "[C+0;H3,H2:1][O+0H:2].[O+0:3]=[O+0:4]>>[*:1](=[*:3])[*:2].[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Oxidation Of Primary Alcohols to Carboxylic Acids, 2-step
    OperatorSmarts(
        "Oxidation Of Primary Alcohols to Carboxylic Acids, 2-step",
        "([C+0;H3,H2:1][O+0H:2].[C+0;H3,H2:5][O+0H:6]).[O+0:3]=[O+0:4]>>([*:1](=[*:3])[*:2].[*:5](=[*:6])[O]).[*:4]",
        (1, 2),
        (1, 2),
        number_of_steps=2,
    ),
    # Formation of Cyclic Acetals from Ketones/Aldehydes with Diols
    OperatorSmarts(
        "Formation of Cyclic Acetals from Ketones/Aldehydes with Diols",
        "[C+0X3!$(*O):1]=[O+0:2].[O+0H:3][C+0X4:4][C+0X4:5][O+0H:6]>>[*:1]1[*:3][*:4][*:5][*:6]1.[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Ethers, Epoxides ####################################################
    # Hydrolysis of Ethers, Esters, Anhydrides
    OperatorSmarts(
        "Hydrolysis of Ethers, Esters, Anhydrides",
        "[C,c;+0:1][O+0:2]-!@[C,c;+0:3].[O+0H2:4]>>[*:1][*:2].[*:3][*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Hydrolysis of Ethers, Esters, Anhydrides, Intramolecular",
        "[C,c;+0:1][O+0:2]-@[C,c;+0:3].[O+0H2:4]>>([*:1][*:2].[*:3][*:4])",
        (1, 1),
        (1,),
    ),
    # Hydrogenolysis of Ethers
    OperatorSmarts(
        "Hydrogenolysis of Ethers",
        "[C,c;+0;!$(*=O):1]!@[O+0:2][C,c;+0;!$(*=O):3].[H][H]>>[*:1].[*:2][*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Hydrogenolysis of Ethers Intramolecular",
        "[C,c;+0;!$(*=O):1]@[O+0:2][C,c;+0;!$(*=O):3].[H][H]>>([*:1].[*:2][*:3])",
        (1, 1),
        (1,),
    ),
    # Williamson Ether Synthesis, Ullmann Condensation
    OperatorSmarts(
        "Williamson Ether Synthesis, Ullmann Condensation",
        "[C,c;!$(*=O);+0:1][O+0H:2].[CX4!H0,c;+0:3][F,Cl,Br,I;+0:4]>>[*:1][*:2][*:3].[*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Williamson Ether Synthesis, Ullmann Condensation, Intramolecular",
        "([C,c;!$(*=O);+0:1][O+0H:2].[CX4!H0,c;+0:3][F,Cl,Br,I;+0:4])>>[*:1][*:2][*:3].[*:4]",
        (1,),
        (1, 1),
    ),
    # Ether Cleavage
    OperatorSmarts(
        "Ether Cleavage",
        "[C,c;+0;!$(*=O):1][O+0:2]!@[C;+0;!$(*=O):3].[F,Cl,Br,I;+0H:4]>>[*:1][*:2].[*:3][*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Ether Cleavage, Intramolecular",
        "[C,c;+0;!$(*=O):1][O+0:2]@[C;+0;!$(*=O):3].[F,Cl,Br,I;+0H:4]>>([*:1][*:2].[*:3][*:4])",
        (1, 1),
        (1,),
    ),
    # Ether Synthesis by Dehydration
    OperatorSmarts(
        "Ether Synthesis by Dehydration",
        "[C,c;!$(*=O);+0:1][O+0H:2].[C,c;!$(*=O);+0:3][O+0H:4]>>[*:1][*:2][*:3].[*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Ether Synthesis by Dehydration, Intramolecular",
        "([C,c;!$(*=O);+0:1][O+0H:2].[C,c;!$(*=O);+0:3][O+0H:4])>>[*:1][*:2][*:3].[*:4]",
        (1,),
        (1, 1),
    ),
    # Epoxides Ring Opening
    OperatorSmarts(
        "Epoxides Ring Opening",
        "[C+0:1]1[C+0:2][O+0:3]1.[O+0!H0;!$(*-C=O):4]>>[*:4][*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    # Aldehydes and Ketones ###############################################
    # Aldehyde Oxidation
    OperatorSmarts(
        "Aldehyde Oxidation",
        "[C+0!H0:1]=[O+0:2].[O+0:3]=[O+0]>>[*:1](=[*:2])[*:3]",
        (1, 0.5),
        (1,),
    ),
    # Aldehyde Oxidation, 2-step
    OperatorSmarts(
        "Aldehyde Oxidation, 2-step",
        "([C+0!H0:1]=[O+0:2].[C+0!H0:5]=[O+0:6]).[O+0:3]=[O+0:4]>>([*:1](=[*:2])[*:3].[*:5](=[*:6])[*:4])",
        (1, 1),
        (1,),
        number_of_steps=2,
    ),
    # Aldehyde & Alcohol Oxidation
    OperatorSmarts(
        "Aldehyde & Alcohol Oxidation, 2-step",
        "([C+0!H0:1]=[O+0:2].[C+0;H3,H2:5][O+0H:6]).[O+0:3]=[O+0:4]>>([*:1](=[*:2])[*:3].[*:5](=[*:6])[O]).[*:4]",
        (1, 1.5),
        (1, 1),
        number_of_steps=2,
    ),
    # Ketone Oxidation
    OperatorSmarts(
        "Ketone Oxidation",
        "[C,c;+0:1][C+0:2](=[O+0:3])!@[C+0H2:4][C,c;+0:5].[O+0:6]=[O+0:7]>>[*:1][*:2](=[*:3])[*:6].[*:5][*:4](=[*:7])[O]",
        (1, 1.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Ketone Oxidation, Intramolecular",
        "[C,c;+0:1][C+0:2](=[O+0:3])@[C+0H2:4][C,c;+0:5].[O+0:6]=[O+0:7]>>([*:1][*:2](=[*:3])[*:6].[*:5][*:4](=[*:7])[O])",
        (1, 1.5),
        (1,),
    ),
    # Baeyer-Villiger Oxidation (ketones and aldehydes)
    OperatorSmarts(
        "Baeyer-Villiger Oxidation (Ketones)",
        "[C,c;+0:5][C+0:1](=[O+0:2])[C,c;+0:3].[O+0:4]=[O+0]>>[*:5][*:1](=[*:2])[*:4][*:3]",
        (1, 0.5),
        (1,),
    ),
    OperatorSmarts(
        "Baeyer-Villiger Oxidation (Aldehydes)",  # with aldehydes doi.org/10.1002/0471264180.or043.03
        "[C+0H:1](=[O+0:2])[C$(*(C)(C)C),C$(*=*),C$(*c),c;+0:3].[O+0:4]=[O+0]>>[*:1](=[*:2])[*:4][*:3]",
        (1, 0.5),
        (1,),
    ),
    # Hydrogenation of Ketones
    OperatorSmarts(
        "Hydrogenation of Ketones",
        "[C+0!$(*-O):1]=[O+0:2].[H][H]>>[*:1][*:2]",
        (1, 1),
        (1,),
    ),
    # Keto-enol Tautomerization        Also covers Imine <-> Enamine, Imine <-> Amine, Amide <-> Imidic Acid Tautomerization
    OperatorSmarts(
        "Keto-enol Tautomerization",
        "[C,N;+0!H0:1][C+0:2]=[O,N;+0:3]>>[*:1]=[*:2][*:3]",
        (1,),
        (1,),
        kekulize_flag=True,
    ),
    OperatorSmarts(
        "Keto-enol Tautomerization Reverse",
        "[C,N;+0:1]=[C+0:2][O,N;+0!H0:3]>>[*:1][*:2]=[*:3]",
        (1,),
        (1,),
        kekulize_flag=True,
    ),
    # Decarbonylation of Aldehydes
    OperatorSmarts(
        "Decarbonylation of Aldehydes",
        "[*+0:1][C+0H:2]=[O+0:3]>>[*:1].[*-:2]#[*+:3]",
        (1,),
        (1, 1),
    ),
    # Ketonization
    OperatorSmarts(
        "Ketonization",
        "[*+0:1][C+0:2](=[O+0:3])[O+0H:4].[*+0:5][C+0:6](=[O+0:7])[O+0H:8]>>[*:1][*:2](=[*:3])[*:5].[*:4]=[*:6]=[*:7].[*:8]",
        (1, 1),
        (1, 1, 1),
    ),
    OperatorSmarts(
        "Ketonization, Intramolecular",
        "([*+0:1][C+0:2](=[O+0:3])[O+0H:4].[*+0:5][C+0:6](=[O+0:7])[O+0H:8])>>[*:1][*:2](=[*:3])[*:5].[*:4]=[*:6]=[*:7].[*:8]",
        (1,),
        (1, 1, 1),
    ),
    # Oxidative Esterification of Aldehydes and Alcohols
    OperatorSmarts(
        "Oxidative Esterification of Aldehydes and Alcohols",
        "[C+0;!$(*=O):5][O+0H:1].[O+0:2]=[C+0!H0:3].[O+0:4]=[O+0]>>[*:5][*:1][*:3]=[*:2].[*:4]",
        (1, 1, 0.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Oxidative Esterification of Aldehydes and Alcohols, Intramolecular",
        "([C+0;!$(*=O):5][O+0H:1].[O+0:2]=[C+0!H0:3]).[O+0:4]=[O+0]>>[*:5][*:1][*:3]=[*:2].[*:4]",
        (1, 0.5),
        (1, 1),
    ),
    # McMurry Reaction
    OperatorSmarts(
        "McMurry Reaction",
        "[*+0:1][C+0;!$(*-O):2](=[O+0:3]).[*+0:5][C+0;!$(*-O):6](=[O+0:7]).[H][H]>>[*:1][*:2]=[*:6][*:5].[*:3].[*:7]",
        (1, 1, 2),
        (1, 1, 1),
    ),
    OperatorSmarts(
        "McMurry Reaction, Intramolecular",
        "([*+0:1][C+0;!$(*-O):2](=[O+0:3]).[*+0:5][C+0;!$(*-O):6](=[O+0:7])).[H][H]>>[*:1][*:2]=[*:6][*:5].[*:3].[*:7]",
        (1, 2),
        (1, 1, 1),
        ring_issue=True,
    ),
    # Hydration of Ketone and Aldehyde
    OperatorSmarts(
        "Hydration of Ketone and Aldehyde",
        "[C+0;!$(*-O):2]=[O+0:3].[O+0H2:4]>>[*:2]([*:3])[*:4]",
        (1, 1),
        (1,),
    ),
    # Hemiacetal Dissociation and Formation, Addition of Alcohols to Carbonyl Groups
    OperatorSmarts(
        "Hemiacetal Dissociation, Addition of Alcohols to Carbonyl Groups Reverse",
        "[C;!$(*=O):1]([O+0H:2])!@[O+0:3][C;!$(*=O):4]>>[*:1]=[*:2].[*:3][*:4]",
        (1,),
        (1, 1),
    ),
    OperatorSmarts(
        "Hemiacetal Dissociation, Addition of Alcohols to Carbonyl Groups Reverse, Intramolecular",
        "[C;!$(*=O):1]([O+0H:2])@[O+0:3][C;!$(*=O):4]>>([*:1]=[*:2].[*:3][*:4])",
        (1,),
        (1,),
    ),
    OperatorSmarts(
        "Hemiacetal Formation, Addition of Alcohols to Carbonyl Groups",
        "[CX3;!$(*O):1]=[O+0:2].[O+0H:3][C;!$(*=O):4]>>[*:1]([*:2])[*:3][*:4]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Hemiacetal Formation, Addition of Alcohols to Carbonyl Groups, Intramolecular",
        "([C;!$(*O):1]=[O+0:2].[O+0H:3][C;!$(*=O):4])>>[*:1]([*:2])[*:3][*:4]",
        (1,),
        (1,),
    ),
    # Reduction of Carbonyl Groups
    OperatorSmarts(
        "Reduction of Carbonyl Groups",
        "[C,c,N;+0:1][C+0:2]=[O+0:3].[H][H]>>[*:1][*:2].[*:3]",
        (1, 2),
        (1, 1),
    ),
    # a-Halogenation
    OperatorSmarts(
        "a-Halogenation",
        "[C+0:2](=[O+0:3])[C+0!H0:4].[F,Cl,Br,I;+0:5][F,Cl,Br,I;+0:6]>>[*:2](=[*:3])[*:4][*:5].[*:6]",
        (1, 1),
        (1, 1),
    ),
    # Wittig Reaction (Combined with Ylide Formation)
    # !NOT BALANCED! Reagents and by-products ignored for simplicity sake, needs to add a correction term for enthalpy calculation
    # Enthalpy correction from DFT: Ph3P + BuLi -> Ph3P=O + Bu + LiBr, delta H = -141.69 kcal/mol, gas phase 298.15 K
    OperatorSmarts(
        "Wittig Reaction (Combined with Ylide Formation)",
        "[C+0!H0;!$(*~[O,S,N]):1][F,Cl,Br,I;+0].[C,c;+0:3][C+0;!$(*-O):4]=[O+0]>>[*:1]=[*:4][*:3]",
        (1, 1),
        (1,),
        enthalpy_correction=-141.69,
    ),
    OperatorSmarts(
        "Wittig Reaction, Intramolecular (Combined with Ylide Formation)",
        "([C+0!H0;!$(*~[O,S,N]):1][F,Cl,Br,I;+0].[C,c;+0:3][C+0;!$(*-O):4]=[O+0])>>[*:1]=[*:4][*:3]",
        (1,),
        (1,),
        enthalpy_correction=-141.69,
        ring_issue=True,
    ),
    # Carboxylic Acids ####################################################
    # Carboxylic Acids Decarboxylation
    OperatorSmarts(
        "Carboxylic Acids Decarboxylation",
        "[*+0:1][C+0:2](=[O+0:3])[O+0H:4]>>[*:1].[*:3]=[*:2]=[*:4]",
        (1,),
        (1, 1),
    ),
    # Acid Anhydrides React with Alcohols to Form Esters
    OperatorSmarts(
        "Acid Anhydrides React with Alcohols to Form Esters",
        "[C+0;!$(*[OH]):1](=[O+0:2])[O+0!R:3][C+0;!$(*[OH]):4](=[O+0:5]).[C,c;+0;!$(*=O):6][O+0H:7]>>[*:1](=[*:2])[*:3][*:6].[*:4](=[*:5])[*:7]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Acid Anhydrides React with Alcohols to Form Esters, Intramolecular",
        "[C+0;!$(*[OH]):1](=[O+0:2])[O+0R:3][C+0;!$(*[OH]):4](=[O+0:5]).[C,c;+0;!$(*=O):6][O+0H:7]>>([*:1](=[*:2])[*:3][*:6].[*:4](=[*:5])[*:7])",
        (1, 1),
        (1,),
    ),
    # Enols and Enolates ##################################################
    # Aldol Condensation
    OperatorSmarts(
        "Aldol Condensation",
        "[C+0:2](=[O+0:3])[C+0!H0:4].[O+0:5]=[C+0X3:6]>>[*:2](=[*:3])[*:4][*:6]([*:5])",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Aldol Condensation, Intramolecular",
        "([C+0:2](=[O+0:3])[C+0!H0:4].[O+0:5]=[C+0X3:6])>>[*:2](=[*:3])[*:4][*:6]([*:5])",
        (1,),
        (1,),
        ring_issue=True,
    ),
    OperatorSmarts(
        "Aldol Condensation (2H)",
        "[C+0:2](=[O+0:3])[C+0;H3,H2:4].[O+0:5]=[C+0X3:6]>>[*:2](=[*:3])[*:4]=[*:6].[*:5]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Aldol Condensation (2H), Intramolecular",
        "([C+0:2](=[O+0:3])[C+0;H3,H2:4].[O+0:5]=[C+0X3:6])>>[*:2](=[*:3])[*:4]=[*:6].[*:5]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Claisen Condensation (including crossed)
    OperatorSmarts(
        "Claisen Condensation",
        "[C+0:2](=[O+0:3])!@[O+0:7][C+0:4].[O+0:5]=[C+0:6][C+0!H0:8]>>[*:2](=[*:3])[*:8][*:6]=[*:5].[*:7][*:4]",
        (1, 1),
        (1, 1),
    ),
    # Dieckmann Condensation (intramolecular version of Claisen Condensation)
    OperatorSmarts(
        "Dieckmann Condensation",
        "([C+0:2](=[O+0:3])!@[O+0:7][C+0:4].[O+0:5]=[C+0:6][C+0!H0:8])>>[*:2](=[*:3])[*:8][*:6]=[*:5].[*:7][*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Michael Reaction
    OperatorSmarts(
        "Michael Reaction",
        "[C+0:1](=[O+0:2])[C+0!H0:3][C+0:4]=[O+0:5].[C+0:6]=[C+0:7][C+0:8]=[O+0:9]>>[*:1](=[*:2])[*:3]([*:6][*:7][*:8]=[*:9])[*:4]=[*:5]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Michael Reaction, Intramolecular",
        "([C+0:1](=[O+0:2])[C+0!H0:3][C+0:4]=[O+0:5].[C+0:6]=[C+0:7][C+0:8]=[O+0:9])>>[*:1](=[*:2])[*:3]([*:6][*:7][*:8]=[*:9])[*:4]=[*:5]",
        (1,),
        (1,),
        ring_issue=True,
    ),
    # with a 6-membered ring
    OperatorSmarts(
        "Michael Reaction with Cyclic Ketones",
        "[C+0:1]1(=[O+0:2])[C+0!H0:3][C+0:10][C+0:11][C+0:12][C+0:13]1.[C+0:6]=[C+0:7][C+0:8]=[O+0:9]>>[*:1]1(=[*:2])[*:3]([*:10][*:11][*:12][*:13]1)[*:6][*:7][*:8]=[*:9]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Michael Reaction with Cyclic Ketones, Intramolecular",
        "([C+0:1]1(=[O+0:2])[C+0!H0:3][C+0:10][C+0:11][C+0:12][C+0:13]1.[C+0:6]=[C+0:7][C+0:8]=[O+0:9])>>[*:1]1(=[*:2])[*:3]([*:10][*:11][*:12][*:13]1)[*:6][*:7][*:8]=[*:9]",
        (1,),
        (1,),
    ),
    # Robinson Annulation (Michael-aldol sequence)
    OperatorSmarts(
        "Robinson Annulation",
        "[C+0:1]1(=[O+0:2])[C+0!H0:3][C+0:14][C+0:15][C+0:16][C+0:17]1.[C+0:8]=[C+0:9][C+0:10](=[O+0:11])[C+0;H2,H3:12]>>[*:1]12[*:3]([*:8][*:9][*:10](=[*:11])[*:12]=2)[*:14][*:15][*:16][*:17]1.[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Esters ##############################################################
    # Hydrolysis of Esters (Covered by hydrolysis of ethers)
    # Esterification, Acid Anhydride Formation
    OperatorSmarts(
        "Esterification, Acid Anhydride Formation",
        "[C+0:2](=[O+0:3])[O+0H:4].[C,c;+0:5][O+0H:6]>>[*:2](=[*:3])[*:6][*:5].[*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Esterification, Acid Anhydride Formation, Intramolecular",
        "([C+0:2](=[O+0:3])[O+0H:4].[C,c;+0:5][O+0H:6])>>[*:2](=[*:3])[*:6][*:5].[*:4]",
        (1,),
        (1, 1),
    ),
    # Esterification of Acid Halides
    OperatorSmarts(
        "Esterification of Acid Halides",
        "[C+0:1](=[O+0:2])[F,Cl,Br,I;+0:3].[C,c;!$(*=O);+0:4][O+0H:5]>>[*:1](=[*:2])[*:5][*:4].[*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Esterification of Acid Halides, Intramolecular",
        "([C+0:1](=[O+0:2])[F,Cl,Br,I;+0:3].[C,c;!$(*=O);+0:4][O+0H:5])>>[*:1](=[*:2])[*:5][*:4].[*:3]",
        (1,),
        (1, 1),
    ),
    # Ester Reduction
    # To Aldehydes
    OperatorSmarts(
        "Ester Reduction to Aldehydes and Alcohols",
        "[O+0H0!R:1][C+0:2]=[O+0:3].[H][H]>>[*:1].[*:2]=[*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Ester Reduction to Aldehydes and Alcohols, Intramolecular",
        "[O+0H0R:1][C+0:2]=[O+0:3].[H][H]>>([*:1].[*:2]=[*:3])",
        (1, 1),
        (1,),
    ),
    # To Alcohols
    OperatorSmarts(
        "Ester Reduction to Alcohols",
        "[O+0H0!R:1][C+0:2]=[O+0:3].[H][H]>>[*:1].[*:2][*:3]",
        (1, 2),
        (1, 1),
    ),
    OperatorSmarts(
        "Ester Reduction to Alcohols, Intramolecular",
        "[O+0H0R:1][C+0:2]=[O+0:3].[H][H]>>([*:1].[*:2][*:3])",
        (1, 2),
        (1,),
    ),
    # Transesterification
    OperatorSmarts(
        "Transesterification",
        "[C+0:1](=[O+0:2])!@[O+0:3][*+0:4].[*+0:6][O+0H;!$(*-C=O):5]>>[*:1](=[*:2])[*:5][*:6].[*:4][*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Transesterification, Intramolecular",
        "([C+0:1](=[O+0:2])!@[O+0:3][*+0:4].[*+0:6][O+0H;!$(*-C=O):5])>>[*:1](=[*:2])[*:5][*:6].[*:4][*:3]",
        (1,),
        (1, 1),
        ring_issue=True,
        kekulize_flag=True,
    ),
    # Benzene and Derivatives #############################################
    # Benzene Hydrogenation
    OperatorSmarts(
        "Benzene Hydrogenation",
        "[c+0:1]1:[c+0:2]:[c+0:3]:[c+0:4]:[c+0:5]:[c+0:6]:1.[H][H]>>[*:1]1[*:2][*:3][*:4][*:5][*:6]1",
        (1, 3),
        (1,),
    ),
    # Kolbe Carboxylation
    OperatorSmarts(
        "Kolbe Carboxylation",
        "[c+0:1]1:[c+0:2]([O+0H:7]):[c+0H:3]:[c+0:4]:[c+0:5]:[c+0:6]:1.[O+0:8]=[C+0:9]=[O+0:10]>>[*:1]1:[*:2]([*:7]):[*:3]([*:9](=[*:8])[*:10]):[*:4]:[*:5]:[*:6]:1",
        (1, 1),
        (1,),
    ),
    # Phenols Oxidation to Quinones
    # 1
    OperatorSmarts(
        "Phenols Oxidation to Quinones 1",
        "[c+0;!$(*~O):1]1:[c+0:2]([O+0H:7]):[c+0;!$(*~O):3]:[c+0;!$(*~O):4]:[c+0H:5]:[c+0;!$(*~O):6]:1.[O+0:8]=[O+0:9]>>[*:1]1[*:2](=[*:7])[*:3]=[*:4][*:5](=[*:8])[*:6]=1.[*:9]",
        (1, 1),
        (1, 1),
    ),
    # 2
    OperatorSmarts(
        "Phenols Oxidation to Quinones 2",
        "[c+0;!$(*~O):1]1:[c+0:2]([O+0H:7]):[c+0H;!$(*~O):3]:[c+0;!$(*~O):4]:[c+0:5]([*;!O+0:10]):[c+0;!$(*~O):6]:1.[O+0:8]=[O+0:9]>>[*:1]1[*:2](=[*:7])[*:3](=[*:8])[*:4]=[*:5]([*:10])[*:6]=1.[*:9]",
        (1, 1),
        (1, 1),
    ),
    # 3 & reverse
    OperatorSmarts(
        "Phenols Oxidation to Quinones 3",
        "[c+0;!$(*~O):1]1:[c+0:2]([O+0H:7]):[c+0:3]([O+0H:8]):[c+0;!$(*~O):4]:[c+0;!$(*~O):5]:[c+0;!$(*~O):6]:1.[O+0:9]=[O+0]>>[*:1]1[*:2](=[*:7])[*:3](=[*:8])[*:4]=[*:5][*:6]=1.[*:9]",
        (1, 0.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Phenols Oxidation to Quinones 3 Reverse",
        "[C+0:1]1[C+0:2](=[O+0:7])[C+0:3](=[O+0:8])[C+0:4]=[C+0:5][C+0:6]=1.[H][H]>>[*:1]1:[*:2]([*:7]):[*:3]([*:8]):[*:4]:[*:5]:[*:6]:1",
        (1, 1),
        (1,),
    ),
    # 4 & reverse
    OperatorSmarts(
        "Phenols Oxidation to Quinones 4",
        "[c+0;!$(*~O):1]1:[c+0:2]([O+0H:7]):[c+0;!$(*~O):3]:[c+0;!$(*~O):4]:[c+0:5]([O+0H:8]):[c+0;!$(*~O):6]:1.[O+0:9]=[O+0]>>[*:1]1[*:2](=[*:7])[*:3]=[*:4][*:5](=[*:8])[*:6]=1.[*:9]",
        (1, 0.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Phenols Oxidation to Quinones 4 Reverse",
        "[C+0:1]1[C+0:2](=[O+0:7])[C+0:3]=[C+0:4][C+0:5](=[O+0:8])[C+0:6]=1.[H][H]>>[*:1]1:[*:2]([*:7]):[*:3]:[*:4]:[*:5]([*:8]):[*:6]:1",
        (1, 1),
        (1,),
    ),
    # Quinones Addition
    OperatorSmarts(
        "Quinones Addition",
        "[C+0:1]1[C+0:2](=[O+0:7])[C+0:3]=[C+0H:4][C+0:5](=[O+0:8])[C+0:6]=1.[F,Cl,Br,I;+0H:9]>>[*:1]1:[*:2]([*:7]):[*:3]:[*:4]([*:9]):[*:5]([*:8]):[*:6]:1",
        (1, 1),
        (1,),
    ),
    # Naphthalene Oxidation
    # with CrO3
    OperatorSmarts(
        "Naphthalene Oxidation with CrO3",
        "[c+0:1]1:[c+0H:2]:[c+0:3]([c+0:9][c+0:10][c+0:11][c+0:12]2):[c+0:4]2:[c+0H:5]:[c+0:6]:1.[O+0:7]=[O+0:8]>>[*:1]1[*:2](=[*:7])[*:3](:[*:9]:[*:10]:[*:11]:[*:12]:2)=[*:4]2[*:5](=[O])[*:6]=1.[*:8]",
        (1, 1.5),
        (1, 1),
    ),
    # with V2O5 catalyst
    OperatorSmarts(
        "Naphthalene Oxidation with V2O5 Catalyst",
        "[c+0:1]1:[c+0:2]:[c+0:3]([c+0H:9][c+0H:10][c+0H][c+0H:12]2):[c+0:4]2:[c+0:5]:[c+0:6]:1.[O+0:7]=[O+0:8]>>[*:1]1:[*:2]:[*:3]([C:9](=[O])[O]):[*:4]([C:12](=[O])[O]):[*:5]:[*:6]:1.[*:7]=[C:10]=[O].[*:8]",
        (1, 4.5),
        (1, 2, 1),
    ),
    # Oxidation of Aromatic Alkanes
    # with -CH3
    OperatorSmarts(
        "Oxidation of Aromatic Alkanes Ar-CH3",
        "[c+0:1][C+0H3:2].[O+0:3]=[O+0:4]>>[*:1][*:2](=[*:3])[O].[*:4]",
        (1, 1.5),
        (1, 1),
    ),
    # with -CH2CH3
    OperatorSmarts(
        "Oxidation of Aromatic Alkanes Ar-CH2CH3",
        "[c+0:1][C+0H2:2][C+0H3:3].[O+0:4]=[O+0:5]>>[*:1][*:2](=[*:4])[O].[O]=[*:3]=[O].[*:5]",
        (1, 2.5),
        (1, 1, 2),
    ),
    # with -(CH3)CH3
    OperatorSmarts(
        "Oxidation of Aromatic Alkanes Ar-(CH3)CH3",
        "[c+0:1][C+0H:2]([C+0H3:3])[C+0H3].[O+0:4]=[O+0:5]>>[*:1][*:2](=[*:4])[O].[O]=[*:3]=[O].[*:5]",
        (1, 4.5),
        (1, 2, 3),
    ),
    # with long chain       forgot why used C+0!H3:4
    OperatorSmarts(
        "Oxidation of Aromatic Alkanes Ar-long_chain",
        "[c+0:1][C+0H2:2]!@[C+0H2:3][C+0!H3:4].[O+0:5]=[O+0:6]>>[*:1][*:2](=[*:5])[O].[*:4][*:3](=[O])[O].[*:6]",
        (1, 2.5),
        (1, 1, 1),
    ),
    # Ring-opening
    OperatorSmarts(
        "Oxidation of Aromatic Alkanes Ring-opening 1",
        "[c+0:1][C+0H2:2]@-[C+0H2:3].[O+0:4]=[O+0:5]>>([*:1][*:2](=[*:4])[O].[*:3](=[O])[O]).[*:5]",
        (1, 2.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Oxidation of Aromatic Alkanes Ring-opening 2",
        "[c+0:1][C+0H:2]@=[C+0H:3].[O+0:4]=[O+0:5]>>([*:1][*:2](=[*:4])[O].[*:3](=[*:5])[O])",
        (1, 2),
        (1,),
    ),
    # Oxidation of Aromatic Alkanes, 2-step
    # with -CH3
    OperatorSmarts(
        "Oxidation of Aromatic Alkanes Ar-CH3, 2-step",
        "([c+0:1][C+0H3:2].[c+0:5][C+0H3:6]).[O+0:3]=[O+0:4]>>([*:1][*:2](=[*:3])[O].[*:5][*:6](=[O])[O]).[*:4]",
        (1, 3),
        (1, 2),
        number_of_steps=2,
    ),
    # Oxidation of Aldehyde & Aromatic Alkanes
    # with -CH3
    OperatorSmarts(
        "Oxidation of Aldehyde & Aromatic Alkanes, 2-step",
        "([C+0!H0:1]=[O+0:2].[c+0:5][C+0H3:6]).[O+0:3]=[O+0:4]>>([*:1](=[*:2])[*:3].[*:5][*:6](=[O])[O]).[*:4]",
        (1, 2),
        (1, 1),
        number_of_steps=2,
    ),
    # Oxidation of Alcohol & Aromatic Alkanes
    # with -CH3
    OperatorSmarts(
        "Oxidation of Alcohol & Aromatic Alkanes, 2-step",
        "([C+0;H3,H2:1][O+0H:2].[c+0:5][C+0H3:6]).[O+0:3]=[O+0:4]>>([*:1](=[*:2])[*:3].[*:5][*:6](=[O])[O]).[*:4]",
        (1, 2.5),
        (1, 2),
        number_of_steps=2,
    ),
    # Halogenation of Aromatic Alkanes
    OperatorSmarts(
        "Halogenation of Aromatic Alkanes",
        "[c+0:1][C+0!H0:2].[F,Cl,Br,I;+0:3][F,Cl,Br,I;+0:4]>>[*:1][*:2][*:3].[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Electrophilic Aromatic Alkylation with Alkenes
    OperatorSmarts(
        "Electrophilic Aromatic Alkylation with Alkenes",
        "[c+0H:1].[C+0:2]=[C+0:3]>>[*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Electrophilic Aromatic Alkylation with Alkenes, Intramolecular",
        "([c+0H:1].[C+0:2]=[C+0:3])>>[*:1][*:2][*:3]",
        (1,),
        (1,),
    ),
    # Electrophilic Aromatic Alkylation with Alcohols
    OperatorSmarts(
        "Electrophilic Aromatic Alkylation with Alcohols",
        "[c+0H:1].[C+0;!$(*=O):2][O+0H:3]>>[*:1][*:2].[*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Electrophilic Aromatic Alkylation with Alcohols, Intramolecular",
        "([c+0H:1].[C+0:2][O+0H:3])>>[*:1][*:2].[*:3]",
        (1,),
        (1, 1),
    ),
    # Furan Carboxylation
    OperatorSmarts(
        "Furan Carboxylation",
        "[o+0:1]1:[c+0H:2]:[c+0:3]:[c+0:4]:[c+0:5]:1.[O+0:8]=[C+0:9]=[O+0:10]>>[*:1]1:[*:2]([*:9](=[*:8])[*:10]):[*:3]:[*:4]:[*:5]:1",
        (1, 1),
        (1,),
    ),
    # Furan Carboxylation, 2-step
    OperatorSmarts(
        "Furan Carboxylation, 2-step",
        "[o+0:1]1:[c+0H:2]:[c+0:3]:[c+0:4]:[c+0H:5]:1.[O+0:8]=[C+0:9]=[O+0:10]>>[*:1]1:[*:2]([*:9](=[*:8])[*:10]):[*:3]:[*:4]:[*:5]([C](=[O])[O]):1",
        (1, 2),
        (1,),
        number_of_steps=2,
    ),
    # Tautomerism Cyclohexadienone?
    # Friedel-Crafts Hydroxyalkylation
    OperatorSmarts(
        "Friedel-Crafts Hydroxyalkylation",
        "[c+0H:1].[CX3!$(*-O)+0:2]=[O+0:3]>>[*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Friedel-Crafts Hydroxyalkylation, Intramolecular",
        "([c+0H:1].[CX3!$(*-O)+0:2]=[O+0:3])>>[*:1][*:2][*:3]",
        (1,),
        (1,),
    ),
    OperatorSmarts(
        "Friedel-Crafts Hydroxyalkylation(BPA Synthesis)",
        "[c+0H:1].[c+0H:2].[CX3!$(*-O)+0:3]=[O+0:4]>>[*:1][*:3][*:2].[*:4]",
        (1, 1, 1),
        (1, 1),
    ),
    # Carbonylation #######################################################
    # Hydroformylation
    OperatorSmarts(
        "Hydroformylation",
        "[C+0:1]=[C+0:2].[C-1]#[O+1].[H][H]>>[*:1][*:2][C]=[O]",
        (1, 1, 1),
        (1,),
    ),
    # Alcohol carbonylation to carboxylic acid, Cativa process
    OperatorSmarts(
        "Cativa Process",
        "[C!$(*=O)+0:1][O+0H:2].[C-1]#[O+1]>>[*:1][C](=[O])[*:2]",
        (1, 1),
        (1,),
    ),
    # Oxidative Carbonylation, Intramolecular?
    OperatorSmarts(
        "Oxidative Carbonylation 1",
        "[C!$(*=O)+0:1][O+0H:2].[C!$(*=O)+0:5][O+0H:6].[C-1]#[O+1].[O+0:7]=[O+0]>>[*:1][*:2][C](=[O])[*:6][*:5].[*:7]",
        (1, 1, 1, 0.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Oxidative Carbonylation 2",
        "[C!$(*=O)+0:1][O+0H:2].[C!$(*=O)+0:5][O+0H:6].[C-1]#[O+1].[O+0:7]=[O+0]>>[*:1][*:2][C](=[O])[C](=[O])[*:6][*:5].[*:7]",
        (1, 1, 2, 0.5),
        (1, 1),
    ),
    # Hydrocarboxylation
    OperatorSmarts(
        "Hydrocarboxylation 1, Hydroesterification(Carboalkoxylation)",
        "[C+0:1]=[C+0:2].[C-1]#[O+1].[O+0!H0;!$(*-C=O):5]>>[*:1][*:2][C](=[O])[*:5]",
        (1, 1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Hydrocarboxylation 2",
        "[C+0:1]#[C+0:2].[C-1]#[O+1].[O+0H2:5]>>[*:1]=[*:2][C](=[O])[*:5]",
        (1, 1, 1),
        (1,),
    ),
    # Carbohydrates #######################################################
    # Oxidation of Pyranosides
    OperatorSmarts(
        "Oxidation of Pyranosides",
        "[C+0:9][C+0H:1]([O+0H:2])@[C+0H:3]([O+0H:4])[C+0H:5]([O+0H:6])[C+0:10].[O+0:7]=[O+0:8]>>([*:9][*:1]=[*:2].[*:10][*:5]=[*:6]).[*:7][*:3]=[*:4].[*:8]",
        (1, 1),
        (1, 1, 1),
    ),
    # Amines, Amides ##############################################################
    # Dehydration of Amides to Form Nitriles
    OperatorSmarts(
        "Dehydration of Amides",
        "[C+0:1](=[O+0:2])[N+0H2:3]>>[*:1]#[*:3].[*:2]",
        (1,),
        (1, 1),
    ),
    # Ketone Reductive Amination
    OperatorSmarts(
        "Ketone Reductive Amination",
        "[CX3!$(*[O,S,N])+0:2](=[O+0:3]).[N+0X3!H0:6].[H][H]>>[*:2][*:6].[*:3]",
        (1, 1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Ketone Reductive Amination, Intramolecular",
        "([CX3!$(*[O,S,N])+0:2](=[O+0:3]).[N+0X3!H0:6]).[H][H]>>[*:2][*:6].[*:3]",
        (1, 1),
        (1, 1),
    ),
    # Alkylation or Acylation of Amines
    OperatorSmarts(
        "Alkylation or Acylation of Amines",
        "[C,c;+0:1][F,Cl,Br,I;+0:2].[NX3!H0,nH;+0:3]>>[*:1][*:3].[*:2]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Alkylation or Acylation of Amines, Intramolecular",
        "([C,c;+0:1][F,Cl,Br,I;+0:2].[NX3!H0,nH;+0:3])>>[*:1][*:3].[*:2]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Alkylation of Tertiary Amines
    OperatorSmarts(
        "Alkylation of Tertiary Amines",
        "[C,c;+0:1][F,Cl,Br,I;+0:2].[NX3H0+0:3]>>[*:1][*+1:3].[*-1:2]",
        (1, 1),
        (1, 1),
    ),
    # Amine Alkylation with Alcohols or Primary Amines
    OperatorSmarts(
        "Amine Alkylation with Alcohols or Primary Amines",
        "[C,c;!$(*=O)+0:1][OH,NH2;+0:2].[NX3!H0,nH;+0:3]>>[*:1][*:3].[*:2]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Amine Alkylation with Alcohols or Primary Amines, Intramolecular",
        "([C,c;!$(*=O)+0:1][OH,NH2;+0:2].[NX3!H0,nH;+0:3])>>[*:1][*:3].[*:2]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Synthesis of Amides with Carboxylic Acid
    OperatorSmarts(
        "Synthesis of Amides with Carboxylic Acid",
        "[CX3+0:2](=[O+0:3])[O+0H:4].[N+0X3!H0:5]>>[*:2](=[*:3])[*:5].[*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Synthesis of Amides with Carboxylic Acid, Intramolecular",
        "([CX3+0:2](=[O+0:3])[O+0H:4].[N+0X3!H0:5])>>[*:2](=[*:3])[*:5].[*:4]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Synthesis of Amides with Acid Anhydrides or Esters
    OperatorSmarts(
        "Synthesis of Amides with Acid Anhydrides or Esters",
        "[CX3+0:2](=[O+0:3])-!@[O+0:4][C,c;+0:6].[N+0X3!H0:5]>>[*:2](=[*:3])[*:5].[*:4][*:6]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Synthesis of Amides with Acid Anhydrides or Esters, Intramolecular",
        "([CX3+0:2](=[O+0:3])-!@[O+0:4][C,c;+0:6].[N+0X3!H0:5])>>[*:2](=[*:3])[*:5].[*:4][*:6]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    # Hydrolysis of Amides
    OperatorSmarts(
        "Hydrolysis of Amides",
        "[CX3!$(*[OH,SH])+0:2](=[O+0:3])-!@[N+0X3:5].[O+0H2:4]>>[*:2](=[*:3])[*:4].[*:5]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Hydrolysis of Amides, Intramolecular",
        "[CX3!$(*[OH,SH])+0:2](=[O+0:3])-@[N+0X3:5].[O+0H2:4]>>([*:2](=[*:3])[*:4].[*:5])",
        (1, 1),
        (1,),
    ),
    # Hofmann Elimination
    OperatorSmarts(
        "Hofmann Elimination R-NH2",
        "[C+0X4!H0:1][CX4+0:2][N+0H2:3].[C+0:4][F,Cl,Br,I;+0:5]>>[*:1]=[*:2].[*:3]([*:4])([*:4])[*:4].[*:5]",
        (1, 3),
        (1, 1, 3),
    ),
    OperatorSmarts(
        "Hofmann Elimination R-NHR",
        "[C+0X4!H0:1][CX4+0:2]-!@[N+0H:3].[C+0:4][F,Cl,Br,I;+0:5]>>[*:1]=[*:2].[*:3]([*:4])[*:4].[*:5]",
        (1, 2),
        (1, 1, 2),
    ),
    OperatorSmarts(
        "Hofmann Elimination R-NHR Ringopening",
        "[C+0X4!H0:1][CX4+0:2]-@[N+0H:3].[C+0:4][F,Cl,Br,I;+0:5]>>([*:1]=[*:2].[*:3]([*:4])[*:4]).[*:5]",
        (1, 2),
        (1, 2),
    ),
    OperatorSmarts(
        "Hofmann Elimination R-NRR",
        "[C+0X4!H0:1][CX4+0:2]-!@[N+0X3H0:3].[C+0:4][F,Cl,Br,I;+0:5]>>[*:1]=[*:2].[*:3][*:4].[*:5]",
        (1, 1),
        (1, 1, 1),
    ),
    OperatorSmarts(
        "Hofmann Elimination R-NRR Ringopening",
        "[C+0X4!H0:1][CX4+0:2]-@[N+0X3H0:3].[C+0:4][F,Cl,Br,I;+0:5]>>([*:1]=[*:2].[*:3][*:4]).[*:5]",
        (1, 1),
        (1, 1),
    ),
    # Hofmann Rearrangement
    OperatorSmarts(
        "Hofmann Rearrangement",
        "[C,c;+0:1][C+0:2](=[O+0:3])[N+0H2:5].[F,Cl,Br,I;+0:6][F,Cl,Br,I;+0].[O+0H2:4]>>[*:1][*:5].[*:3]=[*:2]=[*:4].[*:6]",
        (1, 1, 1),
        (1, 1, 2),
    ),
    OperatorSmarts(
        "Hofmann Rearrangement with Alcohols",
        "[C,c;+0:1][C+0:2](=[O+0:3])[N+0H2:5].[F,Cl,Br,I;+0:6][F,Cl,Br,I;+0].[C+0!$(*=[O,S]):7][O+0H:4]>>[*:1][*:5][*:2](=[*:3])[*:4][*:7].[*:6]",
        (1, 1, 1),
        (1, 2),
    ),
    # Hydroamination of Alkenes
    OperatorSmarts(
        "Hydroamination of Alkenes",
        "[C+0:1]=[C+0:2].[N+0X3!H0:3]>>[*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Hydroamination of Alkenes, Intramolecular",
        "([C+0:1]=[C+0:2].[N+0X3!H0:3])>>[*:1][*:2][*:3]",
        (1,),
        (1,),
        ring_issue=True,
    ),
    # Hydroamination of Alkynes
    OperatorSmarts(
        "Hydroamination of Alkynes",
        "[C+0:1]#[C+0:2].[N+0X3!H0:3]>>[*:1]=[*:2][*:3]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Hydroamination of Alkynes, Intramolecular",
        "([C+0:1]#[C+0:2].[N+0X3!H0:3])>>[*:1]=[*:2][*:3]",
        (1,),
        (1,),
        ring_issue=True,
    ),
    # Hydroamination of Dienes
    OperatorSmarts(
        "Hydroamination of Dienes",
        "[C+0:1]=[C+0:2][C+0:3]=[C+0:4].[N+0X3!H0:5]>>[*:1][*:2]=[*:3][*:4][*:5]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Hydroamination of Dienes, Intramolecular",
        "([C+0:1]=[C+0:2][C+0:3]=[C+0:4].[N+0X3!H0:5])>>[*:1][*:2]=[*:3][*:4][*:5]",
        (1,),
        (1,),
        ring_issue=True,
    ),
    # Ring Opening of Epoxides by Amines
    OperatorSmarts(
        "Ring Opening of Epoxides by Amines",
        "[C+0:1]1[C+0:2][O+0:3]1.[N+0X3!H0:4]>>[*:4][*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    # Aromatic Nitrosation
    OperatorSmarts(
        "Aromatic Nitrosation",
        "[OH,N$(*([#6])([#6])[#6]);+0:10][c+0:1]1[c+0:2][c+0:3][c+0H:4][c+0:5][c+0:6]1.[O+0H:7][N+0:8]=[O+0:9]>>[*:10][*:1]1:[*:2]:[*:3]:[*:4]([*:8]=[*:9]):[*:5]:[*:6]:1.[*:7]",
        (1, 1),
        (1, 1),
    ),
    # Secondary Amines with Nitrous Acid
    OperatorSmarts(
        "Secondary Amines with Nitrous Acid",
        "[*+0:1][N+0H:2][*+0:3].[O+0H:4][N+0:5]=[O+0:6]>>[*:1][*:2]([*:5]=[*:6])[*:3].[*:4]",
        (1, 1),
        (1, 1),
    ),
    # Primary Amines with Nitrous Acid to Alcohols
    OperatorSmarts(
        "Primary Amines with Nitrous Acid to Alcohols",
        "[C,c;+0:1][N+0H2:2].[O+0H:3][N+0:4]=[O+0:5]>>[*:1][*:3].[*:5].[*:2]#[*:4]",
        (1, 1),
        (1, 1, 1),
    ),
    # Primary Amines with Nitrous Acid to Halides
    OperatorSmarts(
        "Primary Amines with Nitrous Acid to Halides",
        "[C,c;+0:1][N+0H2:2].[F,Cl,Br,I;+0H:3].[O+0H][N+0:5]=[O+0:6]>>[*:1][*:3].[*:6].[*:2]#[*:5]",
        (1, 1, 1),
        (1, 2, 1),
    ),
    # Primary Amines with Nitrous Acid to Alkenes
    OperatorSmarts(
        "Primary Amines with Nitrous Acid to Alkenes",
        "[C+0!H0:6][C+0:1][N+0H2:2].[O+0H][N+0:4]=[O+0:5]>>[*:6]=[*:1].[*:5].[*:2]#[*:4]",
        (1, 1),
        (1, 2, 1),
    ),
    # Tiffeneau–Demjanov Rearrangement
    OperatorSmarts(
        "Tiffeneau–Demjanov Rearrangement",
        "[O+0H:1][C+0:2]([C+0:3])([C+0:4])[C+0:5][N+0H2:6].[O+0H][N+0:7]=[O+0:8]>>[*:1]=[*:2]([*:3])[*:5][*:4].[*:8].[*:6]#[*:7]",
        (1, 1),
        (1, 2, 1),
        ring_issue=True,
    ),
    # Sandmeyer Cyanation
    OperatorSmarts(
        "Sandmeyer Cyanation",
        "[c+0:1][N+0H2:2].[C+0H:3]#[N+0:4].[O+0H][N+0:5]=[O+0:6]>>[*:1][*:3]#[*:4].[*:6].[*:2]#[*:5]",
        (1, 1, 1),
        (1, 2, 1),
    ),
    # Reduction of Primary Aromatic Amines
    OperatorSmarts(
        "Reduction of Primary Aromatic Amines",
        "[c+0:1][N+0H2:2].[H][H].[O+0H][N+0:3]=[O+0:4]>>[*:1].[*:4].[*:2]#[*:3]",
        (1, 1, 1),
        (1, 2, 1),
    ),
    # Cope Elimination
    OperatorSmarts(
        "Cope Elimination",
        "[C+0!H0:1][C+0:2]!@[N+0:3]([C+0:4])[C+0:5].[O+0]=[O+0:6]>>[*:1]=[*:2].[*:6][*:3]([*:4])[*:5]",
        (1, 0.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Cope Elimination, Intramolecular",
        "[C+0!H0:1][C+0:2]@[N+0:3]([C+0:4])[C+0:5].[O+0]=[O+0:6]>>([*:1]=[*:2].[*:6][*:3]([*:4])[*:5])",
        (1, 0.5),
        (1,),
    ),
    # Hemiaminal Formation
    OperatorSmarts(
        "Hemiaminal Formation",
        "[CX3!$(*O)+0:1]=[O+0:2].[N,n;+0X3!H0:3]>>[*:1]([*:2])[*:3]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Hemiaminal Formation, Reverse",
        "[CX4!$(*(O)O)+0:1]([O+0H:2])-!@[N+0X3:3]>>[*:1]=[*:2].[*:3]",
        (1,),
        (1, 1),
        kekulize_flag=True,
    ),
    # Hemiaminal Dehydration
    OperatorSmarts(
        "Hemiaminal Dehydration",
        "[C+0:1]([O+0H:2])[N+0X3!H0:3]>>[*:1]=[*:3].[*:2]",
        (1,),
        (1, 1),
    ),
    OperatorSmarts(
        "Hemiaminal Dehydration, Reverse",
        "[C+0:1]=[N+0:3].[O+0H2:2]>>[*:1]([*:2])[*:3]",
        (1, 1),
        (1,),
    ),
    # Nitroso, Nitro Compounds ###############################################
    # Reduction of Nitroso
    OperatorSmarts(
        "Reduction of Nitroso",
        "[C,c;+0:1][N+0:2]=[O+0:3].[H][H]>>[*:1][*:2].[*:3]",
        (1, 2),
        (1, 1),
    ),
    # Reduction of Nitro Compounds
    OperatorSmarts(
        "Reduction of Nitro Compounds",
        "[C,c;+0:1][N+1](=[O+0:3])[O-1H0].[H][H]>>[*:1][N].[*:3]",
        (1, 3),
        (1, 2),
    ),
    # Henry Reaction
    OperatorSmarts(
        "Henry Reaction",
        "[C+0!H0:1][N+1:2](=[O+0:3])[OX1-1:4].[CX3!$(*([OH])[OH])+0:5]=[O+0:6]>>[*:6][*:5][*:1][*:2](=[*:3])[*:4]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Henry Reaction, Intramolecular",
        "([C+0!H0:1][N+1:2](=[O+0:3])[OX1-1:4].[CX3!$(*([OH])[OH])+0:5]=[O+0:6])>>[*:6][*:5][*:1][*:2](=[*:3])[*:4]",
        (1,),
        (1,),
        ring_issue=True,
    ),
    # Benzene Nitration
    OperatorSmarts(
        "Benzene Nitration",
        "[c+0H:1].[O+0H:2][N+1:3](=[O+0:4])[OX1-1:5]>>[*:1][*:3](=[*:4])[*:5].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Oxidation of Nitroso
    OperatorSmarts(
        "Oxidation of Nitroso",
        "[*+0:1][N+0]=[O+0].[O+0]=[O+0]>>[*:1][N+1](=[O])[O-1]",
        (1, 0.5),
        (1,),
    ),
    # Esterification, Acid Anhydride Formation with Nitric Acid, Nitrate Esters Hydrolysis
    OperatorSmarts(
        "Esterification, Acid Anhydride Formation with Nitric Acid",
        "[C,c;+0:1][O+0H:2].[O+0H:3][N+1:4](=[O+0:5])[OX1-1:6]>>[*:1][*:2][*:4](=[*:5])[*:6].[*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Nitrate Esters Hydrolysis",
        "[C,c;+0:1][O+0:2][N+1:4](=[O+0:5])[OX1-1:6].[O+0H2:3]>>[*:1][*:2].[*:3][*:4](=[*:5])[*:6]",
        (1, 1),
        (1, 1),
    ),
    # Imines ###############################################
    # Imines from Aldehydes and Ketones
    OperatorSmarts(
        "Imines from Aldehydes and Ketones",
        "[CX3!$(*[OH])+0:1]=[O+0:2].[N+0X3;H2,H3:3]>>[*:1]=[*:3].[*:2]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Imines from Aldehydes and Ketones, Intramolecular",
        "([CX3!$(*[OH])+0:1]=[O+0:2].[N+0X3;H2,H3:3])>>[*:1]=[*:3].[*:2]",
        (1,),
        (1, 1),
        ring_issue=True,
    ),
    OperatorSmarts(
        "Imines from Aldehydes and Ketones Reverse",
        "[CX3!$(*[OH])+0:1]=!@[N+0:3].[O+0H2:2]>>[*:1]=[*:2].[*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Imines from Aldehydes and Ketones, Intramolecular Reverse",
        "[CX3!$(*[OH])+0:1]=@[N+0:3].[O+0H2:2]>>([*:1]=[*:2].[*:3])",
        (1, 1),
        (1,),
    ),
    # Treatment of Aldehydes and Ketones with a Secondary Amine
    OperatorSmarts(
        "Treatment of Aldehydes and Ketones with a Secondary Amine",
        "[C,N;+0!H0:4][CX3!$(*[OH])+0:1]=[O+0:2].[N+0X3H:3]>>[*:4]=[*:1][*:3].[*:2]",
        (1, 1),
        (1, 1),
        kekulize_flag=True,
    ),
    OperatorSmarts(
        "Treatment of Aldehydes and Ketones with a Secondary Amine, Intramolecular",
        "([C,N;+0!H0:4][CX3!$(*[OH])+0:1]=[O+0:2].[N+0X3H:3])>>[*:4]=[*:1][*:3].[*:2]",
        (1,),
        (1, 1),
        ring_issue=True,
        kekulize_flag=True,
    ),
    OperatorSmarts(
        "Treatment of Aldehydes and Ketones with a Secondary Amine Reverse",
        "[C,N;+0:4]=[CX3!$(*[OH])+0:1]-!@[NX3+0:3].[O+0H2:2]>>[*:4][*:1]=[*:2].[*:3]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Treatment of Aldehydes and Ketones with a Secondary Amine, Intramolecular Reverse",
        "[C,N;+0:4]=[CX3!$(*[OH])+0:1]-@[NX3+0:3].[O+0H2:2]>>([*:4][*:1]=[*:2].[*:3])",
        (1, 1),
        (1,),
    ),
    # Oxime-Nitroso Tautomerization
    OperatorSmarts(
        "Oxime-Nitroso Tautomerization",
        "[C+0:1]=[N+0:2][O+0H:3]>>[*:1][*:2]=[*:3]",
        (1,),
        (1,),
    ),
    OperatorSmarts(
        "Keto-enol Tautomerization Reverse",
        "[C+0!H0:1][N+0:2]=[O+0:3]>>[*:1]=[*:2][*:3]",
        (1,),
        (1,),
    ),
    # Reduction of Imines
    OperatorSmarts(
        "Reduction of Imines", "[C+0:1]=[N+0:2].[H][H]>>[*:1][*:2]", (1, 1), (1,)
    ),
    # Beckmann Rearrangement
    OperatorSmarts(
        "Beckmann Rearrangement with Ketoximes",
        "[C,c;+0:1][C+0:2](=[N+0:3][O+0H:4])[C,c;+0:5]>>[*:1][*:2](=[*:4])[*:3][*:5]",
        (1,),
        (1,),
    ),
    OperatorSmarts(
        "Beckmann Rearrangement with Aldoximes",
        "[C,c;+0:1][C+0H:2](=[N+0:3][O+0H:4])>>[*:1][*:2]#[*:3].[*:4]",
        (1,),
        (1, 1),
    ),
    # Beckmann Rearrangement from Ketones and Aldehydes
    OperatorSmarts(
        "Beckmann Rearrangement from Ketones",
        "[C,c;+0:1][C+0:2](=[O+0:3])[C,c;+0:4].[N+0H2:5][O+0H:6]>>[*:1][*:2](=[*:3])[*:5][*:4].[*:6]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Beckmann Rearrangement from Aldehydes",
        "[C,c;+0:1][C+0H:2]=[O+0:3].[N+0H2:4][O+0H]>>[*:1][*:2]#[*:4].[*:3]",
        (1, 1),
        (1, 2),
    ),
    # Beckmann Fragmentation
    OperatorSmarts(
        "Beckmann Fragmentation",
        "[C,c;+0:1][C+0:2](=[N+0:3][O+0H:4])-!@[C+0:5]([C+0:6])([C+0:7])[*+0!H0:8]>>[*:1][*:2]#[*:3].[*:5]([*:6])([*:7])=[*:8].[*:4]",
        (1,),
        (1, 1, 1),
    ),
    OperatorSmarts(
        "Beckmann Fragmentation, Intramolecular",
        "[C,c;+0:1][C+0:2](=[N+0:3][O+0H:4])-@[C+0:5]([C+0:6])([C+0:7])[*+0!H0:8]>>([*:1][*:2]#[*:3].[*:5]([*:6])([*:7])=[*:8]).[*:4]",
        (1,),
        (1, 1),
    ),
    # Semmler–Wolff Reaction
    OperatorSmarts(
        "Semmler–Wolff Reaction",
        "[C+0!H0:1]1[C+0!H0:2][C+0!H0:3][C+0:4]=[C+0:5][C+0:6]1=[N+0:7][O+0H:8]>>[*:1]1[*:2]=[*:3][*:4]=[*:5][*:6]=1[*:7].[*:8]",
        (1,),
        (1, 1),
        kekulize_flag=True,
    ),
    # Transimination
    OperatorSmarts(
        "Transimination",
        "[C,c;+0:1]=!@[N!$(*[O,S])+0:2].[N!$(*[O,S]);H2,H3;+0:3]>>[*:1]=[*:3].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Imine Metathesis
    OperatorSmarts(
        "Imine Metathesis",
        "[C,c;+0:1]=!@[N!$(*[O,S])+0:2].[C,c;+0:3]=!@[N!$(*[O,S])+0:4]>>[*:1]=[*:4].[*:3]=[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Wolff–Kishner Reduction
    OperatorSmarts(
        "Wolff–Kishner Reduction",
        "[C,c;+0:1][C!$(*[O,S])+0:2](=[N+0:3][N+0H2:4])>>[*:1][*:2].[*:3]#[*:4]",
        (1,),
        (1, 1),
    ),
    # Nitriles, Amino Acids ###############################################
    # Cyanation (Kolbe Nitrile Synthesis)
    OperatorSmarts(
        "Cyanation (Kolbe Nitrile Synthesis)",
        "[C,c+0:1][F,Cl,Br,I;+0:2].[C+0H:3]#[N+0:4]>>[*:1][*:3]#[*:4].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Cyanation of Ketones or Aldehydes (Cyanohydrin Reaction)
    OperatorSmarts(
        "Cyanation of Ketones or Aldehydes (Cyanohydrin Reaction)",
        "[CX3!$(*[O,S,N])+0:1]=[O+0:2].[C+0H:3]#[N+0:4]>>[*:1]([*:3]#[*:4])[*:2]",
        (1, 1),
        (1,),
    ),
    # Alkyl Cyanides from Haloalkanes
    OperatorSmarts(
        "Alkyl Cyanides from Haloalkanes",
        "[C+0:1][F,Cl,Br,I;+0:2].[C+0H:3]#[N+0:4]>>[*:1][*:3]#[*:4].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Nitrile Hydrogenation to Amines
    OperatorSmarts(
        "Nitrile Hydrogenation to Amines",
        "[C+0:1]#[N+0:2].[H][H]>>[*:1][*:2]",
        (1, 2),
        (1,),
    ),
    # Nitrile Hydrogenation to Aldehydes
    OperatorSmarts(
        "Nitrile Hydrogenation to Aldehydes",
        "[C+0:1]#[N+0:2].[O+0H2:3].[H][H]>>[*:1]=[*:3].[*:2]",
        (1, 1, 1),
        (1, 1),
    ),
    # Nitrile Hydrogenation to Imines
    OperatorSmarts(
        "Nitrile Hydrogenation to Imines",
        "[C+0:1]#[N+0:2].[H][H]>>[*:1]=[*:2]",
        (1, 1),
        (1,),
    ),
    # Nitrile Hydrogenation to Secondary Amines
    OperatorSmarts(
        "Nitrile Hydrogenation to Secondary Amines",
        "[C+0:1]#[N+0:2].[C+0:3]#[N+0:4].[H][H]>>[*:1][*:2][*:3].[*:4]",
        (1, 1, 4),
        (1, 1),
    ),
    OperatorSmarts(
        "Nitrile Hydrogenation to Secondary Amines, Intramolecular",
        "([C+0:1]#[N+0:2].[C+0:3]#[N+0:4]).[H][H]>>[*:1][*:2][*:3].[*:4]",
        (1, 4),
        (1, 1),
    ),
    # Nitrile Hydrogenation to Tertiary Amines
    OperatorSmarts(
        "Nitrile Hydrogenation to Tertiary Amines",
        "[C+0:1]#[N+0:2].[C+0:3]#[N+0:4].[C+0:5]#[N+0].[H][H]>>[*:1][*:2]([*:3])[*:5].[*:4]",
        (1, 1, 1, 6),
        (1, 2),
    ),
    # Hydrolysis of Nitriles
    OperatorSmarts(
        "Hydrolysis of Nitriles",
        "[C+0:1]#[N+0:2].[O+0H2:3]>>[*:1](=[O])[*:3].[*:2]",
        (1, 2),
        (1, 1),
    ),
    # Nitriles to Ketones with Grignard Reagents
    # !NOT BALANCED! Reagents and by-products ignored for simplicity sake, needs to add a correction term for enthalpy calculation
    # Enthalpy correction from DFT: Mg -> MgBrOH, delta H = -93.89 kcal/mol, gas phase 298.15 K
    OperatorSmarts(
        "Nitriles to Ketones with Grignard Reagents",
        "[C+0:1]#[N+0:2].[c,C!$(*~O);!$(*#*);+0:3][F,Cl,Br,I;+0].[O+0H2:4]>>[*:1](=[*:4])[*:3].[*:2]",
        (1, 1, 2),
        (1, 1),
        enthalpy_correction=-93.89,
    ),
    # Partial Hydrolysis of Nitriles
    OperatorSmarts(
        "Partial Hydrolysis of Nitriles",
        "[C+0:1]#[N+0:2].[O+0H2:3]>>[*:1](=[*:3])[*:2]",
        (1, 1),
        (1,),
    ),
    # Epoxides Ring Opening with Cyanides
    OperatorSmarts(
        "Epoxides Ring Opening with Cyanides",
        "[C+0:1]1[C+0:2][O+0:3]1.[C+0H:4]#[N+0:5]>>[*:5]#[*:4][*:1][*:2][*:3]",
        (1, 1),
        (1,),
    ),
    # Azides ###############################################
    # Azide-Alkyne, Azide-Nitrile Cycloaddition
    OperatorSmarts(
        "Azide-Alkyne, Azide-Nitrile Cycloaddition",
        "[N+0:1]=[N+1]=[N-1].[C,N;+0:4]#[C+0:5]>>[*:1]1[N]=[N][*:4]=[*:5]1",
        (1, 1),
        (1,),
    ),
    # Azide Salts as Nucleophiles
    OperatorSmarts(
        "Azide Salts as Nucleophiles",
        "[C+0:1][F,Cl,Br,I,OH;+0:2].[N+0H:3]=[N+1:4]=[N-1:5]>>[*:1][*:3]=[*:4]=[*:5].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Reduction of Azides
    OperatorSmarts(
        "Reduction of Azides",
        "[N+0:1]=[N+1]=[N-1].[H][H]>>[*:1].[N]#[N]",
        (1, 1),
        (1, 1),
    ),
    # Curtius Rearrangement
    OperatorSmarts(
        "Curtius Rearrangement",
        "[C,c;+0:1][C+0:2](=[O+0:3])[N+0:4]=[N+1]=[N-1]>>[*:1][*:4]=[*:2]=[*:3].[N]#[N]",
        (1,),
        (1, 1),
    ),
    # Diazo ###############################################
    # Esters from Diazoes and Carboxylic Acids
    OperatorSmarts(
        "Esters from Diazoes and Carboxylic Acids",
        "[CX3!$(*([OH])[OH])+0:1](=[O+0:2])[O+0H:3].[CX3+0:4]=[N+1]=[N-1]>>[*:1](=[*:2])[*:3][*:4].[N]#[N]",
        (1, 1),
        (1, 1),
    ),
    # Isocyanates, Isothiocyanate ###############################################
    # Isocyanates, Isothiocyanate with Nucleophiles
    OperatorSmarts(
        "Isocyanates, Isothiocyanate with Nucleophiles",
        "[C,c;+0:1][N+0:2]=[C+0:3]=[O,S;+0:4].[NX3!H0,O!H0;+0:5]>>[*:1][*:2][*:3](=[*:4])[*:5]",
        (1, 1),
        (1,),
    ),
    # Isocyanates with Water
    OperatorSmarts(
        "Isocyanates with Water",
        "[C,c;+0:1][N+0:2]=[C+0:3]=[O+0:4].[OH2+0:5]>>[*:1][*:2].[*:3](=[*:4])=[*:5]",
        (1, 1),
        (1, 1),
    ),
    # Isocyanates with Carboxylic Acid
    OperatorSmarts(
        "Isocyanates with Carboxylic Acid",
        "[C,c;+0:1][N+0:2]=[C+0:3]=[O+0:4].[O+0H:5][C+0:6](=[O+0:7])[C,c;+0:8]>>[*:1][*:2][*:3](=[*:4])[*:8].[*:5]=[*:6]=[*:7]",
        (1, 1),
        (1, 1),
    ),
    # Synthesis of Isocyanates, Isothiocyanates
    OperatorSmarts(
        "Synthesis of Isocyanates, Isothiocyanates",
        "[C,c;+0:1][N+0H2:2].[Cl+0][C+0:4](=[O,SX1;+0:5])[Cl+0:6]>>[*:1][*:2]=[*:4]=[*:5].[*:6]",
        (1, 1),
        (1, 2),
    ),
    # Synthesis of Isothiocyanates
    OperatorSmarts(
        "Synthesis of Isothiocyanates",
        "[C,c;+0:1][N+0H2:2].[S+0:3]=[C+0:4]=[S+0:5]>>[*:1][*:2]=[*:4]=[*:5].[*:3]",
        (1, 1),
        (1, 1),
    ),
    # Thiols, Thioethers ###############################################
    # Hydrosulfide Anion Substitution with Alkyl Halides
    OperatorSmarts(
        "Hydrosulfide Anion Substitution with Alkyl Halides",
        "[C,c;+0:1][F,Cl,Br,I;+0:2].[S!H0X2+0:3]>>[*:1][*:3].[*:2]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Hydrosulfide Anion Substitution with Alkyl Halides, Intramolecular",
        "([C,c;+0:1][F,Cl,Br,I;+0:2].[S!H0X2+0:3])>>[*:1][*:3].[*:2]",
        (1,),
        (1, 1),
    ),
    # Preparation of Sulfides
    OperatorSmarts(
        "Preparation of Sulfides",
        "[C,c;+0:1][F,Cl,Br,I;+0:2].[C,c;+0:3][F,Cl,Br,I;+0].[SH2X2+0:5]>>[*:1][*:5][*:3].[*:2]",
        (1, 1, 1),
        (1, 2),
    ),
    OperatorSmarts(
        "Preparation of Sulfides, Intramolecular",
        "([C,c;+0:1][F,Cl,Br,I;+0:2].[C,c;+0:3][F,Cl,Br,I;+0]).[SH2X2+0:5]>>[*:1][*:5][*:3].[*:2]",
        (1, 1),
        (1, 2),
    ),
    # Thiols from Thiourea
    OperatorSmarts(
        "Thiols from Thiourea",
        "[C,c;+0:1][F,Cl,Br,I;+0:2].[SX1+0:3]=[C+0:4]([N+0H2:5])[N+0H2:6].[O+0H2:7]>>[*:1][*:3].[*:7]=[*:4]([*:5])[*:6].[*:2]",
        (1, 1, 1),
        (1, 1, 1),
    ),
    # Thioacetals from Thiols and Aldehydes or Ketones
    OperatorSmarts(
        "Thioacetals from Thiols and Aldehydes or Ketones",
        "[C,c;+0:1][SX2+0H:2].[C,c;+0:3][SX2+0H:4].[C!$(*O)X3+0:5]=[O+0:6]>>[*:1][*:2][*:5][*:4][*:3].[*:6]",
        (1, 1, 1),
        (1, 1),
    ),
    # Cyclic Thioacetals from Dithiols
    OperatorSmarts(
        "Cyclic Thioacetals from Dithiols",
        "([C,c;+0:1][SX2+0H:2].[C,c;+0:3][SX2+0H:4]).[C!$(*O)X3+0:5]=[O+0:6]>>[*:1][*:2][*:5][*:4][*:3].[*:6]",
        (1, 1),
        (1, 1),
    ),
    # Oxidation of Thiols, Reduction of Disulfides
    OperatorSmarts(
        "Oxidation of Thiols",
        "[C,c;+0:1][SX2+0H:2].[C,c;+0:3][SX2+0H:4].[O+0:5]=[O+0]>>[*:1][*:2][*:4][*:3].[*:5]",
        (1, 1, 0.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Oxidation of Thiols, Intramolecular",
        "([C,c;+0:1][SX2+0H:2].[C,c;+0:3][SX2+0H:4]).[O+0:5]=[O+0]>>[*:1][*:2][*:4][*:3].[*:5]",
        (1, 0.5),
        (1, 1),
    ),
    OperatorSmarts(
        "Reduction of Disulfides",
        "[C,c;+0:1][SX2+0:2]-!@[SX2+0:4][C,c;+0:3].[H][H]>>[*:1][*:2].[*:3][*:4]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Reduction of Disulfides, Intramolecular",
        "[C,c;+0:1][SX2+0:2]-@[SX2+0:4][C,c;+0:3].[H][H]>>([*:1][*:2].[*:3][*:4])",
        (1, 1),
        (1,),
    ),
    # Sulfoxides ###############################################
    # Oxidation of Sulfides, Sulfoxides
    OperatorSmarts(
        "Oxidation of Sulfides",
        "[*+0:1][SX2!$(*[OH])+0:2][*+0:3].[O+0:4]=[O+0]>>[*:1][*:2](=[*:4])[*:3]",
        (1, 0.5),
        (1,),
    ),
    OperatorSmarts(
        "Oxidation of Sulfoxides",
        "[*+0:1][SX3+0:2](=[O+0:4])[*+0:3].[O+0:5]=[O+0]>>[*:1][*:2](=[*:4])(=[*:5])[*:3]",
        (1, 0.5),
        (1,),
    ),
    # Thiones ###############################################
    # Thione-Thiol Tautomerization
    OperatorSmarts(
        "Thione-Thiol Tautomerization",
        "[N,O;+0!H0:1][C+0:2]=[SX1+0:3]>>[*:1]=[*:2][*:3]",
        (1,),
        (1,),
        kekulize_flag=True,
    ),
    OperatorSmarts(
        "Thione-Thiol Tautomerization Reverse",
        "[N,O;+0:1]=[C+0:2][SX2+0H:3]>>[*:1][*:2]=[*:3]",
        (1,),
        (1,),
        kekulize_flag=True,
    ),
    # Sulfenic, Sulfinic, Sulfonic, Sulfuric Acids ###############################################
    # Sulfenic Acid Tautomerization
    OperatorSmarts(
        "Sulfenic Acid Tautomerization",
        "[C,c;+0:1][SX2+0:2][O+0H:3]>>[*:1][*:2]=[*:3]",
        (1,),
        (1,),
    ),
    # Sulfenic Acid Tautomerization
    OperatorSmarts(
        "Sulfenic Acid Tautomerization Reverse",
        "[C,c;+0:1][SX3+0H:2]=[O+0:3]>>[*:1][*:2][*:3]",
        (1,),
        (1,),
    ),
    # Thiol Oxidation to Sulfenic Acid, Sulfinic Acid, Sulfonic Acid
    OperatorSmarts(
        "Thiol Oxidation to Sulfenic Acid",
        "[C,c;+0:1][SX2+0H:2].[O+0:3]=[O+0]>>[*:1][*:2][*:3]",
        (1, 0.5),
        (1,),
    ),
    OperatorSmarts(
        "Thiol Oxidation to Sulfinic Acid",
        "[C,c;!$(*=O)+0:1][SX2+0H:2].[O+0:3]=[O+0:4]>>[*:1][*:2](=[*:3])[*:4]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Thiol Oxidation to Sulfonic Acid",
        "[C,c;!$(*=O)+0:1][SX2+0H:2].[O+0:3]=[O+0:4]>>[*:1][*:2](=[*:3])(=[O])[*:4]",
        (1, 1.5),
        (1,),
    ),
    # Sulfurous Acid Salts Addition to Oxiranes
    OperatorSmarts(
        "Sulfurous Acid Salts Addition to Oxiranes",
        "[CX4+0:1]1[O+0:2][CX4+0:3]1.[O+0:4]=[SX3+0:5]([O+0H:6])[O+0H:7]>>[*:2][*:1][*:3][*:5](=[*:4])(=[*:6])[*:7]",
        (1, 1),
        (1,),
    ),
    # Esterification, Acid Anhydride Formation with Sulfuric Acid, Organosulfates Hydrolysis
    OperatorSmarts(
        "Esterification, Acid Anhydride Formation with Sulfuric Acid",
        "[C,c;+0:1][O+0H:2].[O+0:3]=[SX4+0:4](=[O+0:5])([OX2+0:6])[O+0H:7]>>[*:1][*:2][*:4](=[*:3])(=[*:5])[*:6].[*:7]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Organosulfates Hydrolysis",
        "[C,c;+0:1][O+0:2]-!@[SX4+0:4](=[O+0:3])(=[O+0:5])[OX2+0:6].[O+0H2:7]>>[*:1][*:2].[*:3]=[*:4](=[*:5])([*:6])[*:7]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Esterification, Acid Anhydride Formation with Sulfuric Acid, Intramolecular",
        "([C,c;+0:1][O+0H:2].[O+0:3]=[SX4+0:4](=[O+0:5])([OX2+0:6])[O+0H:7])>>[*:1][*:2][*:4](=[*:3])(=[*:5])[*:6].[*:7]",
        (1,),
        (1, 1),
    ),
    OperatorSmarts(
        "Organosulfates Hydrolysis, Intramolecular",
        "[C,c;+0:1][O+0:2]-@[SX4+0:4](=[O+0:3])(=[O+0:5])[OX2+0:6].[O+0H2:7]>>([*:1][*:2].[*:3]=[*:4](=[*:5])([*:6])[*:7])",
        (1, 1),
        (1,),
    ),
    # Organosulfates Reduction
    OperatorSmarts(
        "Organosulfates Reduction",
        "[C!$(*=O)+0:1][O+0:2][SX4+0:3](=[O+0:4])(=[O+0:5])[O+0H:6].[H][H]>>[*:1][*:3](=[*:4])(=[*:5])[*:6].[*:2]",
        (1, 1),
        (1, 1),
    ),
    # Alkenes Addition by Sulfuric Acid
    OperatorSmarts(
        "Alkenes Addition by Sulfuric Acid",
        "[C+0:1]=[C+0:2].[O+0:3]=[SX4+0:4](=[O+0:5])([O+0H:6])[O+0H:7]>>[*:1][*:2][*:6][*:4](=[*:3])(=[*:5])[*:7]",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Alkenes Addition by Sulfuric Acid Reverse",
        "[C+0!H0:1][C+0:2][O+0:6][SX4+0:4](=[O+0:3])(=[O+0:5])[O+0H:7]>>[*:1]=[*:2].[*:3]=[*:4](=[*:5])([*:6])[*:7]",
        (1,),
        (1, 1),
    ),
    # Alkenes Addition by Bisulfites
    OperatorSmarts(
        "Alkenes Addition by Bisulfites",
        "[C+0:1]=[C+0:2].[O+0:3]=[SX3+0:4]([O+0H:5])[O+0H:6]>>[*:1][*:2][*:4](=[*:3])(=[*:5])[*:6]",
        (1, 1),
        (1,),
    ),
    # Esterification of Sulfonic Acids, Sulfonic Acid Esters Hydrolysis
    OperatorSmarts(
        "Esterification of Sulfonic Acids",
        "[C,c;+0:1][O+0H:2].[O+0:3]=[SX4+0:4](=[O+0:5])([C,c,N;+0:6])[OH,F,Cl,Br,I;+0:7]>>[*:1][*:2][*:4](=[*:3])(=[*:5])[*:6].[*:7]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Sulfonic Acid Esters Hydrolysis",
        "[C,c;+0:1][O+0:2]-!@[SX4+0:4](=[O+0:3])(=[O+0:5])[C,c,N;+0:6].[O+0H2:7]>>[*:1][*:2].[*:3]=[*:4](=[*:5])([*:6])[*:7]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Esterification of Sulfonic Acids, Intramolecular",
        "([C,c;+0:1][O+0H:2].[O+0:3]=[SX4+0:4](=[O+0:5])([C,c,N;+0:6])[OH,F,Cl,Br,I;+0:7])>>[*:1][*:2][*:4](=[*:3])(=[*:5])[*:6].[*:7]",
        (1,),
        (1, 1),
    ),
    OperatorSmarts(
        "Sulfonic Acid Esters Hydrolysis, Intramolecular",
        "[C,c;+0:1][O+0:2]-@[SX4+0:4](=[O+0:3])(=[O+0:5])[C,c,N;+0:6].[O+0H2:7]>>([*:1][*:2].[*:3]=[*:4](=[*:5])([*:6])[*:7])",
        (1, 1),
        (1,),
    ),
    OperatorSmarts(
        "Sulfonic Esters Nucleophilic Substitution",
        "[C,c;+0:1][O+0:2]-!@[SX4+0:4](=[O+0:3])(=[O+0:5])[C,c;+0:6].[F,Cl,Br,I;+0H:7]>>[*:1][*:7].[*:3]=[*:4](=[*:5])([*:6])[*:2]",
        (1, 1),
        (1, 1),
    ),
    OperatorSmarts(
        "Sulfonic Esters Nucleophilic Substitution, Intramolecular",
        "[C,c;+0:1][O+0:2]-@[SX4+0:4](=[O+0:3])(=[O+0:5])[C,c;+0:6].[F,Cl,Br,I;+0H:7]>>([*:1][*:7].[*:3]=[*:4](=[*:5])([*:6])[*:2])",
        (1, 1),
        (1,),
    ),
    # Sulfonation of Benzene
    OperatorSmarts(
        "Sulfonation of Benzene",
        "[c+0H:1].[O+0:2]=[SX3+0:3](=[O+0:4])=[O+0:5]>>[*:1][*:3](=[*:4])(=[*:5])[*:2]",
        (1, 1),
        (1,),
    ),
    # Disproportionation of Aromatic Sulfinic Acids
    OperatorSmarts(
        "Disproportionation of Aromatic Sulfinic Acids",
        "[c+0:1][SX3+0:2](=[O+0])[O+0H]>>[*:1][*:2](=[O])(=[O])[O].[*:1][*:2](=[O])(=[O])[*:2][*:1].[O]",
        (3,),
        (1, 1, 1),
    ),
    # Alkylation of Sulfinic Acids with Halides
    OperatorSmarts(
        "Alkylation of Sulfinic Acids with Halides",
        "[C,c;+0:1][SX3+0:2](=[O+0:3])[O+0H:4].[C,c;+0:5][F,Cl,Br,I;+0:6]>>[*:1][*:2](=[*:3])(=[*:4])[*:5].[*:6]",
        (1, 1),
        (1, 1),
    ),
    # Reduction of Aromatic Sulfonyl Chlorides to Thiols
    OperatorSmarts(
        "Reduction of Aromatic Sulfonyl Chlorides to Thiols",
        "[c+0:1][SX4+0:2](=[O+0:3])(=[O+0])[F,Cl,Br,I;+0:4].[H][H]>>[*:1][*:2].[*:4].[*:3]",
        (1, 3),
        (1, 1, 2),
    ),
    # Reduction of Sulfonyl Chlorides to Sulfinic Acids
    OperatorSmarts(
        "Reduction of Sulfonyl Chlorides to Sulfinic Acids",
        "[C,c;+0:1][SX4+0:2](=[O+0:3])(=[O+0:4])[F,Cl,Br,I;+0:5].[H][H]>>[*:1][*:2](=[*:3])[*:4].[*:5]",
        (1, 1),
        (1, 1),
    ),
)
