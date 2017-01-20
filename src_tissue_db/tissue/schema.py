# -*- coding: utf-8 -*-

from sqlalchemy import (Boolean, Column, Enum, ForeignKey, Integer, String,
                        Table, UniqueConstraint, Float)


_overall_outcomes = [
    "AliveNED",
    "AliveNEDlost",
    "AliveREC",
    "AliveREClost",
    "DeadNED",
    "DeadREC",
    "DeadLostNED",
    "DeadNEDRx",
]

#class StrictEnum(Enum):
#    """Enum that validates its input."""
#    def __init__(self, *args, **kwargs):
#        kwargs["validate_strings"] = True
#        super().__init__(*args, **kwargs)


#def create_slides_table(metadata):
#    __table_args__ = {'extend_existing': True}
#    slides = Table(
#        "slides",
#        metadata,
#        Column("slide_id", Integer(), primary_key=True),
#        Column("slide_name", String(30)),
#        sqlite_autoincrement=True,
#        extend_existing=True
#    )

#    return slides

def create_spots_table(metadata):
   __table_args__ = {'extend_existing': True}
   spots = Table(
        "spots", metadata,
        Column("spot_id",Integer(), primary_key=True),
        Column("spot_name", String(30), unique = True),
        Column("slide_name", String(30)),
        Column("position", Integer()),
        Column("sex",Boolean(),nullable = True),
        Column("age",Integer(),nullable = True),
        Column("chemo",Boolean(), default =False),
        Column("stage_sumajc",Integer(),nullable=False,index=True),
        Column("grade",Integer(),nullable=False,index=True),        
        Column("overall_outcome", Enum(*_overall_outcomes),nullable=True,index =True),
        Column("recurrent_5yr",Boolean(), default=False),
        sqlite_autoincrement=True, 
        keep_existing=True
   )

   return spots


def create_outcomes_table(metadata):
    __table_args__ = {'extend_existing': True}
    outcomes = Table(
        "outcomes",
        metadata,
        Column("outcome_id", Integer(), primary_key=True),
        Column("outcome", Integer(), index=True),
        sqlite_autoincrement=True,
        keep_existing=True
    )

    return outcomes


def create_ingestion_tables(metadata):
    __table_args__ = {'extend_existing': True}
    ingestion = Table(
        "ingestion",
        metadata,
        Column("spot_id", ForeignKey("spots.spot_id"), primary_key=True),
        Column("ingested", Boolean(), default=False),
        keep_existing=True
    )

    return ingestion


def create_cells_table(metadata):
    __table_args__ = {'extend_existing': True}
    cells = Table(
        "cells",
        metadata,
        Column("cell_id", Integer(), primary_key=True),
        Column("spot_id", ForeignKey("spots.spot_id"), nullable=False),
        Column("ID", Integer()),
        Column("x", Float(),nullable = False),
        Column("y", Float(),nullable = False),
        Column("Epithelial",Boolean(), nullable = False),
        sqlite_autoincrement=True,
        #extend_existing=True
        keep_existing=True
    )

    return cells

#def create_cells_table(metadata):
#    __table_args__ = {'extend_existing': True}
#    cells = Table('cell',metadata,
#          Column('cell_id',Integer(), primary_key = True),
#          Column("spot_id", ForeignKey("spots.spot_id"), nullable=False),
#          Column('slide_name', String(30)),
#    'slide_name', 'POSITION', 'ID', 'CellCentroid.X', 'CellCentroid.Y', 'Epithelium.Stroma', 'Perimeter', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Number.of.Nuclei', 'Area.Cell', 'Area.Nuclei', 'Area.Memb', 'Area.Cyt', 'SS', 'Median.Cell.pERK', 'Median.Cell.CD31', 'Median.Cell.BetaCatenin', 'Median.Cell.S6', 'Median.Cell.pS6235', 'Median.Cell.Fibronectin', 'Median.Cell.beta_actin', 'Median.Cell.pck26', 'Median.Cell.Glut1', 'Median.Cell.NaKATPase', 'Median.Cell.SMA', 'Median.Cell.Albumin', 'Median.Cell.CK19', 'Median.Cell.EGFR', 'Median.Cell.p4EBP1', 'Median.Cell.Wnt5', 'Median.Cell.pNDRG1', 'Median.Cell.FOXO3a', 'Median.Cell.MLH1', 'Median.Cell.E_cad', 'Median.Cell.pGSK3a', 'Median.Cell.LaminA_C', 'Median.Cell.pGSK3beta', 'Median.Cell.EZH2', 'Median.Cell.ALDH1', 'Median.Cell.p21', 'Median.Cell.Claudin1', 'Median.Cell.CD20', 'Median.Cell.IHH', 'Median.Cell.pEGFR', 'Median.Cell.NDRG1', 'Median.Cell.CD68', 'Median.Cell.TKLP1', 'Median.Cell.CD8', 'Median.Cell.CD79', 'Median.Cell.PTEN', 'Median.Cell.pMET_Y1349', 'Median.Cell.pMAPKAPK2', 'Median.Cell.FOXO1', 'Median.Cell.Akt', 'Median.Cell.CA9', 'Median.Cell.CleavedCaspase3', 'Median.Cell.ERK', 'Median.Cell.p_p38MAPK', 'Median.Cell.EPCAM', 'Median.Cell.CD3', 'Median.Cell.MSH2', 'Median.Cell.CD163', 'Median.Cell.cMET_epitomics', 'Median.Cell.PI3Kp110a', 'Median.Cell.4EBP1', 'Median.Cell.COX2', 'Median.Cell.PCNA', 'Median.Cell.p53', 'Median.Cell.DAPI', 'Median.Cell.ColIV', 'Median.Cyt.pERK', 'Median.Cyt.CD31', 'Median.Cyt.BetaCatenin', 'Median.Cyt.S6', 'Median.Cyt.pS6235', 'Median.Cyt.Fibronectin', 'Median.Cyt.beta_actin', 'Median.Cyt.pck26', 'Median.Cyt.Glut1', 'Median.Cyt.NaKATPase', 'Median.Cyt.SMA', 'Median.Cyt.Albumin', 'Median.Cyt.CK19', 'Median.Cyt.EGFR', 'Median.Cyt.p4EBP1', 'Median.Cyt.Wnt5', 'Median.Cyt.pNDRG1', 'Median.Cyt.FOXO3a', 'Median.Cyt.MLH1', 'Median.Cyt.E_cad', 'Median.Cyt.pGSK3a', 'Median.Cyt.LaminA_C', 'Median.Cyt.pGSK3beta', 'Median.Cyt.EZH2', 'Median.Cyt.ALDH1', 'Median.Cyt.p21', 'Median.Cyt.Claudin1', 'Median.Cyt.CD20', 'Median.Cyt.IHH', 'Median.Cyt.pEGFR', 'Median.Cyt.NDRG1', 'Median.Cyt.CD68', 'Median.Cyt.TKLP1', 'Median.Cyt.CD8', 'Median.Cyt.CD79', 'Median.Cyt.PTEN', 'Median.Cyt.pMET_Y1349', 'Median.Cyt.pMAPKAPK2', 'Median.Cyt.FOXO1', 'Median.Cyt.Akt', 'Median.Cyt.CA9', 'Median.Cyt.CleavedCaspase3', 'Median.Cyt.ERK', 'Median.Cyt.p_p38MAPK', 'Median.Cyt.EPCAM', 'Median.Cyt.CD3', 'Median.Cyt.MSH2', 'Median.Cyt.CD163', 'Median.Cyt.cMET_epitomics', 'Median.Cyt.PI3Kp110a', 'Median.Cyt.4EBP1', 'Median.Cyt.COX2', 'Median.Cyt.PCNA', 'Median.Cyt.p53', 'Median.Cyt.DAPI', 'Median.Cyt.ColIV', 'Median.Memb.pERK', 'Median.Memb.CD31', 'Median.Memb.BetaCatenin', 'Median.Memb.S6', 'Median.Memb.pS6235', 'Median.Memb.Fibronectin', 'Median.Memb.beta_actin', 'Median.Memb.pck26', 'Median.Memb.Glut1', 'Median.Memb.NaKATPase', 'Median.Memb.SMA', 'Median.Memb.Albumin', 'Median.Memb.CK19', 'Median.Memb.EGFR', 'Median.Memb.p4EBP1', 'Median.Memb.Wnt5', 'Median.Memb.pNDRG1', 'Median.Memb.FOXO3a', 'Median.Memb.MLH1', 'Median.Memb.E_cad', 'Median.Memb.pGSK3a', 'Median.Memb.LaminA_C', 'Median.Memb.pGSK3beta', 'Median.Memb.EZH2', 'Median.Memb.ALDH1', 'Median.Memb.p21', 'Median.Memb.Claudin1', 'Median.Memb.CD20', 'Median.Memb.IHH', 'Median.Memb.pEGFR', 'Median.Memb.NDRG1', 'Median.Memb.CD68', 'Median.Memb.TKLP1', 'Median.Memb.CD8', 'Median.Memb.CD79', 'Median.Memb.PTEN', 'Median.Memb.pMET_Y1349', 'Median.Memb.pMAPKAPK2', 'Median.Memb.FOXO1', 'Median.Memb.Akt', 'Median.Memb.CA9', 'Median.Memb.CleavedCaspase3', 'Median.Memb.ERK', 'Median.Memb.p_p38MAPK', 'Median.Memb.EPCAM', 'Median.Memb.CD3', 'Median.Memb.MSH2', 'Median.Memb.CD163', 'Median.Memb.cMET_epitomics', 'Median.Memb.PI3Kp110a', 'Median.Memb.4EBP1', 'Median.Memb.COX2', 'Median.Memb.PCNA', 'Median.Memb.p53', 'Median.Memb.DAPI', 'Median.Memb.ColIV', 'Median.Nuc.pERK', 'Median.Nuc.CD31', 'Median.Nuc.BetaCatenin', 'Median.Nuc.S6', 'Median.Nuc.pS6235', 'Median.Nuc.Fibronectin', 'Median.Nuc.beta_actin', 'Median.Nuc.pck26', 'Median.Nuc.Glut1', 'Median.Nuc.NaKATPase', 'Median.Nuc.SMA', 'Median.Nuc.Albumin', 'Median.Nuc.CK19', 'Median.Nuc.EGFR', 'Median.Nuc.p4EBP1', 'Median.Nuc.Wnt5', 'Median.Nuc.pNDRG1', 'Median.Nuc.FOXO3a', 'Median.Nuc.MLH1', 'Median.Nuc.E_cad', 'Median.Nuc.pGSK3a', 'Median.Nuc.LaminA_C', 'Median.Nuc.pGSK3beta', 'Median.Nuc.EZH2', 'Median.Nuc.ALDH1', 'Median.Nuc.p21', 'Median.Nuc.Claudin1', 'Median.Nuc.CD20', 'Median.Nuc.IHH', 'Median.Nuc.pEGFR', 'Median.Nuc.NDRG1', 'Median.Nuc.CD68', 'Median.Nuc.TKLP1', 'Median.Nuc.CD8', 'Median.Nuc.CD79', 'Median.Nuc.PTEN', 'Median.Nuc.pMET_Y1349', 'Median.Nuc.pMAPKAPK2', 'Median.Nuc.FOXO1', 'Median.Nuc.Akt', 'Median.Nuc.CA9', 'Median.Nuc.CleavedCaspase3', 'Median.Nuc.ERK', 'Median.Nuc.p_p38MAPK', 'Median.Nuc.EPCAM', 'Median.Nuc.CD3', 'Median.Nuc.MSH2', 'Median.Nuc.CD163', 'Median.Nuc.cMET_epitomics', 'Median.Nuc.PI3Kp110a', 'Median.Nuc.4EBP1', 'Median.Nuc.COX2', 'Median.Nuc.PCNA', 'Median.Nuc.p53', 'Median.Nuc.DAPI', 'Median.Nuc.ColIV'



#Column("Cell_biomarker", ARRAY(Float)), # how do I make it a vector?
#Column("Cyt_biomarker", ARRAY(Float)),
#Column("Mem_biomarker",ARRAY(Float)),
        
#def create_stages_table(metadata):
#    __table_args__ = {'extend_existing': True}
#    labels = Table(
#        "stages",
#        metadata,
#        Column("stage_id", Integer(), primary_key=True),
#        Column("stage", Integer(), index=True),
#        sqlite_autoincrement=True,
#        extend_existing=True
#    )
#
#    return stages
#
#def create_grades_table(metadata):
#    __table_args__ = {'extend_existing': True}
#    labels = Table(
#        "grades",
#        metadata,
#        Column("grade_id", Integer(), primary_key=True),
#        Column("grade", Integer(), index=True),
#        sqlite_autoincrement=True,
#        extend_existing=True
#    )
#
#    return grades
#_organs = [
#    "adrenal",
#    "anus",
#    "aorta",
#    "appendix",
#    "areolar connective tissue",
#    "bladder",
#    "blood",
#    "blood vessel",
#    "bone",
#    "brain",
#    "breast",
#    "cartilage",
#    "cervix",
#    "colon",
#    "ear",
#    "embryo",
#    "esophagus",
#    "eye",
#    "eyelid",
#    "fallopian",
#    "fat",
#    "fibrocartilage",
#    "fingertip",
#    "gallbladder",
#    "heart",
#    "intestine",
#    "kidney",
#    "larynx",
#    "lip",
#    "liver",
#    "lung",
#    "lymph node",
#    "mesentery",
#    "mouth",
#    "muscle",
#    "nail",
#    "nasal septum",
#    "nerve",
#    "nose",
#    "onion root tip",
#    "ovary",
#    "pancreas",
#    "parathyroid",
#    "parotid",
#    "penis",
#    "pineal",
#    "pituitary",
#    "placenta",
#    "prostate",
#    "salivary gland",
#    "skin",
#    "sperm",
#    "spine",
#    "spleen",
#    "stomach",
#    "tendon",
#    "testis",
#    "thymus",
#    "thyroid",
#    "tongue",
#    "tooth",
#    "trachea",
#    "umbilical cord",
#    "ureter",
#    "urethra",
#    "uterus",
#    "vagina",
#    "vertebra",
#]
#
#_organ_systems = [
#    "blood and lymphatic vessels",
#    "breast",
#    "endocrine system",
#    "eye",
#    "female genital tract",
#    "gastrointestinal tract",
#    "head and neck",
#    "heart",
#    "hematopoietic system",
#    "integument",
#    "liver and biliary tract",
#    "lung and respiratory tract",
#    "male genital tract",
#    "nervous system",
#    "pancreas",
#    "skeletal system",
#    "supporting tissue and muscle",
#    "undetermined",
#    "urinary tract",
#]

#
#_species = [
#    "canine",
#    "feline",
#    "guinea_pig",
#    "human",
#    "monkey",
#    "mouse",
#    "porcine",
#    "rabbit",
#    "rat",
#    "undetermined",
#]
#

