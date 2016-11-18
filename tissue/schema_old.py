# -*- coding: utf-8 -*-

from sqlalchemy import (Boolean, Column, Enum, ForeignKey, Integer, String,
                        Table)


_organs = [
    "adrenal",
    "anus",
    "aorta",
    "appendix",
    "areolar connective tissue",
    "bladder",
    "blood",
    "blood vessel",
    "bone",
    "brain",
    "breast",
    "cartilage",
    "cervix",
    "colon",
    "ear",
    "embryo",
    "esophagus",
    "eye",
    "eyelid",
    "fallopian",
    "fat",
    "fibrocartilage",
    "fingertip",
    "gallbladder",
    "heart",
    "intestine",
    "kidney",
    "larynx",
    "lip",
    "liver",
    "lung",
    "lymph node",
    "mesentery",
    "mouth",
    "muscle",
    "nail",
    "nasal septum",
    "nerve",
    "nose",
    "onion root tip",
    "ovary",
    "pancreas",
    "parathyroid",
    "parotid",
    "penis",
    "pineal",
    "pituitary",
    "placenta",
    "prostate",
    "salivary gland",
    "skin",
    "sperm",
    "spine",
    "spleen",
    "stomach",
    "tendon",
    "testis",
    "thymus",
    "thyroid",
    "tongue",
    "tooth",
    "trachea",
    "umbilical cord",
    "ureter",
    "urethra",
    "uterus",
    "vagina",
    "vertebra",
]

_organ_systems = [
    "blood and lymphatic vessels",
    "breast",
    "endocrine system",
    "eye",
    "female genital tract",
    "gastrointestinal tract",
    "head and neck",
    "heart",
    "hematopoietic system",
    "integument",
    "liver and biliary tract",
    "lung and respiratory tract",
    "male genital tract",
    "nervous system",
    "pancreas",
    "skeletal system",
    "supporting tissue and muscle",
    "undetermined",
    "urinary tract",
]

_species = [
    "canine",
    "feline",
    "guinea_pig",
    "human",
    "monkey",
    "mouse",
    "porcine",
    "rabbit",
    "rat",
    "undetermined",
]


class StrictEnum(Enum):
    """Enum that validates its input."""
    def __init__(self, *args, **kwargs):
        kwargs["validate_strings"] = True
        super().__init__(*args, **kwargs)


def create_slides_table(metadata):
    __table_args__ = {'extend_existing': True}
    slides = Table(
        "slides",
        metadata,
        Column("slide_id", Integer(), primary_key=True),
        Column("slide_name", String(30)),
        Column("organ", StrictEnum(*_organs), index=True),
        Column("organ_system", StrictEnum(*_organ_systems)),
        Column("species", StrictEnum(*_species)),
        Column("contributor", String(100)),
        Column("magnification", Integer()),
        sqlite_autoincrement=True,
        extend_existing=True
    )

    return slides


def create_labels_table(metadata):
    __table_args__ = {'extend_existing': True}
    labels = Table(
        "labels",
        metadata,
        Column("label_id", Integer(), primary_key=True),
        Column("organ", StrictEnum(*_organs), index=True),
        sqlite_autoincrement=True,
        extend_existing=True
    )

    return labels


def create_ingestion_tables(metadata):
    __table_args__ = {'extend_existing': True}
    ingestion = Table(
        "ingestion",
        metadata,
        Column("slide_id", ForeignKey("slides.slide_id"), primary_key=True),
        Column("ingested", Boolean(), default=False),
        keep_existing=True
    )

    return ingestion


def create_patches_table(metadata):
    __table_args__ = {'extend_existing': True}
    patches = Table(
        "patches",
        metadata,
        Column("patch_id", Integer(), primary_key=True),
        Column("slide_id", ForeignKey("slides.slide_id"), nullable=False),
        Column("patch_name", String(12), nullable=False),
        Column("path", String(100), nullable=False),
        Column("height", Integer(), nullable=False),
        Column("width", Integer(), nullable=False),
        Column("channels", Integer(), nullable=False),
        sqlite_autoincrement=True,
        extend_existing=True
    )

    return patches
