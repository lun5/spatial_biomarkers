# -*- coding: utf-8 -*-
"""Ingestion script."""
from collections import namedtuple
import csv
import os

from scipy.misc import imread
from sqlalchemy import (not_, select)
from sqlalchemy.exc import IntegrityError

from tissue.schema import (create_ingestion_tables, create_patches_table,
                           create_labels_table, create_slides_table)


Slide = namedtuple("Slide", "id name label")


def create_labels(engine, metadata):
    labels = create_labels_table(metadata)
    labels.create(engine, checkfirst=True)

    slides = metadata.tables["slides"]
    s = select([slides.c.organ.distinct()])

    with engine.connect() as conn:
        result = conn.execute(s)
        organ_list = [{"organ": r[0]} for r in result.fetchall()]

    ins = labels.insert()
    with engine.connect() as conn:
        result = conn.execute(ins, organ_list)


def ingest_slides(engine, metadata, path):
    # create relation
    slides = create_slides_table(metadata)
    metadata.create_all(engine)
    columns = slides.columns.keys()
    columns.remove("slide_id")

    # collect tuples
    slide_list = []
    with open(path) as f:
        f_csv = csv.DictReader(f)
        for row in f_csv:
            slide_list.append(row)

    # insert tuples into relation
    ins = slides.insert()
    with engine.connect() as conn:
        result = conn.execute(ins, slide_list)


def ingest_ingestion(engine, metadata):
    ingestion = create_ingestion_tables(metadata)
    ingestion.create(engine, checkfirst=True)
    slides = metadata.tables["slides"]

    s = select([slides.c.slide_id])
    ins = ingestion.insert()
    try:
        with engine.connect() as conn:
            result = conn.execute(s)
            slides = []
            for row in sorted(result):
                slides.append(
                    {"slide_id": row["slide_id"],
                    "ingested": False}
                )
            conn.execute(ins, slides)
    except IntegrityError:
        pass


def get_patches(slide, path, suffix):
    dirname = slide.name.split(".")[0] + suffix
    start = os.path.join(path, dirname)
    for relpath, dirs, files in os.walk(start, followlinks=True):
        patches = (f for f in files if f.endswith(".jpeg") and f.find('rotation') == -1)
        for patch in patches:
            im = imread(os.path.join(relpath, patch))
            height, width, channels = im.shape
            p = {
                "slide_id": slide.id,
                "patch_name": patch,
                "path": os.path.relpath(relpath, start=path),
                "height": height,
                "width": width,
                "channels": channels
            }
            yield p


def ingest_patches(engine, metadata, path, suffix):
    # get names of slides that have not been processed yet
    slides = metadata.tables["slides"]
    ingestion = metadata.tables["ingestion"]
    labels = metadata.tables["labels"]
    not_ingested = select([ingestion.c.slide_id]).where(
        not_(ingestion.c.ingested)
    )
    sel = select([slides.c.slide_id, slides.c.slide_name, labels.c.label_id])
    sel = sel.select_from(
        slides.join(labels, slides.c.organ == labels.c.organ)).where(
        slides.c.slide_id.in_(not_ingested)
    )

    with engine.connect() as conn:
        result = conn.execute(sel)
        slide_list = [Slide(*r) for r in result.fetchall()]

    # create relation
    patches = create_patches_table(metadata)
    patches.create(engine, checkfirst=True)

    # batch insert all patches for a given slide
    batches = ((s.id, get_patches(s, path, suffix)) for s in slide_list)
    ins = patches.insert()
    for (slide, batch) in batches:
        patch_list = list(batch)
        if not patch_list:
            continue
        upd = ingestion.update().where(ingestion.c.slide_id == slide).values(
            ingested = True)
        with engine.connect() as conn:
            result = conn.execute(ins, patch_list)
            result = conn.execute(upd)
