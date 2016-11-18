#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Main script."""
import argparse
import os
import sys

from sqlalchemy import (MetaData, create_engine)

from tissue.ingest import (create_outcomes, ingest_ingestion, ingest_cells,
                           ingest_spots)


DATA = "data"
DB = "data/spots.db"
SPOTS = "input/clinical_data_with_chemo.csv"
CELLS = "input/Filtered.QCed.NoTcorr.SegExclude.HighDAPI.log2.Median.Norm.Cells.csv"
BIOMARKERS = "data/biomarkers.csv"

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--ingest-spots", action="store_true")
    parser.add_argument("--ingest-cells", action="store_true")
    #parser.add_argument("--suffix", default="_nonoverlap_256_files")
    return parser.parse_args(argv)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = parse_args(argv)
    
    if not os.path.isdir(DATA):
       os.makedirs(DATA)
    engine = create_engine("sqlite:///{}".format(DB))
    metadata = MetaData()
    metadata.reflect(bind=engine)

    if args.ingest_spots:
        ingest_spots(engine, metadata, SPOTS)
        create_outcomes(engine, metadata)

    if args.ingest_cells:
        ingest_ingestion(engine, metadata)
        ingest_cells(engine, metadata, CELLS, BIOMARKERS)

if __name__ == "__main__":
    sys.exit(main())
