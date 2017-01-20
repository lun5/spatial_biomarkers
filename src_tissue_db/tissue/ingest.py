# -*- coding: utf-8 -*-
"""Ingestion script."""
from collections import namedtuple
import csv
import os

import numpy as np
from sqlalchemy import (not_, select, and_)
from sqlalchemy.exc import IntegrityError

from tissue.schema import (create_ingestion_tables, create_cells_table,
                           create_outcomes_table, create_spots_table)

fields = ['case_ID','slide_name','POSITION_ORDINAL','SEX','AGE.DX','chemo.final',
          'SUMAJC','GRADE','Overall.Outcome','recurevent5yrs']
not_NA_fields = ['GRADE', 'SUMAJC','Overall.Outcome','recurevent5yrs']
Spot = namedtuple("Spot", 'id name')


def create_outcomes(engine, metadata):
    outcomes = create_outcomes_table(metadata)
    outcomes.create(engine, checkfirst=True)

    spots = metadata.tables["spots"]
    s = select([spots.c.overall_outcome.distinct()])

    with engine.connect() as conn:
        result = conn.execute(s)
        outcome_list = [{"overall_outcome": r[0]} for r in result.fetchall()]

    ins = outcomes.insert()
    with engine.connect() as conn:
        result = conn.execute(ins, outcome_list)


def ingest_spots(engine, metadata, path):
    # create relation
    spots = create_spots_table(metadata)
    metadata.create_all(engine)
    columns = spots.columns.keys()
    columns.remove("spot_id")

    spot_list = []
    with open(path) as f:
        f_csv = csv.DictReader(f)
        for row in f_csv:
            has_NA = [('NA' in row[r]) for r in fields]
            if sum(has_NA) == 0:
                spot_list.append({'spot_name': row['case_ID'], 
                                  'slide_name': row['slide_name'],
                                  'position': int(row['POSITION_ORDINAL']),
                                  'sex': int(row['SEX'])==1,
                                  'age': int(row['AGE.DX']),
                                  'chemo': int(row['chemo.final']) ==1,
                                  'stage_sumajc': int(row['SUMAJC']),
                                  'grade': int(row['GRADE']),
                                  'overall_outcome': row['Overall.Outcome'],
                                  'recurrent_5yr': int(row['recurevent5yrs']) == 1})

    # insert tuples into relation
    ins = spots.insert()
    with engine.connect() as conn:
        result = conn.execute(ins, spot_list[:700])


def ingest_ingestion(engine, metadata):
    ingestion = create_ingestion_tables(metadata)
    ingestion.create(engine, checkfirst=True)
    spots = metadata.tables["spots"]

    s = select([spots.c.spot_id])
    ins = ingestion.insert()
    try:
        with engine.connect() as conn:
            result = conn.execute(s)
            spots = []
            for row in sorted(result):
                spots.append(
                    {"spot_id": row["spot_id"],
                    "ingested": False}
                )
            conn.execute(ins, spots)
    except IntegrityError:
        pass


def ingest_cells(engine, metadata, path, data_save_path):
    ingestion = metadata.tables["ingestion"]
    spots = metadata.tables['spots']
    not_ingested = select([ingestion.c.spot_id]).where(
        not_(ingestion.c.ingested)
    )
        
    # create relations
    cells = create_cells_table(metadata)
    cells.create(engine, checkfirst=True)
    
    cell_list = []
    data = []
    spots_NA = []
    with open(path) as f:
        f_csv = csv.DictReader(f)
        data_fields = f_csv.fieldnames[3:15] + f_csv.fieldnames[16:]
        #print data_fields
        for row in f_csv:
            has_NA = [('NA' in row[r]) for r in f_csv.fieldnames]
            if sum(has_NA) == 0:
                #slide_name = row['slide_name']
                #position = row['POSITION']
                spot_name = row['SS']
                sel = select([spots.c.spot_id, spots.c.spot_name]).where(
                    spots.c.spot_name == spot_name)
                #sel = sel.where(and_(spots.c.slide_name == slide_name,
                #                     spots.c.position == position))
                with engine.connect() as conn:
                    result = conn.execute(sel)
                    spot_list = [Spot(*r) for r in result.fetchall()]
                 
                if len(spot_list) == 0:
                    #print 'there is no such spot with spot_name ', spot_name
                    if not (spot_name in spots_NA):
                        spots_NA.append(spot_name)
                    continue
                elif len(spot_list) > 1:
                    print 'ERROR: there are ', str(len(spot_list)), ' with slide_name ', slide_name,' and postion ', position
                    return
                else:
                    #continue
                    spot = spot_list[0]
                    spot_id = spot.id
                    data.append([float(row[r]) for r in data_fields])
                
                    p = {
                        'spot_id': spot_id,
                        'ID': int(row['ID']),
                        'x': float(row['CellCentroid.X']),
                        'y': float(row['CellCentroid.Y']),
                        'Epithelial': int(row['Epithelium.Stroma']) == 1,
                    }
                    cell_list.append(p)
                    upd = ingestion.update().where(and_(ingestion.c.spot_id == spot_list[0].id,
                                                       ingestion.c.ingested == False))
                    upd = upd.values(ingested = True)
    
                    with engine.connect() as conn:
                        result = conn.execute(upd)
    print spots_NA
    with open('spot_NA.csv','wb') as myfile:
        wr = csv.writer(myfile, quoting = csv.QUOTE_ALL)
        wr.writerow(spots_NA)

    print 'Done with cell list'
    
    # save the data
    data = np.asarray(data)
    #np.savetxt(data_save_path, data, delimiter=',') 
    np.save(data_save_path,data)
    # insert tuples into relation
    ins = cells.insert()
    with engine.connect() as conn:
        result = conn.execute(ins, cell_list)
    
