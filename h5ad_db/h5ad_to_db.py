from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import scanpy as sc
import structlog
import duckdb
import lancedb
import argparse
import time

logger = structlog.get_logger(__name__)

def setup_argparse():
    parser = argparse.ArgumentParser(description='Generate or update database from h5ad files')
    
    parser.add_argument('--create', 
                        type=str,  # Can change to int, float, etc.
                        nargs='+',  # '+' means one or more arguments
                        help='Create the initial database')
    
    parser.add_argument('--update',
                        type=str,  # Can change to int, float, etc.
                        nargs = '+',
                        help='Update the datbase')
    
    return parser.parse_args()

def create_duckdb(initial_h5ads, db_name="cell_database", db_path="db/"):
    logger.info(f"Creating database {db_name} from the following paths: {initial_h5ads}")
    count = 0
    
    for path in initial_h5ads:
        logger.info(f"Path processing: {path}")
        if count == 0:
            try:
                adata = sc.read_h5ad(path)
            except Exception as e:
                logger.error(f"Reading path: {path} failed, trying next path")
                logger.error(f"{e}")
                continue
            start_time = time.perf_counter()
            MakeDb(adata=adata, db_name=db_name, db_path=db_path, make_buffer_file=True)
            end_time = time.perf_counter()
            logger.info(f"Database creation completed, took {(start_time - end_time):.4f} seconds")
        else:
            try:
                adata = sc.read_h5ad(path)
            except Exception as e:
                logger.error(f"Reading path: {path} failed, trying next path")
                logger.error(f"{e}")
                continue
            start_time = time.perf_counter()
            new_asql = AnnSQL(adata=adata)
            existing_asql = AnnSQL(db= db_path + db_name + ".asql")
            existing_asql.query("INSERT INTO existing_asql * FROM new_asql")
            end_time = time.perf_counter()
            logger.info(f"Database update completed, took {(start_time - end_time):.4f} seconds")
        count += 1

# Not really convinced this will work depending on how in-memory is implemented
# Might need to create a temporary on disk datbase and try using duck db bindings
# to update the db
def update_duckdb(h5ads, db_name="cell_database", db_path="db/"):
    for path in h5ads:
        adata = sc.read_h5ad(path)
        new_asql = AnnSQL(adata=adata)
        existing_asql = AnnSQL(db= db_path + db_name + ".asql")
        existing_asql.query("INSERT INTO existing_asql * FROM new_asql")


if __name__ == "__main__":
    args = setup_argparse()
    if args.create:
        h5ad_files = args.create
        logger.info(f"h5ad path received: {h5ad_files}")
        create_duckdb(h5ad_files)

