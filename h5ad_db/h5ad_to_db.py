import scanpy as sc
import structlog
from qdrant_client import QdrantClient
import argparse
import time
import psycopg
import asyncio

# NEW PLAN: Just go ahead and ingest the h5ad files into postres using psycopg and qdrant
# Ideally do this asynchronously, so we can eventually have multiple clients connect and add 
# rows to the DB
logger = structlog.get_logger(__name__)

# connection establishment
# obv will need to change on production
conn = psycopg2.connect(
    database="cell_database",
    user='postgres',
    password='password',
    host='localhost',
    port= '5432'
)
cursor = conn.cursor()

client = QdrantClient(host="localhost", port=6333)

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

def create_dbs(initial_h5ad, db_name="cell_database", db_path="db/"):
    # A problem we'll run into is that h5ad files have different observations. 
    # Ig this is not the biggest problem, users should just be aware that 
    # they need to make sure they are happy with their var and obs before creating the database
    adata = read_h5ad(initial_h5ad)
    obs_names = adata.obs_names
    # This will make up most of the columns of X
    var_names = adata.var_names
    # For now only 4 tables, X (cell_id as primary key + var_names as columns), var (var column names as columns + var_id as foreign key),
    # and a UMAP embedding table (cell_id as foreign key). A problem is that X is sparse, so we can get ~ 4-10x better space efficiency
    # by not representing null. One solution is to store 3 columns, the data (csr.data), the row pointers (csr.indptr) and the indicies (csr.indicies) 
    # This does mean that there's more overhead  and effort to represent gene expression dist from a column, but for our main use cases this is
    # acceptable, since we expect to mainly be working with med sized cell popualtions (between 10-100)
    # SOME IMPORTANT NOTES!!!
    # Because of how X is stored, it is important that data is preprocessed to fit the gene # and ordering. We are gonna use tabula sapiens gene # and
    # ordering for this. You can easily create a new table to do this, since yo ucan just change the max array size for indicies, but you will get
    # wonky results if you use the existing database without this info. i.e. please go by the var name foreign key for how new things must beformtated.
    # you should enforce this rule for update_dbs
    # Also prefer unscaled data, but it's actually not the end of the world if some rows are scaled provided we check that later
    sql = '''CREATE TABLE X {
            cell_id varchar(40) CONSTRAINT firstkey PRIMARY KEY,
            data double[],
            indptr int[],
            indicies int[],
            };'''
    sql = '''CREATE TABLE var {
            };'''
    sql = '''CREATE TABLE obs {
            };'''
    sql = '''CREATE TABLE umap {
            cell_id varchar(40) FOREIGN KEY,
            x double,
            y double,
            };'''
    cursor.execute(sql)
    conn.commit()

async def update_dbs(h5ads, db_name="cell_database", db_path="db/"):
    async for path in h5ads:
        adata = await sc.read_h5ad(path)
        # awa
    pass


if __name__ == "__main__":
    args = setup_argparse()
    if args.create:
        h5ad_file = args.create
        logger.info(f"h5ad path received: {h5ad_file}")
        create_dbs(h5ad_file)

