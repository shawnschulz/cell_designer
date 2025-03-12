import scanpy as sc
import structlog
from qdrant_client import QdrantClient
from qdrant_client.models import Distance, VectorParams
import time
import argparse

logger = structlog.get_logger(__name__)
client = QdrantClient(host="localhost", port=6333)

def setup_argparse():
    parser = argparse.ArgumentParser(description='Generate or update vector database fromr h5ad files')
    
    parser.add_argument('--create', 
                        type=str,  # Can change to int, float, etc.
                        nargs='+',  # '+' means one or more arguments
                        help='Create the initial database')
    
    parser.add_argument('--update',
                        type=str,  # Can change to int, float, etc.
                        nargs = '+',
                        help='Update the datbase')
    
    return parser.parse_args()

# Idea here is to add hte gene sequences before adding to the vector database,
# ideally this should make clustering a breeze and add search functionality
# both for our receptor classificaiton task and for cell type classification
# among other parts of the problem space
def add_gene_sequences_to_h5ad(h5ad):
    pass

def create_vector_database(initial_h5ads, db_name="cell_database", vector_size=768):
    logger.info(f"Creating database {db_name} from the following paths: {initial_h5ads}")
    count = 0
    client.create_collection(collection_name=db_name, vectors_config=VectorParams(size=vector_size, distance=Distance.DOT))
    
    for path in initial_h5ads:
        logger.info(f"Path processing: {path}")
        try:
            adata = sc.read_h5ad(path)
        except Exception as e:
            logger.error(f"Reading path: {path} failed, trying next path")
            logger.error(f"{e}")
            continue
        start_time = time.perf_counter()
        ### INSERT CODE TO INSERT INTO VECTOR DATABASE HERE ###
        end_time = time.perf_counter()
        logger.info(f"Database creation completed, took {(start_time - end_time):.4f} seconds")

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

