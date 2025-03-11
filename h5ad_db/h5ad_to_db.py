from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import scanpy as sc
import structlog
import duckdb

logger = structlog.get_logger(__name__)

def setup_argparse():
    parser = argparse.ArgumentParser(description='Crawl or scrape urls into the qdrant database')
    
    parser.add_argument('--scrape', 
                        type=str,  # Can change to int, float, etc.
                        nargs='+',  # '+' means one or more arguments
                        help='Scrape any number of urls.')
    
    parser.add_argument('--crawl',
                        type=str,  # Can change to int, float, etc.
                        help='Crawl from the specified url. Takes only one url')
    parser.add_argument('--class_grep',
                        type=str,  # Can change to int, float, etc.curl -LsSf https://astral.sh/uv/install.sh | sh
                        help='Use to specify a specific keyword in webpages to enter links under')
    parser.add_argument('--llm_extraction',
                        action='store_true',
                        help="Use a large language model to extract information from web pages before storing in database. NOTE: I personally feel this doesn't work as well as just scraping the whole thing")
    
    return parser.parse_args()

def create_duckdb(initial_h5ad, db_name="cell_database", db_path="db/"):
    adata = sc.read_h5ad(initial_h5ad)
    MakeDb(adata, db_name, db_path)

# Not really convinced this will work depending on how in-memory is implemented
def update_duckdb(h5ad, db_name="cell_database", db_path="db/"):
    adata = sc.read_h5ad(h5ad)
    new_asql = AnnSQL(adata=adata)
    existing_asql = AnnSQL(db= db_path + db_name + ".asql")
    asql.query("INSERT INTO existing_asql * FROM new_asql")


