
MAX_VAL = 950
MIN_VAL = 750
STEP = -50

DEFAULT_THRESHOLD = 750

# MAX_RECORDS = 10000

MAX_HITS_FOR_IMAGE = 0

VALID_ORGANISMS = (
    "Homo sapiens", 
    "Rattus norvegicus", 
    "Mus musculus",
)

def filter_columns(cols, cols_to_remove):
    return list(filter(lambda col: col not in cols_to_remove, cols))