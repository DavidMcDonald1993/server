
MAX_VAL = 950
MIN_VAL = 700
STEP = -50

DEFAULT_THRESHOLD = 750
DEFAULT_TARGET_COVERAGE = 0.1

THRESHOLDS = range(MAX_VAL, MIN_VAL, STEP)

# MAX_RECORDS = 10000

MAX_HITS_FOR_IMAGE = 1000

VALID_ORGANISMS = (
    "Homo sapiens", 
    "Rattus norvegicus", 
    "Mus musculus",
)

def filter_columns(cols, cols_to_remove):
    return list(filter(lambda col: col not in cols_to_remove, cols))