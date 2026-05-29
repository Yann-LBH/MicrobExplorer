from utils_qc_Reads import run_full_qc as run_reads
from utils_qc_Contigs import run_full_qc as run_contigs
#from utils_qc_Kegg import run_full_qc as run_kegg

source = getattr(snakemake.wildcards, "source", None)
if source == "reads":
    run_reads()
elif source == "contigs":
    run_contigs()
#elif source == "kegg":
#    run_kegg()
else:
    raise ValueError(f"Unsupported source wildcard: {source}")
