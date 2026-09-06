# Reference databases

Not committed: freely available, but redistribution is not ours to grant.

| File | Source | How |
|---|---|---|
| `chebi_lite_3_stars.sdf.gz` | ChEBI, 3-star (manually curated) | fetched automatically by `scripts/database_matches.py` |
| `hmdb_structures.sdf` | HMDB | download **"Structures"** from <https://www.hmdb.ca/downloads> and unzip here. HMDB refuses scripted download, so this step is manual |

KEGG is not included. Its bulk data is licensed and cannot be redistributed or
re-derived here; the earlier counts in `ProcessedData/DatabaseMatches/` came
from a previous release and are kept as deposited.
