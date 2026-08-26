import requests
from .info import AFDBEntry


def uniprot_to_entries(uniprot_id: str):
    # If input is full AFDB ID like AF-P0AFB5-F1-model-v4, extract UniProt ID
    if uniprot_id.startswith("AF-"):
        # Format: AF-{uniprot}-F1-model-v{...}
        parts = uniprot_id.split("-")
        if len(parts) >= 2:
            uniprot_id = parts[1]
    r = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}")
    entries = r.json()
    # Handle API errors - API returns dict with "error" instead of list
    if isinstance(entries, dict) and "error" in entries:
        raise ValueError(f"AlphaFold DB API error: {entries['error']}")
    if not isinstance(entries, list):
        raise ValueError(f"AlphaFold DB API returned unexpected format: {type(entries)}")
    for i in range(len(entries)):
        entries[i] = AFDBEntry(data=entries[i])
    return entries
