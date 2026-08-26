import requests
from pathlib import Path
from dataclasses import dataclass


@dataclass
class AFDBEntry:

    data: dict = None

    @property
    def id(self):
        return self.data["entryId"]

    @property
    def file_types(self):
        return ["bcif", "cif", "pdb", "pdbImage", "plddtDoc", "paeDoc"]

    @property
    def file_Urls(self):
        result = {}
        for file_type in self.file_types:
            key = f"{file_type}Url"
            result[key] = self.data.get(key, None)
        return result

    def download(self, file_type: str, folder_dir=".", requested_version=None):
        # If requested specific version, construct URL manually
        if requested_version is not None:
            # AlphaFold filename convention: AF-{uniprot}-F1-model_v{version}.{ext}
            # entryId is already AF-{uniprot}-F1
            # Note: it's model_v{version} with underscore, not hyphen
            url = f"https://alphafold.ebi.ac.uk/files/{self.id}-model_v{requested_version}.{file_type}"
        else:
            url = self.file_Urls[f"{file_type}Url"]
            if url is None:
                raise ValueError(f"File type {file_type} not found for entry {self.id}")
        
        r = requests.get(url)
        if r.status_code != 200:
            raise ValueError(f"Failed to download {file_type} from {url} (status code: {r.status_code})")
        
        # Get filename from URL always (includes version)
        filename = url.split("/")[-1]
        
        with open(Path(folder_dir) / filename, "wb") as f:
            f.write(r.content)
