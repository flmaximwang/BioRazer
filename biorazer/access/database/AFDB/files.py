"""AlphaFold DB 预测结构文件下载。"""
import requests
from pathlib import Path

from ....logger import initialize_logger
from .info import AFDBEntry
from .query import uniprot_to_entries

# AlphaFold DB 提供的文件类型 (与 AFDBEntry.file_types 一致)
SUPPORTED_FMT = ["bcif", "cif", "pdb", "pdbImage", "plddtDoc", "paeDoc"]


def fetch(
    uniprot_id: str,
    fmt: str = "cif",
    download_dir: str | Path = ".",
    overwrite=False,
    logger=None,
):
    """下载 AlphaFold DB 中该 UniProt 登录号对应的预测结构文件。

    Parameters
    ----------
    uniprot_id : str
        UniProt 登录号, 例如 "P0DP23"。也支持完整 AlphaFold ID 如 "AF-P0DP23-F1-model-v4"。
    fmt : str
        文件类型, 见 SUPPORTED_FMT。
    download_dir : str | Path
        下载目录。
    overwrite : bool
        已存在同名文件时是否覆盖。
    logger : logging.Logger
        日志器, 缺省自动初始化。

    Returns
    -------
    Path
        下载得到的文件路径。
    """
    if not logger:
        logger = initialize_logger(__name__)
    if fmt not in SUPPORTED_FMT:
        raise ValueError(
            f"Unsupported format: {fmt} (supported: {', '.join(SUPPORTED_FMT)})"
        )
    
    # Parse requested version from full AFDB ID
    requested_version = None
    if uniprot_id.startswith("AF-") and "-model-v" in uniprot_id:
        # Extract version: AF-P0AFB5-F1-model-v4 -> v4 -> 4
        version_part = uniprot_id.split("-model-v")[-1]
        if version_part.isdigit():
            requested_version = int(version_part)
    
    entries = uniprot_to_entries(uniprot_id)
    if not isinstance(entries, list) or not entries or not hasattr(entries[0], 'data'):
        raise ValueError(f"No AFDB entry found for {uniprot_id}")
    
    # Select entry based on requested version
    entry = entries[0]  # Default to first entry (latest version if multiple fragments)
    if requested_version is not None:
        # Check if requested version is available
        all_versions = entry.data.get("allVersions", [])
        if requested_version not in all_versions:
            raise ValueError(
                f"Requested version v{requested_version} not available for {uniprot_id}. "
                f"Available versions: {all_versions}"
            )
        # Check if the file exists - AlphaFold DB only keeps latest version available for download
        url = f"https://alphafold.ebi.ac.uk/files/{entry.id}-model_v{requested_version}.{fmt}"
        r = requests.head(url)
        if r.status_code != 200:
            raise ValueError(
                f"Version v{requested_version} is listed in API but no download available. "
                f"AlphaFold DB only keeps the latest version (v{entry.data['latestVersion']}) "
                f"publicly available."
            )
    else:
        entry.requested_version = None
    
    # The entry already has correct URLs from API (for latest version or requested version)
    # The filename from URL already includes version, so we don't need to construct it
    url_key = f"{fmt}Url"
    if requested_version is None:
        # Get filename from the API URL which already has version
        url = entry.data.get(url_key, "")
        if url:
            file_name = url.split("/")[-1]
        else:
            file_name = f"{entry.id}.{fmt}"
    else:
        # We'll construct URL in download(), filename will be correct there
        # entryId is AF-P0AFB5-F1, filename will be AF-P0AFB5-F1-model_vX.fmt
        file_name = f"{entry.id}-model_v{requested_version}.{fmt}"
    
    file_path = Path(download_dir) / file_name
    if file_path.exists() and not overwrite:
        logger.warning(f"{file_path} already exists, skipping")
        return file_path
    Path(download_dir).mkdir(parents=True, exist_ok=True)
    entry.download(fmt, folder_dir=str(download_dir), requested_version=requested_version)
    logger.info(f"{file_name} downloaded to {download_dir}")
    return file_path
