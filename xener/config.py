from pathlib import Path

from .utils.logger import get_logger

logger = get_logger('xener.config')

_DATA_URL = {
    "blastdb.zip":"https://xenor.dcs.cloud/api/public/download?file=blastdb.zip",
    "default_run_conf.yaml":"https://xenor.dcs.cloud/api/public/download?file=default_run_conf.yaml"}
_XENER_DATA_DIR = Path.home() / ".xener" / "data"

def _ensure_data(DATA_KEY:str=None, force_update:bool=False):
    assert DATA_KEY is None or DATA_KEY in _DATA_URL.keys(), "Invalid data key"

    _XENER_DATA_DIR.mkdir(parents=True, exist_ok=True)
    if DATA_KEY is None:
        if _XENER_DATA_DIR.exists() and any(_XENER_DATA_DIR.iterdir()):
            logger.debug('_ensure_data: data dir %s already populated, skip download.', _XENER_DATA_DIR)
            return
        logger.info('Dependence on missing data, currently downloading dependent data...')
        for name, url in _DATA_URL.items():
            data_path = _XENER_DATA_DIR / name
            __download_data(name, data_path, url)
    else:
        name, url = DATA_KEY, _DATA_URL[DATA_KEY]
        data_path = _XENER_DATA_DIR / name
        if data_path.exists() and not force_update:
            logger.debug('_ensure_data: %s already exists, skip download.', data_path)
            return
        __download_data(name, data_path, url)
def __download_data(name:str, data_path:Path, url:str):
    try:
        import requests

        resp = requests.get(url, stream=True, timeout=60)
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        downloaded = 0
        with open(data_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        percent = min(downloaded * 100 // total, 100)
                        logger.info('Download %s progress: %s%%', name, percent)
        if total == 0:
            logger.info('Download %s progress: 100%%', name)
        logger.info('Download %s complete: %s bytes to %s', name, downloaded, data_path)
    except Exception as e:
        if data_path.exists():
            data_path.unlink()
        logger.error('Data download %s failed: %s', name, e)
        raise RuntimeError(f"[xener] Data download {name} failed: {e}") from e

    if name.endswith(".zip"):
        try:
            import zipfile

            with zipfile.ZipFile(data_path, "r") as z:
                z.extractall(_XENER_DATA_DIR)
            data_path.unlink()
            logger.info('The data dependency has been configured in %s', _XENER_DATA_DIR)
        except Exception as e:
            if data_path.exists():
                data_path.unlink()
            logger.error('Data decompression failed for %s: %s', name, e)
            raise RuntimeError(f"[xener] Data decompression failed: {e}") from e