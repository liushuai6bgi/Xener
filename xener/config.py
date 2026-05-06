from pathlib import Path

_DATA_URL = {
    "blastdb.zip":"https://xenor.dcs.cloud/api/public/download?file=blastdb.zip",
    "default_run_conf.yaml":"https://xenor.dcs.cloud/api/public/download?file=default_run_conf.yaml"}
_XENER_DATA_DIR = Path.home() / ".xener" / "data"

def _ensure_data(DATA_KEY:str=None, force_update:bool=False):
    assert DATA_KEY is None or DATA_KEY in _DATA_URL.keys(), "Invalid data key"
    
    if DATA_KEY is None:
        if _XENER_DATA_DIR.exists() and any(_XENER_DATA_DIR.iterdir()):
            return
        _XENER_DATA_DIR.mkdir(parents=True, exist_ok=True)
        print(f"[xener] Dependence on missing data, currently downloading dependent data...")
        for name, url in _DATA_URL.items():
            data_path = _XENER_DATA_DIR / name
            __download_data(name, data_path, url)
    else:
        name, url = DATA_KEY, _DATA_URL[DATA_KEY]
        data_path = _XENER_DATA_DIR / name
        if data_path.exists() and not force_update: return
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
                        print(f"\r[xener] Download {name} progress: {percent}%", end="", flush=True)
        if total == 0:
            print(f"\r[xener] Download {name} progress: 100%", end="")
        print()
    except Exception as e:
        if data_path.exists():
            data_path.unlink()
        raise RuntimeError(f"[xener] Data download {name} failed: {e}") from e

    if name.endswith(".zip"):
        try:
            import zipfile

            with zipfile.ZipFile(data_path, "r") as z:
                z.extractall(_XENER_DATA_DIR)
            data_path.unlink()
            print(f"[xener] The data dependency has been configured in {_XENER_DATA_DIR}")
        except Exception as e:
            if data_path.exists():
                data_path.unlink()
            raise RuntimeError(f"[xener] Data decompression failed: {e}") from e