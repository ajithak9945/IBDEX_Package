from importlib.resources import files
from pathlib import Path


def artifact_path(name: str) -> Path:
    return Path(files("ibdex_projector").joinpath("artifacts", name))
