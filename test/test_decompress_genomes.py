import subprocess
import gnfish.decompress_genomes as dg
import shutil
from pathlib import Path


def test_main(tmp_path):
    # dg.main()
    data_dir = tmp_path / "data"
    shutil.copytree(Path(__file__).parent / "data", data_dir)
    subprocess.call(f"decompress_genomes --directory {data_dir}", shell=True)

    with open(tmp_path / "data" / "genome.fa") as f:
        assert f.readlines() == [">1\n", "CAGTCGATGCTAGTCGTAGCTGATGTCGATGCTA\n"]
