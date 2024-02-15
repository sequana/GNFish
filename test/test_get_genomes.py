import subprocess
import shutil
from pathlib import Path


def test_main(tmp_path):
    # dg.main()
    # data_dir = tmp_path / "data"
    # shutil.copytree(Path(__file__).parent / "data", data_dir)
    assert subprocess.call(
        "get_genomes hlorente ../gnfish/query.txt --genomic",
        shell=True)
