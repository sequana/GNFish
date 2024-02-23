import subprocess
import shutil
import os
from pathlib import Path


def test_main(tmp_path):
    dir_to_copy = "data"
    data_dir = tmp_path / dir_to_copy
    shutil.copytree(Path(__file__).parent / dir_to_copy, data_dir)
    subprocess.run(
        f"get_query_sequences hlorente {data_dir}/query_seqs_3.txt",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    output_file = "../Data/Query_seqs/aquaporin_query_seqs.fas"
    assert os.path.exists(output_file)
    assert os.path.getsize(output_file)
