import subprocess
import shutil
import os
from pathlib import Path
import pdb


def test_main(tmp_path):
    dir_to_copy = "data"
    data_dir = tmp_path / dir_to_copy
    shutil.copytree(Path(__file__).parent / dir_to_copy, data_dir)
    result = subprocess.run(
        f"get_combined_seqs 'Bole_AQP11_12' {data_dir}/ --directory {data_dir} --pattern 'RAW_3_translated.fas'",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    output_file = str(data_dir) + "/Bole_AQP11_12_all_combined.fas"
    assert os.path.exists(output_file)
    assert os.path.getsize(output_file)
