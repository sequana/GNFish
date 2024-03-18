import subprocess
import shutil
import os
from pathlib import Path


def test_main(tmp_path):
    dir_to_copy = "data"
    data_dir = tmp_path / dir_to_copy
    shutil.copytree(Path(__file__).parent / dir_to_copy, data_dir)
    result = subprocess.run(
        f"alignment_trimming '' {data_dir}/Bole_AQP11_12_all_combined_3_ali.fas {data_dir}/Bole_AQP11_12_all_combined_3_ali_trim.fas",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    output_file = str(data_dir) + "/Bole_AQP11_12_all_combined_3_ali_trim.fas"
    assert os.path.exists(output_file)
    assert os.path.getsize(output_file)
