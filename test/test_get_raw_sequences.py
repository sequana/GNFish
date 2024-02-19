import subprocess
import shutil
import os
from pathlib import Path


def test_main(tmp_path):
    dir_to_copy = "data"
    data_dir = tmp_path / dir_to_copy
    shutil.copytree(Path(__file__).parent / dir_to_copy, data_dir)
    result = subprocess.run(
        f"get_raw_sequences --directory {data_dir}",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    output_file = (
        str(data_dir)
        + "/Boleophthalmus_pectinirostris_14439011_rna_Boleophthalmus_pectinirostris_AQP11_XM_020938655_1_extraction_RAW.fas"
    )
    assert os.path.exists(output_file)
    assert os.path.getsize(output_file)
