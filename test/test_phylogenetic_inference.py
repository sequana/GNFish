import subprocess
import shutil
import os
from pathlib import Path


def test_main(tmp_path):
    dir_to_copy = "data"
    data_dir = tmp_path / dir_to_copy
    shutil.copytree(Path(__file__).parent / dir_to_copy, data_dir)
    result = subprocess.run(
        f"phylogenetic_inference '' {data_dir}/AQP16_8_ali_extract.fasta",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    output_file = str(data_dir) + "/AQP16_8_ali_extract.fasta_TEST_UFBS_alrt.treefile"
    assert os.path.exists(output_file)
    assert os.path.getsize(output_file)
