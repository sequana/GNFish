import subprocess
import shutil
import os
from pathlib import Path


def test_main(tmp_path):
    dir_to_copy = "data"
    data_dir = tmp_path / dir_to_copy
    shutil.copytree(Path(__file__).parent / dir_to_copy, data_dir)
    result = subprocess.run(
        f"blast '' {data_dir}/protein_query_seqs.fas prot --directory {data_dir} --genomic",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    output_file = str(data_dir) + "/Bole_pseudo_genome_out.tsv"
    assert os.path.exists(output_file)
    assert os.path.getsize(output_file)
