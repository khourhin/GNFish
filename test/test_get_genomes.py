import subprocess
import shutil
from pathlib import Path


def test_main(tmp_path):
    # dg.main()
    data_dir = tmp_path / "data"
    shutil.copytree(Path(__file__).parent / "data", data_dir)
    result = subprocess.run(
        f"get_genomes {data_dir} hlorente {data_dir}/query.txt ",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert result.returncode == 0, f"Command failed with output:\n"
    "{result.stdout}\n{result.stderr}"
