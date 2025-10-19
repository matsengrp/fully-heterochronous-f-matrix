from os import system
from tempfile import TemporaryDirectory
import json
from python.f_matrix import sample_heterochronous_f_matrix
from collections import Counter
import unittest
from time import time
import numpy as np
import pytest
import shutil
import subprocess


def requires_r_with_ape():
    """Skip test if R with ape package is not available."""
    if not shutil.which("Rscript"):
        pytest.skip("R comparison tests require R with ape package")
    
    try:
        result = subprocess.run(['Rscript', '-e', 'library("ape")'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode != 0:
            pytest.skip("R comparison tests require R with ape package installed")
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pytest.skip("R comparison tests require R with ape package")


def sample_with_R(n, r):
    with TemporaryDirectory() as temp_dir:
        outpath = f"{temp_dir}/r.npy"
        system(f"Rscript Rscript/generate_fmatrices_optimized.R {n} {r} {outpath}")
        f_mats = np.load(outpath)
    f_mats = list(map(lambda m: tuple(map(tuple, m)), f_mats))
    return f_mats


def tensor_to_int_tuple(the_tensor):
    return tuple(map(lambda x: tuple(map(int, x)), the_tensor.tolist()))


def sample_with_python(n, r):
    f_mats = [
        tensor_to_int_tuple(sample_heterochronous_f_matrix(n)[0]) for _ in range(r)
    ]
    return f_mats


class TestFMatrixSampling(unittest.TestCase):

    @pytest.mark.r_comparison
    def test_samples(self):
        requires_r_with_ape()
        
        for n in range(3, 6):
            r = 10 ** (n + 2)
            t0 = -time()
            f_mats_from_r = Counter(sample_with_R(n, r))
            t0 += time()
            print(f"R sampling took {int(t0)} seconds")
            t0 = -time()
            f_mats_from_python = Counter(sample_with_python(n, r))
            t0 += time()
            print(f"python sampling took {int(t0)} seconds")

            missing_in_R = sum(m not in f_mats_from_r for m in f_mats_from_python)
            missing_in_python = sum(m not in f_mats_from_python for m in f_mats_from_r)
            self.assertEqual(0, missing_in_R)
            self.assertEqual(0, missing_in_python)

            max_diff = max(
                abs(python_count - f_mats_from_r[mat])
                for mat, python_count in f_mats_from_python.items()
                if mat in f_mats_from_r
            )
            max_diff /= r
            self.assertLessEqual(max_diff, 0.005)


if __name__ == "__main__":
    unittest.main()
