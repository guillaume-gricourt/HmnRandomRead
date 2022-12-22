import gzip
import logging
import os
import subprocess
import unittest
from typing import List


class MainTests(unittest.TestCase):
    dataset_path = os.path.join(os.path.dirname(__file__), "dataset")
    exec_path = os.path.join(os.path.dirname(__file__), "HmnRandomRead")
    # Gold input
    toy_fa = os.path.join(dataset_path, "toy.fa")
    toy_fai = os.path.join(dataset_path, "toy.fa.fai")
    toy_scaff = os.path.join(dataset_path, "toy.fa.scaff")
    toy_one_r1 = os.path.join(dataset_path, "one.r1.fastq.gz")
    toy_one_r2 = os.path.join(dataset_path, "one.r2.fastq.gz")

    @classmethod
    def cmd(
        cls, args: List[str], show_output: bool = True
    ) -> subprocess.CompletedProcess:
        """Run a command line.

        Parameters
        ----------
        args: List[str]
            A list of argument
        show_output: bool (default: True)
            Output command line

        Return
        ------
        subprocess.CompletedProcess
            Return result obtained with subprocess
        """
        code = subprocess.run(args, capture_output=True, encoding="utf8")
        if show_output and code.stdout is not None:
            logging.info(code.stdout)
        if show_output and code.stderr is not None:
            logging.warning(code.stderr)
        return code

    @classmethod
    def load_fastq(cls, path: str) -> str:
        data = ""
        if path.endswith("gz") or path.endswith("gzip"):
            with gzip.open(path) as fid:
                data = fid.read().decode("utf-8")
        else:
            with open(path) as fid:
                data = fid.read()
        return data

    @classmethod
    def compare_fastq(cls, one: str, two: str) -> bool:
        first = cls.load_fastq(path=one)
        second = cls.load_fastq(path=two)
        return first == second
