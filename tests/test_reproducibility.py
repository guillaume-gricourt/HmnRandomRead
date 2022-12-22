import filecmp
import os
import shutil
import tempfile

from tests import main_tests


class TestReproducibility(main_tests.MainTests):
    def test_run(self):
        with tempfile.TemporaryDirectory() as outdir:
            shutil.copyfile(self.toy_fa, os.path.join(outdir, "toy.fa"))

            args = [self.exec_path]
            args += ["-r", os.path.join(outdir, "toy.fa") + ",100"]
            args += ["-r1", os.path.join(outdir, "one.r1.fastq.gz")]
            args += ["-r2", os.path.join(outdir, "one.r2.fastq.gz")]
            args += ["-s", "0"]

            main_tests.MainTests.cmd(args=args)

            # Compare files
            self.assertTrue(
                filecmp.cmp(os.path.join(outdir, "toy.fa.fai"), self.toy_fai)
            )
            self.assertTrue(
                filecmp.cmp(os.path.join(outdir, "toy.fa.scaff"), self.toy_scaff)
            )

            self.assertTrue(
                main_tests.MainTests.compare_fastq(
                    os.path.join(outdir, "one.r1.fastq.gz"), self.toy_one_r1
                )
            )
            self.assertTrue(
                main_tests.MainTests.compare_fastq(
                    os.path.join(outdir, "one.r2.fastq.gz"), self.toy_one_r2
                )
            )
