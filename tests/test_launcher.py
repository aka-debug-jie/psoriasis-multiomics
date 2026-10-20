"""Behavior tests for the bundled launcher's environment and R kernel spec."""
import importlib.util
import json
import hashlib
import os
from pathlib import Path
import shutil
import tempfile
import subprocess
import sys
import tarfile
import unittest
from unittest.mock import patch


SOURCE = Path(__file__).resolve().parents[1] / "runtime" / "launch.py"


class MakeEnvironmentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        # Keep non-ASCII characters and spaces in the simulated resource root.
        self.resource_root = Path(self.temp.name) / "课程 包"
        self.resource_root.mkdir()
        launcher = self.resource_root / "launch.py"
        shutil.copy2(SOURCE, launcher)
        spec = importlib.util.spec_from_file_location("student_bundle_launch", launcher)
        self.launch = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.launch)
        self.runtime = self.resource_root / "work" / "runtime"
        self.r_binary = self.runtime / "bin" / "R"
        self.r_binary.parent.mkdir(parents=True)
        self.r_binary.touch()
        (self.runtime / "course-library").mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def test_environment_paths_are_inside_resource_root(self):
        env, binary = self.launch.make_environment()
        root = self.resource_root.resolve()
        library = self.runtime / "course-library"
        work = root / "work"

        self.assertEqual(binary, self.r_binary)
        self.assertEqual(Path(env["TEACHING_LIBRARY"]), library)
        self.assertEqual((root / "course" / env["TEACHING_NOTEBOOK_OUTPUT"]).resolve(), work / "results")
        self.assertEqual(Path(env["JUPYTER_PATH"]), work / "jupyter")
        self.assertEqual(Path(env["JUPYTER_CONFIG_DIR"]), work / "jupyter-config")
        self.assertEqual(Path(env["JUPYTER_RUNTIME_DIR"]), work / "jupyter-runtime")
        self.assertTrue((work / "jupyter/kernels/ir-psoriasis/kernel.json").is_file())
        self.assertTrue((work / "tmp").is_dir())

    def test_user_r_python_and_conda_settings_are_removed_and_path_is_rebuilt(self):
        global_values = {
            "PYTHONPATH": "/user/python/site",
            "PYTHONHOME": "/user/python/home",
            "R_LIBS": "/user/R/site",
            "R_LIBS_USER": "/user/R/library",
            "R_LIBS_SITE": "/user/R/site-library",
            "R_HOME": "/user/R/home",
            "RENV_PROJECT": "/user/project",
            "CONDA_PREFIX": "/user/conda",
            "TEACHING_SMR": "/user/smr",
            "TEACHING_PLINK": "/user/plink",
            "JUPYTER_ALLOW_INSECURE_WRITES": "True",
            "PATH": "/user/bin:/user/tools",
        }
        with patch.dict(os.environ, global_values):
            env, _ = self.launch.make_environment()

        for name in set(global_values) - {"PATH", "JUPYTER_ALLOW_INSECURE_WRITES"}:
            self.assertNotIn(name, env)
        self.assertNotIn("/user/", env["PATH"])
        self.assertEqual(
            env["PATH"],
            os.pathsep.join((str(self.runtime / "bin"), "/usr/bin:/bin:/usr/sbin:/sbin")),
        )
        self.assertEqual(env["PYTHONNOUSERSITE"], "1")
        self.assertEqual(env["JUPYTER_ALLOW_INSECURE_WRITES"], "False")
        self.assertEqual(env["TEACHING_LIBRARY"], str(self.runtime / "course-library"))

    def test_kernel_json_has_r_arguments_and_isolates_library(self):
        env, binary = self.launch.make_environment()
        kernel_file = Path(env["JUPYTER_PATH"]) / "kernels/ir-psoriasis/kernel.json"
        spec = json.loads(kernel_file.read_text(encoding="utf-8"))
        argv = spec["argv"]

        self.assertEqual(argv[:3], [str(binary), "--vanilla", "--slave"])
        self.assertEqual(argv[3], "-e")
        self.assertIn(".libPaths(c(", argv[4])
        self.assertIn(json.dumps(str(self.runtime / "course-library")), argv[4])
        self.assertIn(".Library", argv[4])
        self.assertIn("IRkernel::main()", argv[4])
        self.assertEqual(argv[5:], ["--args", "{connection_file}"])
        self.assertEqual(spec["language"], "R")
        self.assertEqual(spec["display_name"], "R · 多组学课程")

    def test_r_shell_assignments_accept_space_and_unicode_paths(self):
        import subprocess
        self.r_binary.write_text("#!/bin/sh\nR_HOME_DIR=/old/lib/R\nR_SHARE_DIR=/old/lib/R/share\nprintf '%s\\n' \"$R_HOME_DIR\" \"$R_SHARE_DIR\"\n",encoding="utf-8")
        self.launch.make_environment()
        output=subprocess.check_output(["sh",str(self.r_binary)],text=True).splitlines()
        self.assertEqual(output,[str(self.runtime / "lib/R"),str(self.runtime / "lib/R/share")])

    def test_unicode_and_spaces_survive_environment_and_kernel_configuration(self):
        env, binary = self.launch.make_environment()
        root = str(self.resource_root.resolve())
        spec_path = Path(env["JUPYTER_PATH"]) / "kernels/ir-psoriasis/kernel.json"
        spec = json.loads(spec_path.read_text(encoding="utf-8"))

        self.assertIn("课程 包", root)
        self.assertTrue(env["TEACHING_LIBRARY"].startswith(root + os.sep))
        self.assertEqual((self.resource_root / "course" / env["TEACHING_NOTEBOOK_OUTPUT"]).resolve(), self.resource_root / "work/results")
        self.assertEqual(spec["argv"][0], str(binary))
        self.assertIn("课程 包", spec["argv"][0])
        expression = spec["argv"][4]
        encoded_library = expression.removeprefix(".libPaths(c(")
        library, _ = json.JSONDecoder().raw_decode(encoded_library)
        self.assertEqual(library, str(self.runtime / "course-library"))
        self.assertTrue(spec_path.is_file())

    def test_entry_first_unpack_ignores_external_python_home(self):
        bundle=self.resource_root.parent / '入口 空格'
        resources=bundle/'resources';seed=self.resource_root/'seed'
        (resources/'runtime').mkdir(parents=True);(seed/'bin').mkdir(parents=True)
        (seed/'bin/python').symlink_to(sys.executable)
        (seed/'bin/conda-unpack').write_text('from pathlib import Path\nPath(__file__).with_name("unpacked").touch()\n')
        with tarfile.open(resources/'runtime/seed.tar.gz','w:gz')as archive:
            archive.add(seed/'bin',arcname='bin')
        (resources/'runtime/platform.txt').write_text('test-platform\n')
        (resources/'runtime/language.txt').write_text('en\n')
        (resources/'runtime/min_free_kb.txt').write_text('1024\n')
        (resources/'runtime/sha256.txt').write_text(hashlib.sha256((resources/'runtime/seed.tar.gz').read_bytes()).hexdigest())
        (resources/'launch.py').write_text('print(1+1)\n')
        entry=bundle/'运行课程.command';shutil.copy2(SOURCE.parent/'start.command',entry)
        env=os.environ.copy();env.update(PYTHONHOME='/nonexistent/course-python',PYTHONPATH='/nonexistent/course-site')
        result=subprocess.check_output(['bash',str(entry)],env=env,text=True)
        self.assertEqual(result.splitlines()[-1],'2')
        self.assertTrue((resources/'work/runtime/bin/unpacked').exists())

    def test_parallel_first_launch_prepares_once_and_preserves_student_work(self):
        bundle=self.resource_root.parent/'parallel course'
        resources=bundle/'resources';seed=self.resource_root/'parallel-seed'
        (resources/'runtime').mkdir(parents=True);(seed/'bin').mkdir(parents=True)
        (seed/'bin/python').symlink_to(sys.executable)
        (seed/'bin/conda-unpack').write_text('from pathlib import Path\nimport time\ntime.sleep(.4)\nwith Path(__file__).with_name("preparations").open("a") as f:f.write("ready\\n")\n')
        with tarfile.open(resources/'runtime/seed.tar.gz','w:gz')as archive:archive.add(seed/'bin',arcname='bin')
        for name,value in [('platform.txt','test-platform'),('language.txt','en'),('min_free_kb.txt','1024'),('sha256.txt',hashlib.sha256((resources/'runtime/seed.tar.gz').read_bytes()).hexdigest())]:
            (resources/'runtime'/name).write_text(value+'\n')
        (resources/'launch.py').write_text('print("opened")\n')
        entry=bundle/'Run course.command';shutil.copy2(SOURCE.parent/'start.command',entry)
        saved=resources/'work/results/student.txt';saved.parent.mkdir(parents=True);saved.write_text('keep')
        processes=[subprocess.Popen(['bash',str(entry)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)for _ in range(2)]
        for process in processes:
            output,error=process.communicate(timeout=15)
            self.assertEqual(process.returncode,0,error)
            self.assertEqual(output.splitlines()[-1],'opened')
        self.assertEqual((resources/'work/runtime/bin/preparations').read_text(),'ready\n')
        self.assertEqual(saved.read_text(),'keep')


if __name__ == "__main__":
    unittest.main()
