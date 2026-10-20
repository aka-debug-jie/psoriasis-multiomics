"""Regression: multi-command modules must execute children, not report an empty command as success."""
from pathlib import Path
import shutil,subprocess,tempfile,unittest

ROOT=Path(__file__).resolve().parents[2]
RSCRIPT=ROOT/'research/envs/r442_sc/bin/Rscript'


class DispatchCommands(unittest.TestCase):
    def run_fixture(self,module,fail=None):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);research=root/'research';research.mkdir()
            shutil.copyfile(ROOT/'research/execute.R',research/'execute.R')
            stub=("#!/usr/bin/env python3\nimport pathlib,sys\n"+
                  f"root=pathlib.Path({str(root)!r})\n"+
                  "name=pathlib.Path(sys.argv[1]).name\n"+
                  "with (root/'calls.txt').open('a') as f:f.write(name+'\\n')\n"+
                  "failure=root/'fail.txt'\n"+
                  "sys.exit(17 if failure.exists() and failure.read_text()==name else 0)\n")
            for rel in ['envs/data/bin/python','envs/r442_sc/bin/Rscript','envs/r442/bin/Rscript']:
                p=research/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(stub);p.chmod(0o755)
            if fail:(root/'fail.txt').write_text(fail)
            result=subprocess.run([str(RSCRIPT),'--vanilla',str(research/'execute.R'),'--module',module],text=True,capture_output=True)
            calls=(root/'calls.txt').read_text().splitlines() if (root/'calls.txt').exists() else []
            return result,calls

    def test_covariance_runs_all_three_commands(self):
        result,calls=self.run_fixture('covariance_resolution')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(calls,['resolve_shared_outcome.py','diagnose_mvp_conventions.py','assess_givw_heterogeneity.py'])

    def test_failure_stops_later_commands_and_propagates(self):
        result,calls=self.run_fixture('covariance_resolution','diagnose_mvp_conventions.py')
        self.assertEqual(result.returncode,17)
        self.assertEqual(calls,['resolve_shared_outcome.py','diagnose_mvp_conventions.py'])

    def test_final_analysis_and_plot_both_execute(self):
        result,calls=self.run_fixture('final_research')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(calls,['finalize_research.py','plot_final_research.R'])

    def test_single_command_module_still_executes(self):
        result,calls=self.run_fixture('mr_substitute')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(calls,['run_substitute_mr.R'])


if __name__=='__main__':unittest.main()
