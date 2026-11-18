import json
import os
import shlex
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from gosh_cli.main import cli
from gosh_cli.settings import ALL_TOOLS


class PresetCliFixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.runner = CliRunner()
        self.cwd = os.getcwd()
        os.chdir(self.root)
        self.addCleanup(os.chdir, self.cwd)
        self.env = patch.dict(os.environ, {"XDG_CONFIG_HOME": str(self.root / "config")})
        self.env.start()
        self.addCleanup(self.env.stop)
        self.config = self.root / "config" / "gosh" / "presets.json"
        self.sheet = self.root / "samplesheet.csv"
        self.sheet.write_text("pair,patient,sample,status\nP,P,T,1\nP,P,N,0\n")

    def invoke(self, *args, success=True):
        result = self.runner.invoke(cli, list(args))
        if success:
            self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        else:
            self.assertNotEqual(result.exit_code, 0, result.output)
        return result

    def render(self, preset="ffpe", output="params.json", success=True):
        return self.invoke(
            "run", "params", "--preset", preset,
            "--samplesheet", str(self.sheet), "--outdir", "./results/",
            "--reference", "hg19", "--output", output, success=success,
        )

    def read_json(self, path):
        return json.loads(Path(path).read_text())

    def cli_process(self, *args):
        env = dict(os.environ)
        repo = str(Path(__file__).resolve().parents[2])
        env["PYTHONPATH"] = repo + os.pathsep + env.get("PYTHONPATH", "")
        return subprocess.run(
            [sys.executable, "-B", "-c", "from gosh_cli.main import cli; cli()", *args],
            env=env, text=True, capture_output=True, check=False,
        )


class PresetsTest(PresetCliFixture):
    def test_names_builtin_values_and_typed_process_persistence(self):
        names = self.invoke("preset", "list").output.splitlines()
        self.assertEqual(names, ["default", "ffpe", "tumor-only", "heme", "jabba", "hrd", "taps"])
        expected_skips = {
            "default": "", "ffpe": "", "tumor-only": "",
            "heme": "msisensorpro,hrdetect,onenesstwoness",
            "jabba": "sage,snpeff,snv_multiplicity,signatures,hrdetect",
            "hrd": "non_integer_balance,lp_phased_balance,events,fusions",
            "taps": "non_integer_balance,lp_phased_balance",
        }
        for name, skip in expected_skips.items():
            with self.subTest(name=name):
                shown = json.loads(self.invoke("preset", "show", name).output)
                self.assertEqual(shown["skip_tools"], skip)
        result = self.cli_process("preset", "set", "ffpe", "filter_ffpe_chimera", "false")
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.cli_process("preset", "show", "ffpe")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIs(json.loads(result.stdout)["filter_ffpe_chimera"], False)
        self.assertEqual(self.read_json(self.config), {"ffpe": {"filter_ffpe_chimera": False}})
        self.assertEqual(stat.S_IMODE(self.config.stat().st_mode), 0o600)
        self.config.chmod(0o640)
        self.invoke("preset", "set", "heme", "purple_min_purity", "0.3")
        self.invoke("preset", "set", "heme", "custom_count", "2")
        self.invoke("preset", "set", "heme", "custom_enabled", "true")
        self.invoke("preset", "set", "heme", "custom_values", '[1,true,null,"a"]')
        shown = json.loads(self.invoke("preset", "show", "heme").output)
        self.assertEqual(shown["purple_min_purity"], 0.3)
        self.assertIs(type(shown["custom_count"]), int)
        self.assertIs(shown["custom_enabled"], True)
        self.assertEqual(shown["custom_values"], [1, True, None, "a"])
        self.assertEqual(stat.S_IMODE(self.config.stat().st_mode), 0o640)
        self.invoke("preset", "unset", "ffpe", "filter_ffpe_chimera")
        self.invoke("preset", "unset", "ffpe", "absent")
        self.assertIs(json.loads(self.invoke("preset", "show", "ffpe").output)["filter_ffpe_chimera"], True)

    def test_invalid_writes_and_hand_edits_never_rewrite_config(self):
        self.invoke("preset", "set", "ffpe", "filter_ffpe_chimera", "false")
        original = self.config.read_bytes()
        invalid = [
            ("missing", "x", "true"), ("ffpe", "", "true"),
            ("ffpe", "filter_ffpe_chimera", "1"),
            ("ffpe", "filter_ffpe_chimera", '"true"'),
            ("heme", "purple_min_purity", "true"),
            ("ffpe", "skip_tools", '"not_a_tool"'),
            ("ffpe", "skip_tools", "[]"),
            ("ffpe", "x", "True"), ("ffpe", "x", "NaN"),
            ("ffpe", "x", '{"api_key":"value"}'),
        ]
        invalid += [("ffpe", key, '"value"') for key in (
            "input", "outdir", "genome", "email", "tumor_only",
            "ONCOKB_API_KEY", "SecretName", "PASSWORD", "accessToken",
        )]
        for name, key, value in invalid:
            with self.subTest(key=key, value=value):
                self.invoke("preset", "set", name, key, value, success=False)
                self.assertEqual(self.config.read_bytes(), original)
        for content in ("{broken", "[]", '{"ffpe":[]}', '{"unknown":{}}',
                        '{"ffpe":{"tumor_only":true}}',
                        '{"ffpe":{"filter_ffpe_chimera":1}}',
                        '{"ffpe":{"ONCOKB_API_KEY":"private"}}'):
            with self.subTest(content=content):
                self.config.write_text(content)
                self.invoke("preset", "show", "ffpe", success=False)
                self.invoke("preset", "set", "ffpe", "filter_ffpe_chimera", "false", success=False)
                self.invoke("preset", "unset", "ffpe", "filter_ffpe_chimera", success=False)
                self.render(output="invalid.json", success=False)
                self.assertEqual(self.config.read_text(), content)
                self.assertFalse(Path("invalid.json").exists())

    def test_params_snapshots_modes_and_no_overwrite(self):
        self.render(output="snapshots/A.json")
        original = Path("snapshots/A.json").read_bytes()
        a = self.read_json("snapshots/A.json")
        self.assertIs(a["tumor_only"], False)
        self.assertIs(a["purple_use_svs"], True)
        self.assertEqual(a["genome"], "GATK.GRCh37")
        self.assertEqual(a["email"], "")
        self.assertEqual(a["outdir"], "./results/")
        for key in ("filter_ffpe_impact", "filter_ffpe_chimera", "sv_filter_ffpe_chimera"):
            self.assertIs(a[key], True)
        self.invoke("preset", "set", "ffpe", "filter_ffpe_chimera", "false")
        self.render(output="snapshots/B.json")
        self.assertIs(self.read_json("snapshots/B.json")["filter_ffpe_chimera"], False)
        self.assertEqual(Path("snapshots/A.json").read_bytes(), original)
        self.render(output="snapshots/A.json", success=False)
        self.assertEqual(Path("snapshots/A.json").read_bytes(), original)
        self.sheet.write_text("patient,sample,status\nH,H,1\n")
        self.render("heme", "heme.json")
        h = self.read_json("heme.json")
        self.assertIs(h["tumor_only"], True)
        self.assertIs(h["is_heme"], True)
        self.assertIs(h["is_retier_whitelist_junctions"], True)
        self.assertIs(h["purple_use_svs"], False)
        self.assertIs(h["purple_use_smlvs"], False)
        self.assertEqual(h["purple_highly_diploid_percentage"], 1.0)
        self.assertEqual(h["purple_min_purity"], 0.25)
        self.assertEqual(h["purple_ploidy_penalty_factor"], 0.6)
        self.assertIs(h["purple_revise_purity_ploidy"], True)
        self.assertNotIn("filter_ffpe_impact", h)
        self.render("tumor-only", "tumor.json")
        t = self.read_json("tumor.json")
        self.assertIs(t["tumor_only"], True)
        self.assertNotIn("is_heme", t)
        self.assertNotIn("filter_ffpe_impact", t)
        self.render("ffpe", "singleton-ffpe.json")
        self.assertIs(self.read_json("singleton-ffpe.json")["tumor_only"], True)

    def test_every_case_is_validated_and_pair_precedes_patient(self):
        cases = [
            "patient,sample,status\n",
            "patient,sample,status\nP,T,2\n",
            "patient,sample,status\nP,T,true\n",
            "patient,sample,status\nP,N,0\n",
            "patient,sample,status\nP,T,1\nP,N,0\nS,S,1\n",
            "pair,patient,sample,status\nP,shared,T,1\nP,shared,N,0\nS,shared,S,1\n",
            "patient,sample,status\nP,T,1\nP,N,0\nQ,N2,0\n",
            "sample,status\nT,1\n",
            "patient,sample,status\n,T,1\n",
        ]
        for contents in cases:
            with self.subTest(contents=contents):
                self.sheet.write_text(contents)
                self.render(output="bad.json", success=False)
                self.assertFalse(Path("bad.json").exists())
        self.sheet.write_text("pair,patient,sample,status\n,P,T,1\n,P,N,0\n")
        self.render(output="fallback.json")
        self.assertIs(self.read_json("fallback.json")["tumor_only"], False)
        self.render("tumor-only", "bad.json", success=False)
        self.assertFalse(Path("bad.json").exists())

    def test_wizard_reuses_preset_then_applies_custom_values(self):
        from gosh_cli.core.params_wizard import create_params_file
        Path("results").mkdir()
        self.invoke("preset", "set", "ffpe", "filter_ffpe_chimera", "false")
        answers = ["", "", "y", "filter_ffpe_chimera", "true", "custom_number", "3", "", "y"]
        with patch("builtins.input", side_effect=answers):
            create_params_file(preset="ffpe", samplesheet=str(self.sheet))
        data = self.read_json("params.json")
        self.assertIs(data["filter_ffpe_chimera"], True)
        self.assertIs(data["sv_filter_ffpe_chimera"], True)
        self.assertEqual(data["custom_number"], 3)


class PipelinePresetTest(PresetCliFixture):
    def setUp(self):
        super().setUp()
        self.defaults = patch("gosh_cli.cli.run.get_environment_defaults", return_value={"profile": "nyu"})
        self.defaults.start()
        self.addCleanup(self.defaults.stop)
        self.modules = patch("gosh_cli.core.module_loader.load_required_modules", return_value="")
        self.load_modules = self.modules.start()
        self.addCleanup(self.modules.stop)
        self.execution = patch("gosh_cli.core.nextflow.NextflowRunner.run")
        self.run_pipeline = self.execution.start()
        self.addCleanup(self.execution.stop)
        self.secret_env = patch.dict(os.environ)
        self.secret_env.start()
        self.addCleanup(self.secret_env.stop)
        os.environ.pop("ONCOKB_API_KEY", None)
        self.render()
        self.args = [
            "run", "pipeline", "--non-interactive", "--pipeline-dir", "/probe pipeline",
            "--samplesheet", str(self.sheet), "--params-file", "params.json",
            "--profile", "singularity,nyu",
        ]

    def pipeline(self, *extra, success=True):
        return self.invoke(*self.args, *extra, success=success)

    def list_secret(self, output="ONCOKB_API_KEY\nONCOKB_TOKEN\n"):
        return patch("gosh_cli.cli.run.run", return_value=subprocess.CompletedProcess([], 0, output, ""))

    def test_existing_snapshot_ignores_default_and_changed_preset_config(self):
        self.invoke("preset", "set", "ffpe", "filter_ffpe_chimera", "false")
        self.invoke("preset", "set", "default", "filter_ffpe_impact", "false")
        # Default execution need not even read the mutable preset store.
        self.config.write_text("{broken")
        with self.list_secret():
            result = self.pipeline()
        saved = self.read_json("params.json")
        self.assertIs(saved["filter_ffpe_chimera"], True)
        self.assertIs(saved["filter_ffpe_impact"], True)
        command = self.run_pipeline.call_args_list[0].args[0]
        tokens = shlex.split(command)
        self.assertEqual(tokens[tokens.index("-profile") + 1], "singularity,nyu")
        self.assertEqual(tokens[tokens.index("run") + 1], "/probe pipeline")
        self.assertIn("-resume", tokens)
        self.assertNotIn("secrets set", result.output)
        self.assertEqual([call.args[0] for call in self.run_pipeline.call_args_list[1:]], ["gosh run post", "gosh protect"])

    def test_explicit_presets_and_tool_flags_override_snapshot(self):
        self.invoke("preset", "set", "ffpe", "filter_ffpe_chimera", "false")
        with self.list_secret():
            self.pipeline("--preset", "ffpe", "--skip-tools", "sage")
        saved = self.read_json("params.json")
        self.assertIs(saved["filter_ffpe_chimera"], False)
        self.assertEqual(saved["skip_tools"], "sage")
        with self.list_secret():
            self.pipeline("--preset", "heme", "--include-tools", "sage,oncokb")
        saved = self.read_json("params.json")
        self.assertIs(saved["is_heme"], True)
        self.assertEqual(set(saved["skip_tools"].split(",")), set(ALL_TOOLS) - {"sage", "oncokb"})
        with self.list_secret():
            self.pipeline("--preset", "heme", "--skip-tools", "")
        self.assertEqual(self.read_json("params.json")["skip_tools"], "")
        self.pipeline("--skip-tools", "", "--include-tools", "sage", success=False)

    def test_invalid_explicit_preset_config_blocks_execution(self):
        self.invoke("preset", "set", "ffpe", "filter_ffpe_chimera", "false")
        self.config.write_text('{"ffpe":{"input":"other.csv"}}')
        original = Path("params.json").read_bytes()
        self.pipeline("--preset", "ffpe", success=False)
        self.pipeline("--preset", "tumor-only", success=False)
        self.assertEqual(Path("params.json").read_bytes(), original)
        self.run_pipeline.assert_not_called()

    def test_noninteractive_requires_explicit_inputs_and_never_falls_back(self):
        with patch("gosh_cli.core.params_wizard.create_params_file") as wizard:
            for args in (
                ["run", "pipeline", "--non-interactive"],
                ["run", "pipeline", "--non-interactive", "--samplesheet", str(self.sheet)],
                ["run", "pipeline", "--non-interactive", "--params-file", "params.json"],
                self.args + ["--params-file", "missing.json"],
                self.args + ["--samplesheet", "missing.csv"],
                self.args + ["--processes", "SAGE"],
                self.args + ["--samples", "T"],
            ):
                with self.subTest(args=args):
                    self.invoke(*args, success=False)
            wizard.assert_not_called()
        self.run_pipeline.assert_not_called()

    def test_noninteractive_revalidates_all_statuses_and_stored_mode(self):
        original = Path("params.json").read_bytes()
        self.sheet.write_text("patient,sample,status\nP,T,1\nP,N,0\nS,S,1\n")
        self.pipeline(success=False)
        self.assertEqual(Path("params.json").read_bytes(), original)
        self.sheet.write_text("patient,sample,status\nP,T,1\n")
        result = self.pipeline(success=False)
        self.assertIn("inconsistent", result.output)
        self.assertEqual(Path("params.json").read_bytes(), original)
        data = self.read_json("params.json")
        data["tumor_only"] = 1
        Path("params.json").write_text(json.dumps(data))
        self.pipeline(success=False)
        self.run_pipeline.assert_not_called()

    def test_secret_requires_exact_name_and_sanitizes_failures(self):
        for output in ("ONCOKB_TOKEN\n", "NOT_ONCOKB_API_KEY\n", "ONCOKB_API_KEY_EXTRA\n", ""):
            with self.subTest(output=output), self.list_secret(output):
                result = self.pipeline(success=False)
                self.assertIn("Nextflow secret ONCOKB_API_KEY is required", result.output)
        with patch("gosh_cli.cli.run.run", side_effect=subprocess.CalledProcessError(1, ["nextflow", "secrets", "list"], stderr="sensitive details")):
            result = self.pipeline(success=False)
            self.assertNotIn("sensitive details", result.output)
        with patch("gosh_cli.cli.run.run", side_effect=subprocess.TimeoutExpired(["nextflow", "secrets", "list"], 60, output="private output")):
            result = self.pipeline(success=False)
            self.assertIn("Nextflow secret ONCOKB_API_KEY is required", result.output)
            self.assertNotIn("private output", result.output)
        key = "private-value; $(touch should-not-exist)"
        failure = subprocess.CalledProcessError(1, ["nextflow", "secrets", "set", "ONCOKB_API_KEY", key], stderr=key)
        with patch("gosh_cli.cli.run.run", side_effect=failure):
            result = self.pipeline("--oncokb-api-key", key, success=False)
        self.assertIn("Unable to provision", result.output)
        self.assertNotIn(key, result.output)
        self.assertNotIn(key, Path("params.json").read_text())
        self.run_pipeline.assert_not_called()

    def test_secret_provisioning_uses_positional_args_in_module_shell(self):
        # Use a real local process at this boundary: the key contains shell syntax.
        executable = self.root / "nextflow"
        executable.write_text(
            "#!" + sys.executable + "\n"
            "import json, sys\n"
            "from pathlib import Path\n"
            "if sys.argv[1:3] == ['secrets', 'set']:\n"
            "    Path('secret-call.json').write_text(json.dumps(sys.argv[1:]))\n"
            "elif sys.argv[1:] == ['secrets', 'list']:\n"
            "    print('ONCOKB_API_KEY')\n"
            "else:\n"
            "    sys.exit(2)\n"
        )
        executable.chmod(0o700)
        self.load_modules.return_value = "true && "
        key = "secret ' value; $(touch INJECTION)"
        with patch.dict(os.environ, {"PATH": str(self.root) + os.pathsep + os.environ["PATH"]}):
            result = self.pipeline("--oncokb-api-key", key)
            self.assertEqual(self.read_json("secret-call.json"), ["secrets", "set", "ONCOKB_API_KEY", key])
            self.assertFalse(Path("INJECTION").exists())
            self.assertNotIn(key, result.output)
            self.assertNotIn(key, Path("params.json").read_text())
            self.assertNotIn(key, self.run_pipeline.call_args_list[0].args[0])
            Path("secret-call.json").unlink()
            with patch.dict(os.environ, {"ONCOKB_API_KEY": key}):
                result = self.pipeline()
            self.assertEqual(self.read_json("secret-call.json")[-1], key)
            self.assertNotIn(key, result.output)
            Path("secret-call.json").unlink()
            self.pipeline()
            self.assertFalse(Path("secret-call.json").exists())

    def test_monitor_host_and_duplicate_destination_fail_before_effects(self):
        original = Path("params.json").read_bytes()
        for host in ("http://host", "user@host", "a/b", "a:80", "a b", "$(id)", "-host", "a;id"):
            with self.subTest(host=host):
                self.pipeline("--gos-eye", "8998", "--gos-eye-host", host, success=False)
        self.pipeline("--gos-eye-host", "monitor.cluster", success=False)
        for flag in ("-with-weblog", "-with-weblog=http://other/webhook"):
            self.pipeline("--gos-eye", "8998", "--", flag, success=False)
        self.assertEqual(Path("params.json").read_bytes(), original)
        self.assertFalse(Path("hg19").exists())
        self.run_pipeline.assert_not_called()

    def test_preview_rejects_execution_options_before_effects(self):
        original = Path("params.json").read_bytes()
        options = (["--bg"], ["--processes", "X"], ["--samples", "S"],
                   ["--preset", "heme"], ["--skip-tools", ""], ["--include-tools", "sage"],
                   ["--gos-eye", "8998"], ["--oncokb-api-key", "fake"])
        with patch("gosh_cli.core.params_wizard.create_params_file") as wizard, self.list_secret() as secrets:
            for flags in options:
                self.pipeline(*flags, "--", "-preview", success=False)
            for flag in ("-resume", "-bg", "-with-weblog", "-with-report", "-with-trace", "-with-timeline", "-with-notification", "-N"):
                for token in (flag, flag + "=value"):
                    self.pipeline("--", "-preview", token, success=False)
            wizard.assert_not_called()
            secrets.assert_not_called()
        self.assertEqual(Path("params.json").read_bytes(), original)
        self.assertFalse(Path("hg19").exists())
        self.run_pipeline.assert_not_called()

    def test_preview_failure_preserves_snapshot_and_removes_temporary_params(self):
        original = Path("params.json").read_bytes()
        sheet_bytes = self.sheet.read_bytes()
        executable = self.root / "nextflow"
        executable.write_text("#!/bin/sh\nexit 37\n")
        executable.chmod(0o700)
        with patch.dict(os.environ, {"PATH": str(self.root) + os.pathsep + os.environ["PATH"]}), patch("gosh_cli.core.params_wizard.create_params_file") as wizard:
            result = self.pipeline("--outdir", "./override/", "--", "-preview", "-with-dag", "dag.dot", success=False)
            self.assertEqual(result.exit_code, 37, result.output)
            wizard.assert_not_called()
        self.assertEqual(Path("params.json").read_bytes(), original)
        self.assertEqual(self.sheet.read_bytes(), sheet_bytes)
        self.assertFalse(list(self.root.glob(".gosh-preview-*")))
        for name in ("hg19", "dag.dot", "outputs.csv", "override"):
            self.assertFalse(Path(name).exists(), name)
        self.run_pipeline.assert_not_called()

    def test_preview_without_noninteractive_still_validates_mode_and_never_prompts(self):
        Path("params.json").unlink()
        with patch("gosh_cli.core.params_wizard.create_params_file") as wizard:
            self.invoke("run", "pipeline", "--", "-preview", success=False)
            wizard.assert_not_called()
        self.render()
        self.sheet.write_text("patient,sample,status\nP,T,1\n")
        original = Path("params.json").read_bytes()
        result = self.invoke("run", "pipeline", "--", "-preview", success=False)
        self.assertIn("inconsistent", result.output)
        self.assertEqual(Path("params.json").read_bytes(), original)
        self.run_pipeline.assert_not_called()


if __name__ == "__main__":
    unittest.main()
