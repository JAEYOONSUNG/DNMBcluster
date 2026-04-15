"""Smoke tests for the M0 CLI scaffold."""
from click.testing import CliRunner

from dnmbcluster.cli import main


def test_help() -> None:
    result = CliRunner().invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "DNMBcluster" in result.output


def test_list_engines() -> None:
    result = CliRunner().invoke(main, ["list-engines"])
    assert result.exit_code == 0
    assert "mmseqs2" in result.output
    assert "usearch12" in result.output
