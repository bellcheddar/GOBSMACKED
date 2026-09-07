"""How run.py leaves.

One test, for one expensive regression. The runner imports torch (via pandadock)
and OpenMM at module scope whether or not their stages run, and both leave
non-daemon threads behind. `sys.exit` waits for every one of them: on a real run
that wait was 25 minutes between the last line printed and the process actually
ending, which every watchdog outside reads as a hung stage rather than a
finished one.

Checked on the syntax tree rather than by grepping for a string, so a comment
mentioning os._exit cannot satisfy it and neither can a docstring.
"""

from __future__ import annotations

import ast
from pathlib import Path

RUN_PY = Path(__file__).resolve().parents[1] / "bundle_template" / "run.py"


def _main_block(tree: ast.Module) -> ast.If:
    for node in tree.body:
        if (isinstance(node, ast.If) and isinstance(node.test, ast.Compare)
                and isinstance(node.test.left, ast.Name)
                and node.test.left.id == "__name__"):
            return node
    raise AssertionError("run.py has no `if __name__ == '__main__'` block")


def _calls(node: ast.AST) -> set[str]:
    found = set()
    for child in ast.walk(node):
        if isinstance(child, ast.Call):
            found.add(ast.unparse(child.func))
    return found


def test_the_runner_leaves_without_waiting_for_library_threads():
    block = _main_block(ast.parse(RUN_PY.read_text(encoding="utf-8")))
    calls = _calls(block)
    assert "os._exit" in calls, (
        "run.py must leave with os._exit: sys.exit waits for torch's and "
        "OpenMM's non-daemon threads, which cost 25 minutes on a real run")
    assert "sys.exit" not in calls
    # os._exit skips Python's own buffers, so the flush has to be explicit or
    # the last lines of the run never reach the log.
    assert any(c.endswith(".flush") for c in calls), "flush before os._exit"


def test_the_exit_code_is_the_one_main_returned():
    """A hard exit that always reported 0 would be worse than the hang: every
    supervisor keys off the return code."""
    block = _main_block(ast.parse(RUN_PY.read_text(encoding="utf-8")))
    exits = [c for c in ast.walk(block)
             if isinstance(c, ast.Call) and ast.unparse(c.func) == "os._exit"]
    assert len(exits) == 1
    assert [ast.unparse(a) for a in exits[0].args] == ["code"]
    assigned = [ast.unparse(n.value) for n in ast.walk(block)
                if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == "code"]
    assert assigned == ["main()"]
