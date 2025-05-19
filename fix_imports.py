"""Script to fix imports in Python files."""

import os
from pathlib import Path


def fix_imports(directory):
    """Fix imports in all Python files in the directory."""
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".py"):
                filepath = Path(root) / file
                fix_file_imports(filepath)


def fix_file_imports(filepath):
    """Fix imports in a single Python file."""
    with open(filepath) as f:
        content = f.read()

    # Replace direct _helpers imports with relative imports
    if "from _helpers import" in content:
        # Calculate relative path to utils/helpers.py
        rel_path = Path(filepath).relative_to(Path("src/pypsa_usa"))
        dots = len(rel_path.parents) - 1
        prefix = "." * dots if dots > 0 else "."

        # Replace the import
        new_content = content.replace(
            "from _helpers import",
            f"from {prefix}utils.helpers import",
        )

        with open(filepath, "w") as f:
            f.write(new_content)
        print(f"Fixed imports in {filepath}")


if __name__ == "__main__":
    fix_imports("src/pypsa_usa")
