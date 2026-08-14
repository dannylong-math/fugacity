"""Sphinx configuration for the Fugacity documentation."""

from pathlib import Path
import shutil
import subprocess


DOCS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = DOCS_DIR.parent
INCLUDE_DIR = PROJECT_ROOT / "include"
PUBLIC_HEADERS = sorted(INCLUDE_DIR.rglob("*.h")) + sorted(INCLUDE_DIR.rglob("*.hpp"))
API_INPUT = "".join(
    f'#include "{header.relative_to(INCLUDE_DIR).as_posix()}"\n'
    for header in PUBLIC_HEADERS
)
CLANGXX = shutil.which("clang++")
if CLANGXX is None:
    raise RuntimeError("clang++ is required to parse Fugacity's C++ API")
CLANG_RESOURCE_DIR = subprocess.run(
    [CLANGXX, "-print-resource-dir"],
    check=True,
    capture_output=True,
    text=True,
).stdout.strip()

project = "Fugacity"
author = "The Fugacity contributors"
copyright = "2026, The Fugacity contributors"
version = "0.1"
release = "0.1.0"

extensions = [
    "sphinx_immaterial",
    "sphinx_immaterial.apidoc.cpp.apigen",
]

html_theme = "sphinx_immaterial"
html_title = "Fugacity documentation"
html_theme_options = {
    "font": False,
    "features": [
        "navigation.expand",
        "navigation.sections",
        "navigation.top",
        "navigation.footer",
        "toc.follow",
        "toc.sticky",
        "content.code.copy",
        "content.tooltips",
    ],
    "palette": [
        {
            "media": "(prefers-color-scheme: light)",
            "scheme": "default",
            "primary": "blue-grey",
            "accent": "deep-orange",
            "toggle": {
                "icon": "material/brightness-7",
                "name": "Switch to dark mode",
            },
        },
        {
            "media": "(prefers-color-scheme: dark)",
            "scheme": "slate",
            "primary": "blue-grey",
            "accent": "amber",
            "toggle": {
                "icon": "material/brightness-4",
                "name": "Switch to light mode",
            },
        },
    ],
    "toc_title_is_page_title": True,
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
highlight_language = "cpp"
primary_domain = "cpp"

cpp_apigen_configs = [
    {
        "document_prefix": "api/generated/",
        "api_parser_config": {
            "input_content": API_INPUT,
            "compiler_flags": [
                "-std=c++23",
                "-x",
                "c++",
                "-DFUGACITY_DOCUMENTATION=1",
                "-resource-dir",
                CLANG_RESOURCE_DIR,
                "-I",
                str(INCLUDE_DIR),
            ],
            "include_directory_map": {f"{INCLUDE_DIR}/": ""},
            "allow_paths": [r"^fugacity/"],
            "allow_symbols": [r"^fugacity(::.*)?$"],
            "disallow_symbols": [r"^fugacity::detail(::.*)?$"],
            "disallow_namespaces": [r"^std$"],
            "type_replacements": {
                "FUGACITY_ALWAYS_INLINE": "",
                "FUGACITY_RESTRICT": "",
            },
        },
    }
]

cpp_apigen_rst_prolog = """
.. default-role:: cpp:expr

.. default-literal-role::
"""

# The generated signatures intentionally refer to standard-library concepts and
# types that are outside this project's API inventory.
nitpicky = False
