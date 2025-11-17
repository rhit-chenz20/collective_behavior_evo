#!/usr/bin/env python3
"""
Make 4 SLiM scripts by replacing label lines in a base script.

Labels expected (as their own lines):
    ######PHENOTYPE_CALCULATION######
    ######FITNESS_ADJUSTMENT######

Replacement preserves the leading indentation of the label line.
"""

from pathlib import Path
import re
from typing import Tuple

# --- Internal helpers --------------------------------------------------------

PHENO_LABEL = r"######PHENOTYPE_CALCULATION######"
FIT_LABEL   = r"######FITNESS_ADJUSTMENT######"

# Regex to match the entire label line, capturing leading indentation.
PHENO_PATTERN = re.compile(rf"^([ \t]*){re.escape(PHENO_LABEL)}[ \t]*$", re.MULTILINE)
FIT_PATTERN   = re.compile(rf"^([ \t]*){re.escape(FIT_LABEL)}[ \t]*$",   re.MULTILINE)

def _indent_block(indent: str, text: str) -> str:
    """Indent every non-empty line of `text` by `indent`."""
    lines = text.splitlines()
    return "\n".join((indent + ln if ln.strip() != "" else ln) for ln in lines)

def replace_labels(
    content: str,
    phenotype_code: str,
    fitness_code: str,
) -> str:
    """Replace label lines with provided code blocks, keeping indentation."""
    # Replace PHENOTYPE
    def _pheno_repl(m: re.Match) -> str:
        indent = m.group(1)
        return _indent_block(indent, phenotype_code)

    # Replace FITNESS
    def _fit_repl(m: re.Match) -> str:
        indent = m.group(1)
        return _indent_block(indent, fitness_code)

    # Do replacements, but also verify presence
    if not PHENO_PATTERN.search(content):
        raise ValueError(f"Label not found: {PHENO_LABEL}")
    if not FIT_PATTERN.search(content):
        raise ValueError(f"Label not found: {FIT_LABEL}")

    content = PHENO_PATTERN.sub(_pheno_repl, content)
    content = FIT_PATTERN.sub(_fit_repl, content)
    return content

def make_variant(
    base_path: Path,
    out_path: Path,
    phenotype_code: str,
    fitness_code: str,
) -> None:
    """Create a single variant file from base + provided replacement code."""
    base_txt = base_path.read_text(encoding="utf-8")
    new_txt = replace_labels(base_txt, phenotype_code, fitness_code)
    out_path.write_text(new_txt, encoding="utf-8")
    print(f"Wrote: {out_path}")

# --- Fill these 4 variants and run -------------------------------------------

if __name__ == "__main__":
    # Point to your base SLiM file:
    BASE = Path("base_script.slim")

    # ↓↓↓ Fill these blocks with your actual SLiM code (multiline strings OK) ↓↓↓
    # Variant 1
    PHENO_1 = """// phenotype for ave model
phenotype[group.index] = ((n-1)*(1-psi)*a[group.index]+psi*sum(as)) / (n-1+psi);
"""
    FIT_1 = """// no fitness adjustment needed
adj_fitness = abs_fitness;
"""

    # Variant 2
    PHENO_2 = """// phenotype for best model
fit = dnorm(as, optimum, sqrt(Vs));
best_fit_i = whichMax(fit);
best_fit_a = as[best_fit_i];
phenotype[group.index] = (1-psi) * as + psi * best_fit_a;
"""
    FIT_2 = """// no fitness adjustment needed
adj_fitness = abs_fitness;
"""

    # Variant 3
    PHENO_3 = """// phenotype for extreme model
mean_a = mean(as);
dist = abs(as-mean_a);
extreme_i = whichMax(dist);
extreme_a = as[extreme_i];
phenotype[group.index] = (1-psi) * as + psi * extreme_a;
"""
    FIT_3 = """// no fitness adjustment needed
adj_fitness = abs_fitness;
"""

    # Variant 4
    PHENO_4 = """// no phenotype for inverse-variance model; record group variance to tagF
phenotype[group.index] = a[group.index];
a_var = var(as);
group.tagF = a_var;
"""
    FIT_4 = """// adjust the fitness of individuals
adj_fitness = abs_fitness / (1 + psi * inds.tagF);
"""

    # Choose output filenames:
    OUT_1 = Path("ave.slim")
    OUT_2 = Path("fit.slim")
    OUT_3 = Path("ext.slim")
    OUT_4 = Path("inv.slim")

    # Generate all 4
    make_variant(BASE, OUT_1, PHENO_1, FIT_1)
    make_variant(BASE, OUT_2, PHENO_2, FIT_2)
    make_variant(BASE, OUT_3, PHENO_3, FIT_3)
    make_variant(BASE, OUT_4, PHENO_4, FIT_4)
