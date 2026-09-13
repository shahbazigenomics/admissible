"""admissible - audit whether exome data can support an interpretation.

This tool does not interpret variants.  It audits the evidence base before
interpretation, and answers one question: *is this evidence admissible?*
"""

from .model import CheckResult, Finding, Report, Severity, Status, Verdict

__version__ = "0.1.0"
__all__ = ["Report", "CheckResult", "Finding", "Status", "Severity", "Verdict", "__version__"]
