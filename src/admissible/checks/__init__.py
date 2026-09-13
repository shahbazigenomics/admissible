"""The five checks.  Each is importable and runnable on its own."""

from .callability import check_callability
from .genotype import check_genotype
from .identity import IdentityConfig, check_identity
from .models import check_models
from .provenance import check_provenance

__all__ = [
    "IdentityConfig",
    "check_identity",
    "check_genotype",
    "check_provenance",
    "check_callability",
    "check_models",
]
