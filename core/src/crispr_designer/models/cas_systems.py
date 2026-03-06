'''
CRISPR Cas systems definitions
Each nuclease system is defined with it's PAM sequence, guide length constraints, and cleavage characteristics.
'''

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import ClassVar

