# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:18:14 2024

@author: ivanp
"""

from .funcs import calculate_telomere_lengths
from .cli import command_line_target, validate_targets_target

__all__ = ["calculate_telomere_lengths",
           "command_line_target", "validate_targets_target"]
