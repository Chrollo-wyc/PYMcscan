# [file name]: mcscanx/cli/__init__.py
"""
MCScanX命令行接口模块
包含所有命令行工具的主程序。
"""

from .mcscan import main as mcscan_main
from .dup_classifier import main as dup_classifier_main

__all__ = [
    'mcscan_main',
    'dup_classifier_main',
]