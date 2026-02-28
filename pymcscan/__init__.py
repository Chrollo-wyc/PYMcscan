# 首先，创建项目的__init__.py文件

# [file name]: mcscanx/__init__.py
"""
MCScanX - Multiple Collinearity Scan Toolbox
Python implementation of the MCScanX genomic collinearity analysis toolkit.
"""

__version__ = "1.0.0"
__author__ = "MCScanX Team"

from pymcscan.core.data_structures import MCScanContext, GeneFeature, BlastRecord, SegFeature
from pymcscan.core.reader import DataReader
from pymcscan.core.dagchainer import DAGChainer
from pymcscan.core.msa import MSA
from pymcscan.core.classifier import GeneClassifier
from pymcscan.output.output_utils import OutputUtils

# 核心模块导出
__all__ = [
    'MCScanContext',
    'GeneFeature',
    'BlastRecord',
    'SegFeature',
    'DataReader',
    'DAGChainer',
    'MSA',
    'GeneClassifier',
    'OutputUtils',
]