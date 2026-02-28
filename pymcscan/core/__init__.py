# [file name]: mcscanx/core/__init__.py
"""
MCScanX核心模块
包含主要的数据结构和算法实现。
"""

from .data_structures import *
from .reader import DataReader
from .dagchainer import DAGChainer
from .msa import MSA
from .classifier import GeneClassifier

__all__ = [
    'DataReader',
    'DAGChainer',
    'MSA',
    'GeneClassifier',
    'MCScanContext',
    'GeneFeature',
    'BlastRecord',
    'SegFeature',
    'ScoreTuple',
    'PathTuple',
    'MoreFeature',
    'OrthoStat',
]