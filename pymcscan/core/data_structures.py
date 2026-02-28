from dataclasses import dataclass, field
from typing import List, Set, Dict, Optional, Tuple
import math
from collections import defaultdict


@dataclass
class GeneFeature:
    """基因特征"""
    name: str
    mol: str  # 分子/染色体标识
    mid: int  # 中间位置
    gene_id: int = 0
    cursor: List[int] = field(default_factory=list)
    in_blocks: int = 0
    cr_blocks: int = 0
    sp: Set[str] = field(default_factory=set)

    def __lt__(self, other: 'GeneFeature') -> bool:
        """排序比较：先按mol，再按mid，最后按name"""
        return (self.mol < other.mol) or \
            (self.mol == other.mol and self.mid < other.mid) or \
            (self.mol == other.mol and self.mid == other.mid and self.name < other.name)


@dataclass
class BlastRecord:
    """BLAST比对记录"""
    gene1: str
    gene2: str
    mol_pair: str
    pair_id: int
    score: float


@dataclass
class SegFeature:
    """片段特征（共线性块）"""
    pids: List[int] = field(default_factory=list)
    index: int = 0
    s1: Optional[GeneFeature] = None
    t1: Optional[GeneFeature] = None
    s2: Optional[GeneFeature] = None
    t2: Optional[GeneFeature] = None
    score: float = 0.0
    e_value: float = 0.0
    mol_pair: str = ""
    same_strand: bool = True


@dataclass
class ScoreTuple:
    """评分元组"""
    pair_id: int
    x: int  # gene1在基因列表中的索引
    y: int  # gene2在基因列表中的索引
    x_mid: int = 0
    y_mid: int = 0
    score: float = 0.0
    gene1: str = ""
    gene2: str = ""

    def __lt__(self, other: 'ScoreTuple') -> bool:
        """排序比较：先按x，再按y"""
        return (self.x < other.x) or (self.x == other.x and self.y < other.y)


@dataclass
class PathTuple:
    """路径元组"""
    score: float
    rc: int  # 行列和
    sub: int


@dataclass
class MoreFeature:
    """附加特征"""
    tandem: int = 0
    depth: int = 0


@dataclass
class OrthoStat:
    """直系同源统计"""
    all_num: int = 0
    syn_num: int = 0


class MCScanContext:
    """MCScan上下文，管理所有全局状态"""

    def __init__(self):
        # 全局数据
        self.gene_map: Dict[str, GeneFeature] = {}
        self.match_list: List[BlastRecord] = []
        self.seg_list: List[SegFeature] = []
        self.mol_pairs: Dict[str, int] = defaultdict(int)
        self.allg: List[GeneFeature] = []
        self.cmp_sp: Dict[str, OrthoStat] = defaultdict(OrthoStat)
        self.gene_more: List[MoreFeature] = []
        self.score: List[ScoreTuple] = []

        # 算法参数（默认值）
        self.MATCH_SCORE = 50
        self.MATCH_SIZE = 5
        self.GAP_PENALTY = -1
        self.GAP_SIZE = 5
        self.OVERLAP_WINDOW = 5
        self.E_VALUE = 1e-5
        self.MAX_GAPS = 25
        self.N_PROXIMAL = 10
        self.CUTOFF_SCORE = self.MATCH_SCORE * self.MATCH_SIZE
        self.IN_SYNTENY = 0
        self.e_mode = 0

    def fill_allg(self):
        """填充allg基因列表"""
        self.allg = sorted(self.gene_map.values())
        for i, gene in enumerate(self.allg):
            gene.gene_id = i

    def set_cutoff(self):
        """设置截断分数"""
        self.CUTOFF_SCORE = self.MATCH_SCORE * self.MATCH_SIZE