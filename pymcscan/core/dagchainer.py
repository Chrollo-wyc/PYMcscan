from pymcscan.core.data_structures import *
from pymcscan.utils.math_utils import ln_perm


class DAGChainer:
    """DAG链式算法实现"""

    def __init__(self, context: MCScanContext):
        self.context = context
        self.max_y = 0

    def check_overlap(self, xx: List[int], yy: List[int]) -> bool:
        """检查比对是否重叠"""
        if not xx or not yy:
            return False
        xmin, xmax = min(xx), max(xx)
        ymin, ymax = min(yy), max(yy)
        return xmin <= ymax and ymin <= xmax

    def retrieve_pos(self, pid: int) -> Tuple[int, int]:
        """获取BLAST对的位置"""
        match_rec = self.context.match_list[pid]
        pos1 = self.context.gene_map[match_rec.gene1].mid
        pos2 = self.context.gene_map[match_rec.gene2].mid
        return pos1, pos2

    def is_significant(self, sf: SegFeature, score_list: List[ScoreTuple]) -> bool:
        """判断共线性块是否显著"""
        if not sf.pids:
            return False

        # 获取起止坐标
        s1_a, s1_b = sf.s1.mid, sf.t1.mid
        s2_a, s2_b = sf.s2.mid, sf.t2.mid

        # 计算锚点数
        m = len(sf.pids)

        # 计算区域内的匹配数
        N = 0
        for score_item in score_list:
            if (s1_a <= score_item.x_mid <= s1_b and
                    s2_a <= score_item.y_mid <= s2_b):
                N += 1

        # 计算连续锚点间的距离
        summation = 0.0
        l1_pos1, l2_pos1 = self.retrieve_pos(sf.pids[0])

        for i in range(1, m):
            l1_pos2, l2_pos2 = self.retrieve_pos(sf.pids[i])
            l1 = abs(l1_pos2 - l1_pos1)
            l2 = abs(l2_pos2 - l2_pos1)
            l1_pos1, l2_pos1 = l1_pos2, l2_pos2

            if l1 > 0 and l2 > 0:
                summation += math.log(l1) + math.log(l2)

        # 计算区域长度
        L1 = s1_b - s1_a
        L2 = s2_b - s2_a

        if L1 <= 0 or L2 <= 0:
            sf.e_value = 1.0
            return False

        # 计算e值
        sf.e_value = math.exp(math.log(2) + ln_perm(N, m) +
                              summation - (m - 1) * (math.log(L1) + math.log(L2)))

        return sf.e_value < self.context.E_VALUE

    def check_self(self, s: str) -> bool:
        """检查是否是自比较"""
        if '&' not in s:
            return False
        parts = s.split('&')
        return len(parts) == 2 and parts[0] == parts[1]

    def print_chains(self, score_list: List[ScoreTuple], mol_pair: str):
        """寻找并输出最高评分链"""
        is_self = self.check_self(mol_pair)
        score_list.sort()

        done = False
        while not done:
            n = len(score_list)
            path_score = [0.0] * n
            from_idx = [-1] * n

            # 初始化路径分数
            for i in range(n):
                path_score[i] = score_list[i].score

            # 动态规划寻找最佳路径
            for j in range(1, n):
                for i in range(j - 1, -1, -1):
                    del_x = score_list[j].x - score_list[i].x - 1
                    del_y = score_list[j].y - score_list[i].y - 1

                    if del_x >= 0 and del_y >= 0:
                        if del_x > self.context.MAX_GAPS:
                            break
                        if del_y > self.context.MAX_GAPS:
                            continue

                        num_gaps = max(del_x, del_y)
                        x = path_score[i] + score_list[j].score

                        # 缺口惩罚
                        if num_gaps > 0:
                            x += num_gaps * self.context.GAP_PENALTY

                        if x > path_score[j]:
                            path_score[j] = x
                            from_idx[j] = i

            # 收集高分路径
            high_paths = []
            for i in range(n):
                if path_score[i] >= self.context.CUTOFF_SCORE:
                    high_paths.append(PathTuple(
                        score=path_score[i],
                        rc=score_list[i].x + score_list[i].y,
                        sub=i
                    ))

            # 按分数降序排序
            high_paths.sort(key=lambda x: (-x.score, -x.rc))

            done = True
            for path in high_paths:
                if from_idx[path.sub] != -2:
                    # 重建路径
                    ans = []
                    j = path.sub
                    while from_idx[j] >= 0:
                        ans.append(j)
                        j = from_idx[j]
                    ans.append(j)

                    if from_idx[j] == -2:
                        done = False
                        break

                    # 反转路径
                    ans.reverse()

                    if is_self:
                        xx = []
                        yy = []
                        for idx in ans:
                            from_idx[idx] = -2
                            xx.append(score_list[idx].x)
                            yy.append(score_list[idx].y)

                    # 创建共线性片段
                    if not (is_self and self.check_overlap(xx, yy)):
                        sf = SegFeature()
                        sf.score = path_score[path.sub]

                        for idx in ans:
                            from_idx[idx] = -2
                            pid = score_list[idx].pair_id
                            sf.pids.append(pid)

                        # 设置起止位置
                        br = self.context.match_list[sf.pids[0]]
                        sf.s1 = self.context.gene_map[br.gene1]
                        sf.s2 = self.context.gene_map[br.gene2]

                        br = self.context.match_list[sf.pids[-1]]
                        sf.t1 = self.context.gene_map[br.gene1]
                        sf.t2 = self.context.gene_map[br.gene2]

                        # 确定方向
                        sf.same_strand = sf.s2.mid < sf.t2.mid
                        if not sf.same_strand:
                            sf.s2, sf.t2 = sf.t2, sf.s2

                        sf.mol_pair = mol_pair

                        # 显著性检验
                        if self.is_significant(sf, score_list):
                            self.context.seg_list.append(sf)

            # 更新分数列表
            if not done:
                new_scores = []
                for i in range(n):
                    if from_idx[i] != -2:
                        new_scores.append(score_list[i])
                score_list = new_scores

    def dag_main(self, score_list: List[ScoreTuple], mol_pair: str):
        """DAG主算法"""
        n = len(score_list)
        if n == 0:
            return

        # 按y坐标排序
        score_list.sort(key=lambda x: x.y)
        self.max_y = score_list[-1].y

        # 正向方向
        self.print_chains(score_list.copy(), mol_pair)

        # 反向互补
        for i in range(n):
            score_list[i].y = self.max_y - score_list[i].y + 1

        # 反向方向
        self.print_chains(score_list, mol_pair)

        # 清空分数列表
        self.context.score.clear()