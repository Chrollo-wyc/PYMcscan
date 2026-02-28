import math


def ln_fact(x: int) -> float:
    """计算ln(x!)，使用斯特林公式"""
    if x < 12:
        # 小数字直接计算
        result = 0.0
        for i in range(2, x + 1):
            result += math.log(i)
        return result
    else:
        dx = float(x)
        invx = 1.0 / dx
        invx2 = invx * invx
        invx3 = invx2 * invx
        invx5 = invx3 * invx2
        invx7 = invx5 * invx2

        sum_val = ((dx + 0.5) * math.log(dx)) - dx
        sum_val += math.log(2 * math.pi) / 2.0
        sum_val += invx / 12.0 - invx3 / 360.0
        sum_val += invx5 / 1260.0 - invx7 / 1680.0

        return sum_val


def ln_perm(n: int, r: int) -> float:
    """计算排列数的自然对数：ln(P(n,r)) = ln(n!/(n-r)!)"""
    if r > n or r <= 0:
        return 0.0
    return ln_fact(n) - ln_fact(n - r)


def ln_comb(n: int, k: int) -> float:
    """计算组合数的自然对数：ln(C(n,k)) = ln(n!/(k!(n-k)!))"""
    if k <= 0 or k >= n:
        return 0.0
    return ln_fact(n) - ln_fact(k) - ln_fact(n - k)