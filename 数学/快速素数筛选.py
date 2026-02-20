#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
================================================================================
五大素数生成算法对比分析（含内存占用与结果验证）
Comparative Analysis of Five Prime Generation Algorithms (with memory & validation)

算法列表：
1. 试除法 (Trial Division) - 基础验证法
2. 标准筛 (Standard Sieve) - Eratosthenes筛
3. 分段筛 (Segmented Sieve) - 内存优化版
4. 标准Wheel (Standard Wheel) - 数组实现的30-wheel
5. 自适应wheel算法 (Adaptive Wheel) - bytearray实现的30030-wheel

对比维度：
- 时间复杂度
- 空间复杂度
- 实际运行时间
- 峰值内存占用
- 结果正确性（素数个数一一对比验证）
- 适用场景
================================================================================
"""

import time
import math
import array
import tracemalloc
from typing import List, Tuple


# ==================== 算法1：试除法 ====================
def trial_division(n: int) -> List[int]:
    """试除法 - O(n√n)时间，O(π(n))空间"""
    if n < 2:
        return []
    primes = []
    for x in range(2, n + 1):
        is_prime = True
        limit = int(math.sqrt(x))
        for p in primes:
            if p > limit:
                break
            if x % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(x)
    return primes


# ==================== 算法2：标准筛 ====================
def standard_sieve(n: int) -> List[int]:
    """Eratosthenes筛 - O(n log log n)时间，O(n)空间"""
    if n < 2:
        return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(n)) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False
    return [i for i in range(n + 1) if is_prime[i]]


# ==================== 算法3：分段筛 ====================
def segmented_sieve(n: int) -> List[int]:
    """分段筛 - O(n log log n)时间，O(√n)空间"""
    if n < 2:
        return []
    sqrt_n = int(math.sqrt(n)) + 1

    base_primes = []
    is_comp = [False] * (sqrt_n + 1)
    for i in range(2, sqrt_n):
        if not is_comp[i]:
            base_primes.append(i)
            for j in range(i * i, sqrt_n, i):
                is_comp[j] = True

    if n <= sqrt_n:
        return [p for p in base_primes if p <= n]

    primes = [p for p in base_primes if p <= n]
    seg_size = 32768

    for low in range(sqrt_n, n + 1, seg_size):
        high = min(low + seg_size - 1, n)
        sieve = array.array("b", [1]) * (high - low + 1)
        for p in base_primes:
            if p * p > high:
                break
            start = ((low + p - 1) // p) * p
            for j in range(start, high + 1, p):
                if j >= low:
                    sieve[j - low] = 0
        for i in range(high - low + 1):
            if sieve[i]:
                primes.append(low + i)
    return primes


# ==================== 算法4：标准Wheel ====================
def standard_wheel_sieve(n: int) -> List[int]:
    """30-wheel筛 - O(n log log n)时间，O(n)空间，标记减少73%"""
    if n < 2:
        return []
    small = [2, 3, 5]
    if n <= 5:
        return [p for p in small if p <= n]

    sqrt_n = int(math.sqrt(n)) + 1
    base = []
    is_comp = [False] * (sqrt_n + 1)
    for i in range(2, sqrt_n):
        if not is_comp[i]:
            base.append(i)
            for j in range(i * i, sqrt_n, i):
                is_comp[j] = True

    primes = small[:]
    wheel = [1, 7, 11, 13, 17, 19, 23, 29]
    is_prime = [False] * (n + 1)

    for k in range(n // 30 + 1):
        for r in wheel:
            x = 30 * k + r
            if x <= n and x > 5:
                is_prime[x] = True

    for p in base:
        if p <= 5:
            continue
        for j in range(p * p, n + 1, p):
            if is_prime[j]:
                is_prime[j] = False

    for x in range(6, n + 1):
        if is_prime[x]:
            primes.append(x)
    return primes


# ==================== 算法5：自适应wheel算法 ====================
def adaptive_wheel_sieve(n: int) -> List[int]:
    """30030-wheel分段筛 - O(n log log n)时间，O(√n + seg)空间"""
    if n < 2:
        return []
    small_primes = [2, 3, 5, 7, 11, 13]
    if n <= 13:
        return [p for p in small_primes if p <= n]

    sqrt_n = int(math.sqrt(n)) + 1
    max_small = max(small_primes)

    base_primes = []
    is_comp = [False] * (sqrt_n + 1)
    for i in range(2, sqrt_n):
        if not is_comp[i]:
            base_primes.append(i)
            for j in range(i * i, sqrt_n, i):
                is_comp[j] = True

    primes = [p for p in small_primes if p <= n]
    primes.extend([p for p in base_primes if p > max_small and p <= n])

    modulus = 30030

    # 处理sqrt_n到modulus间隙
    gap_start = sqrt_n
    gap_end = min(n, modulus)
    if gap_start <= gap_end:
        for x in range(gap_start, gap_end + 1):
            is_prime = True
            for p in base_primes:
                if p * p > x:
                    break
                if x % p == 0:
                    is_prime = False
                    break
            if is_prime:
                primes.append(x)

    if n <= modulus:
        return sorted(set(primes))

    residues = [i for i in range(1, modulus) if all(i % p != 0 for p in small_primes)]
    seg_size = modulus * 156
    start = modulus + 1

    for low in range(start, n + 1, seg_size):
        high = min(low + seg_size - 1, n)
        length = high - low + 1
        sieve = bytearray(length)

        first_block = (low // modulus) * modulus
        for k in range(first_block, high + 1, modulus):
            for r in residues:
                x = k + r
                if low <= x <= high:
                    sieve[x - low] = 1

        for p in base_primes:
            if p <= max_small or p * p > high:
                continue
            start_mark = ((low + p - 1) // p) * p
            if start_mark < p * p:
                start_mark = p * p
            for j in range(start_mark, high + 1, p):
                if j >= low:
                    sieve[j - low] = 0

        for k in range(first_block, high + 1, modulus):
            for r in residues:
                x = k + r
                if low <= x <= high and sieve[x - low]:
                    primes.append(x)

    return sorted(set(primes))


# ==================== 测试框架 ====================
def benchmark(algorithms: List[Tuple[str, callable]], test_cases: List[int]) -> None:
    """性能对比测试 - 时间、内存、结果一一验证"""
    print("=" * 150)
    print("素数生成算法性能对比（时间 + 内存 + 结果验证）")
    print("=" * 150)

    header = f"{'范围':<12}"
    for name, _ in algorithms:
        header += f"{name:<28}"
    print(header)

    subheader = f"{'':<12}"
    for _ in algorithms:
        subheader += f"{'时间(s)':>8}{'内存(MB)':>10}{'素数个数':>10}"
    print(subheader)
    print("-" * 150)

    for n in test_cases:
        line = f"1-{n:<8}"
        results = []
        reference_count = None

        for name, func in algorithms:
            try:
                tracemalloc.start()
                t0 = time.time()
                primes = func(n)
                elapsed = time.time() - t0
                current, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()

                peak_mb = peak / (1024 * 1024)
                count = len(primes)

                # 结果验证
                if reference_count is None:
                    reference_count = count
                    status = f"{count:,}"
                else:
                    status = f"{count:,}" if count == reference_count else f"{count:,}✗"

                results.append((name, elapsed, peak_mb, count))
                line += f"{elapsed:>8.2f}{peak_mb:>10.2f}{status:>10}"

            except MemoryError:
                tracemalloc.stop()
                line += f"{'[OOM]':>8}{'':>10}{'[ERR]':>10}"
            except Exception as e:
                tracemalloc.stop()
                line += f"{'[ERR]':>8}{'':>10}{'[ERR]':>10}"

        print(line)
        if reference_count:
            print(f"  参考素数个数: {reference_count:,} (π({n}))")

    print()


def theoretical_comparison() -> None:
    """理论复杂度对比"""
    print("=" * 150)
    print("理论复杂度分析")
    print("=" * 150)

    data = [
        ("试除法", "O(n√n / log n)", "O(π(n))", "O(n^1.5)", "最直观，但指数级慢"),
        ("标准筛", "O(n log log n)", "O(n)", "O(n log log n)", "线性对数级，内存线性"),
        (
            "分段筛",
            "O(n log log n)",
            "O(√n)",
            "O(n log log n)",
            "时间同标准筛，内存次线性",
        ),
        (
            "标准Wheel",
            "O(n log log n)",
            "O(n)",
            "O(n log log n)",
            "常数优化，理论标记减少73%",
        ),
        (
            "自适应wheel",
            "O(n log log n)",
            "O(√n + seg)",
            "O(n log log n)",
            "Wheel+分段，常数折中",
        ),
    ]

    print(
        f"{'算法':<12} {'时间复杂度':<25} {'空间复杂度':<15} {'渐进性能':<20} {'理论评价'}"
    )
    print("-" * 150)

    for name, time_c, space_c, asymptotic, note in data:
        print(f"{name:<12} {time_c:<25} {space_c:<15} {asymptotic:<20} {note}")

    print(
        "\n注：所有筛法（除试除法）时间复杂度均为O(n log log n)，差异仅在于常数因子。"
    )
    print()


def feature_comparison() -> None:
    """特性对比矩阵"""
    print("=" * 150)
    print("特性对比矩阵")
    print("=" * 150)

    features = {
        "试除法": ["✓", "✗", "✓", "✓", "教学/验证", "n < 10万"],
        "标准筛": ["✗", "✓★", "✓", "✗", "内存充足", "10万 < n < 1000万"],
        "分段筛": ["✓★", "✓", "✓", "✓★", "通用/嵌入式", "n > 1000万"],
        "标准Wheel": ["✗", "✓★", "✗", "✗", "理论最优", "n < 1000万且内存足"],
        "自适应wheel": ["✓", "✓", "✓", "✓", "折中方案", "n > 1000万通用场景"],
    }

    headers = ["算法", "低内存", "高速度", "易实现", "可扩展", "定位", "适用n范围"]
    col_widths = [12, 8, 8, 8, 8, 15, 20]

    for h, w in zip(headers, col_widths):
        print(f"{h:<{w}}", end="")
    print()
    print("-" * 150)

    for algo, vals in features.items():
        print(f"{algo:<12}", end="")
        for v, w in zip(vals, col_widths[1:]):
            print(f"{v:<{w}}", end="")
        print()

    print("注：★表示该维度最优或接近最优。")
    print()


if __name__ == "__main__":
    theoretical_comparison()
    feature_comparison()

    algorithms = [
        ("试除法", trial_division),
        ("标准筛", standard_sieve),
        ("分段筛", segmented_sieve),
        ("标准Wheel", standard_wheel_sieve),
        ("自适应wheel", adaptive_wheel_sieve),
    ]

    test_cases = [100_000, 500_000, 1_000_000, 10_000_000]
    benchmark(algorithms, test_cases)
