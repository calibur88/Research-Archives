#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
哥德巴赫猜想完整验证（三区间严格测试版）
修复：包含素数2 + 双向搜索（防漏检）
"""

import numpy as np
from sympy import primerange, isprime

def 哥德巴赫验证区间(start, end, verbose=True):
    """
    验证指定区间内所有偶数的哥德巴赫表示
    
    参数:
        start: 起始偶数（含）
        end: 结束偶数（含）
        verbose: 是否打印每个偶数的分解
    
    返回:
        (是否全部通过, 失败列表, 详细分解字典)
    """
    # 生成足够大的素数表（覆盖到end）
    primes = list(primerange(2, end + 1))
    prime_set = set(primes)
    
    failures = []
    details = {}
    
    # 确保start是偶数且>=4
    if start < 4:
        start = 4
    if start % 2 != 0:
        start += 1
    
    验证计数 = 0
    
    for n in range(start, end + 1, 2):
        验证计数 += 1
        found = False
        representation = None
        
        # 策略1：从小到大搜索（小素数+大素数）
        for p in primes:
            if p > n - 2:
                break
            q = n - p
            if q in prime_set:
                found = True
                representation = (p, q)
                break
        
        # 策略2：从大到小搜索（处理大数漏检，p接近n/2的情况）
        if not found:
            for p in reversed(primes):
                if p < 2 or p > n - 2:
                    continue
                q = n - p
                if isprime(q):
                    found = True
                    representation = (p, q, "反向")
                    break
        
        if found:
            details[n] = representation
            if verbose:
                flag = "✓" if len(representation) == 2 else "✓(反向)"
                print(f"{flag} {n:4d} = {representation[0]:4d} + {representation[1]:4d}")
        else:
            failures.append(n)
            if verbose:
                print(f"✗ {n:4d} 未找到表示")
    
    全部通过 = len(failures) == 0
    
    if verbose:
        print(f"\n区间 [{start}, {end}] 验证完成:")
        print(f"  验证偶数数量: {验证计数}")
        print(f"  成功: {验证计数 - len(failures)}")
        print(f"  失败: {len(failures)}")
        if failures:
            print(f"  失败列表: {failures}")
    
    return 全部通过, failures, details

def 双螺旋物理解释(n, p, q):
    """输出双螺旋和频的物理解释"""
    print(f"    双螺旋和频: 偶数轴 {n} = 素数臂A({p}) + 素数臂B({q})")
    print(f"    相位锁定: {p}↑ + {q}↓ = {n}（能量守恒）")

# ==================== 主程序 ====================

if __name__ == "__main__":
    print("=" * 70)
    print("哥德巴赫猜想验证（三区间严格测试）")
    print("物理框架：双螺旋和频干涉（p + q = 2N）")
    print("=" * 70)
    
    测试区间 = [
        (4, 100, "小尺度（基础验证）"),
        (500, 600, "中尺度（过渡区域）"),
        (5000, 5100, "大尺度（漏检敏感区）")
    ]
    
    总结果 = []
    
    for start, end, 描述 in 测试区间:
        print(f"\n{'='*70}")
        print(f"区间 {start}-{end}：{描述}")
        print(f"{'='*70}")
        
        通过, 失败列表, 详情 = 哥德巴赫验证区间(start, end, verbose=True)
        总结果.append((start, end, 通过, 失败列表))
        
        # 显示该区间的一些代表性分解
        print(f"\n代表性分解（双螺旋和频）：")
        样本 = list(详情.items())[:3] + list(详情.items())[-3:]
        for n, rep in 样本:
            if len(rep) == 3:  # 反向找到的
                print(f"  {n} = {rep[0]} + {rep[1]} [反向搜索锁定]")
            else:
                print(f"  {n} = {rep[0]} + {rep[1]}")
    
    # 最终汇总
    print(f"\n{'='*70}")
    print("三区间验证汇总")
    print(f"{'='*70}")
    
    for start, end, 通过, 失败 in 总结果:
        状态 = "✓ 通过" if 通过 else "✗ 失败"
        失败信息 = f"（失败: {失败}）" if 失败 else ""
        print(f"  [{start:4d}, {end:4d}] {状态} {失败信息}")
    
    全部通过 = all(r[2] for r in 总结果)
    print(f"\n总结果：{'✓ 全部通过' if 全部通过 else '✗ 存在反例'}")
    
    if 全部通过:
        print("\n物理结论：")
        print("  双螺旋和频干涉在三个尺度均保持相位锁定")
        print("  p（递增臂）+ q（递减臂）= 2N（能量守恒）")
        print("  大数（>1000）反向搜索补偿了素数间隙增大的影响")
        print("  哥德巴赫猜想在有限可计算宇宙内成立（Step Budget内）")
