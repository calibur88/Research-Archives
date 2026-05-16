#!/usr/bin/env python3
"""
DyVeRT  — 构图修复 + 骨骼动画 GIF + 评估指标
修复:
  1. 相机更近更高, 构图居中
  2. Y轴反转修复 (y=1 对应图像顶部)
  3. 地面提亮, 环境光增强, 人物可见
  4. 骨骼动画后自动重建人物 BVH
  5. 8帧 GIF 输出
  6. GT 16x 累积 + SSIM/MSE/梯度评估
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, uniform_filter
from PIL import Image
import time, warnings
warnings.filterwarnings("ignore")

# ============================================================
# 0. SSIM
# ============================================================
def ssim_numpy(img1, img2, K1=0.01, K2=0.03, data_range=1.0):
    C1 = (K1 * data_range) ** 2
    C2 = (K2 * data_range) ** 2
    mu1 = uniform_filter(img1, 7, mode="reflect")
    mu2 = uniform_filter(img2, 7, mode="reflect")
    mu1_sq, mu2_sq, mu1_mu2 = mu1**2, mu2**2, mu1*mu2
    s1_sq = uniform_filter(img1**2, 7, mode="reflect") - mu1_sq
    s2_sq = uniform_filter(img2**2, 7, mode="reflect") - mu2_sq
    s12 = uniform_filter(img1*img2, 7, mode="reflect") - mu1_mu2
    ssim_map = ((2*mu1_mu2+C1)*(2*s12+C2)) / ((mu1_sq+mu2_sq+C1)*(s1_sq+s2_sq+C2))
    return np.mean(ssim_map)

# ============================================================
# 1. AABB
# ============================================================
class AABB:
    __slots__ = ('min', 'max', 'center')
    def __init__(self, min_pt, max_pt):
        self.min = np.asarray(min_pt, dtype=np.float32)
        self.max = np.asarray(max_pt, dtype=np.float32)
        self.center = (self.min + self.max) * 0.5

    def intersect_batch(self, o, d, t_max):
        inv_d = 1.0 / (d + 1e-9)
        t1 = (self.min - o) * inv_d
        t2 = (self.max - o) * inv_d
        t_near = np.maximum(np.maximum(np.minimum(t1[:,0], t2[:,0]),
                                        np.minimum(t1[:,1], t2[:,1])),
                                        np.minimum(t1[:,2], t2[:,2]))
        t_far  = np.minimum(np.minimum(np.maximum(t1[:,0], t2[:,0]),
                                        np.maximum(t1[:,1], t2[:,1])),
                                        np.maximum(t1[:,2], t2[:,2]))
        hit = (t_far > t_near) & (t_far > 0.001)
        if isinstance(t_max, np.ndarray):
            hit &= t_near < t_max
        t_hit = np.where(t_near > 0.001, t_near, np.where(t_far > 0.001, t_far, np.inf))
        t_hit = t_hit.astype(np.float32)
        t_hit[~hit] = np.inf
        return t_hit, hit

# ============================================================
# 2. 图元
# ============================================================
class Prim:
    pass

class BoxPrim(Prim):
    __slots__ = ('aabb','color','rough','metal')
    def __init__(self, bmin, bmax, color, rough=0.9, metal=0.0):
        self.aabb = AABB(np.asarray(bmin, np.float32), np.asarray(bmax, np.float32))
        self.color = np.asarray(color, np.float32)
        self.rough = float(rough)
        self.metal = float(metal)

    def intersect_batch(self, o, d, t_max):
        N = o.shape[0]
        t, hit = self.aabb.intersect_batch(o, d, t_max)
        if not np.any(hit):
            return t, hit, None, None, None, None, None
        pos = o[hit] + d[hit] * t[hit][:,None]
        n = np.zeros((hit.sum(), 3), dtype=np.float32)
        eps = 0.001
        p = pos
        bmin, bmax = self.aabb.min, self.aabb.max
        near_x  = np.abs(p[:,0] - bmin[0]) < eps
        near_x2 = np.abs(p[:,0] - bmax[0]) < eps
        near_y  = np.abs(p[:,1] - bmin[1]) < eps
        near_y2 = np.abs(p[:,1] - bmax[1]) < eps
        near_z  = np.abs(p[:,2] - bmin[2]) < eps
        near_z2 = np.abs(p[:,2] - bmax[2]) < eps
        n[near_x]  = [-1,0,0]; n[near_x2] = [1,0,0]
        n[near_y]  = [0,-1,0]; n[near_y2] = [0,1,0]
        n[near_z]  = [0,0,-1]; n[near_z2] = [0,0,1]
        mask = ~(near_x | near_x2 | near_y | near_y2 | near_z | near_z2)
        n[mask] = [0,0,1]
        n = n / (np.linalg.norm(n, axis=1, keepdims=True) + 1e-9)
        alb = np.broadcast_to(self.color, (hit.sum(), 3))
        rgh = np.full(hit.sum(), self.rough, dtype=np.float32)
        met = np.full(hit.sum(), self.metal, dtype=np.float32)
        return t, hit, n, pos, alb, rgh, met

class SpherePrim(Prim):
    __slots__ = ('center','radius','color','rough','metal','aabb')
    def __init__(self, center, radius, color, rough=0.5, metal=0.0):
        self.center = np.asarray(center, np.float32)
        self.radius = float(radius)
        self.color  = np.asarray(color, np.float32)
        self.rough  = float(rough)
        self.metal  = float(metal)
        r = self.radius
        self.aabb = AABB(self.center - r, self.center + r)

    def intersect_batch(self, o, d, t_max):
        N = o.shape[0]
        oc = o - self.center
        a = np.sum(d*d, axis=1)
        b_half = 2.0*np.sum(oc*d, axis=1)
        c = np.sum(oc*oc, axis=1) - self.radius**2
        disc = b_half**2 - 4*a*c
        valid = disc >= 0
        if not np.any(valid):
            return np.full(N, np.inf, np.float32), np.zeros(N, dtype=bool), None,None,None,None,None
        sd = np.sqrt(np.maximum(disc[valid], 0))
        av = a[valid]; bv = b_half[valid]
        t1 = (-bv - sd) / (2*av)
        t2 = (-bv + sd) / (2*av)
        t_hit = np.where((t1>0.001)&(t1<t2), t1, np.where(t2>0.001, t2, np.inf))
        t_full = np.full(N, np.inf, np.float32)
        t_full[valid] = t_hit
        hit = (t_full < t_max) & (t_full > 0.001)
        t_full[~hit] = np.inf
        if not np.any(hit):
            return t_full, hit, None,None,None,None,None
        pos = o[hit] + d[hit] * t_full[hit][:,None]
        n = pos - self.center
        n = n / (np.linalg.norm(n, axis=1, keepdims=True) + 1e-9)
        alb = np.broadcast_to(self.color, (hit.sum(), 3))
        rgh = np.full(hit.sum(), self.rough, dtype=np.float32)
        met = np.full(hit.sum(), self.metal, dtype=np.float32)
        return t_full, hit, n, pos, alb, rgh, met

class CapsulePrim(Prim):
    __slots__ = ('a','b','radius','color','rough','metal','axis','axis_len2','aabb')
    def __init__(self, a, b, radius, color, rough=0.7, metal=0.0):
        self.a = np.asarray(a, np.float32)
        self.b = np.asarray(b, np.float32)
        self.radius = float(radius)
        self.color  = np.asarray(color, np.float32)
        self.rough  = float(rough)
        self.metal  = float(metal)
        self.update_pose(a, b)

    def update_pose(self, a, b):
        self.a = np.asarray(a, np.float32)
        self.b = np.asarray(b, np.float32)
        self.axis = self.b - self.a
        self.axis_len2 = float(np.dot(self.axis, self.axis))
        r = self.radius
        self.aabb = AABB(np.minimum(self.a, self.b) - r,
                         np.maximum(self.a, self.b) + r)

    def intersect_batch(self, o, d, t_max):
        N = o.shape[0]
        ao = o - self.a
        t_proj = np.sum(ao * self.axis, axis=1) / (self.axis_len2 + 1e-9)
        t_proj = np.clip(t_proj, 0.0, 1.0)
        closest = self.a + t_proj[:,None] * self.axis
        oc = o - closest
        a = np.sum(d*d, axis=1)
        b_half = 2.0*np.sum(oc*d, axis=1)
        c = np.sum(oc*oc, axis=1) - self.radius**2
        disc = b_half**2 - 4*a*c
        valid = disc >= 0
        if not np.any(valid):
            return np.full(N, np.inf, np.float32), np.zeros(N, dtype=bool), None,None,None,None,None
        sd = np.sqrt(np.maximum(disc[valid], 0))
        av = a[valid]; bv = b_half[valid]
        t1 = (-bv - sd) / (2*av)
        t2 = (-bv + sd) / (2*av)
        t_hit = np.where((t1>0.001)&(t1<t2), t1, np.where(t2>0.001, t2, np.inf))
        t_full = np.full(N, np.inf, np.float32)
        t_full[valid] = t_hit
        hit = (t_full < t_max) & (t_full > 0.001)
        t_full[~hit] = np.inf
        if not np.any(hit):
            return t_full, hit, None,None,None,None,None
        pos = o[hit] + d[hit] * t_full[hit][:,None]
        n = pos - closest[hit]
        n = n / (np.linalg.norm(n, axis=1, keepdims=True) + 1e-9)
        alb = np.broadcast_to(self.color, (hit.sum(), 3))
        rgh = np.full(hit.sum(), self.rough, dtype=np.float32)
        met = np.full(hit.sum(), self.metal, dtype=np.float32)
        return t_full, hit, n, pos, alb, rgh, met

class PlanePrim(Prim):
    __slots__ = ('point','normal','color','rough','metal','aabb')
    def __init__(self, point, normal, color, rough=0.95):
        self.point  = np.asarray(point, np.float32)
        self.normal = np.asarray(normal, np.float32)
        self.normal = self.normal / (np.linalg.norm(self.normal) + 1e-9)
        self.color  = np.asarray(color, np.float32)
        self.rough  = float(rough)
        self.metal  = 0.0
        big = 1000.0
        y = self.point[1]
        self.aabb = AABB(np.array([-big, y-0.01, -big], np.float32),
                         np.array([ big, y+0.01,  big], np.float32))

    def intersect_batch(self, o, d, t_max):
        N = o.shape[0]
        denom = np.sum(d * self.normal, axis=1)
        valid = np.abs(denom) > 1e-6
        if not np.any(valid):
            return np.full(N, np.inf, np.float32), np.zeros(N, dtype=bool), None,None,None,None,None
        t = np.sum((self.point - o) * self.normal, axis=1) / denom
        t_full = np.full(N, np.inf, np.float32)
        t_full[valid] = t[valid]
        hit = (t_full > 0.001) & (t_full < t_max)
        t_full[~hit] = np.inf
        if not np.any(hit):
            return t_full, hit, None,None,None,None,None
        pos = o[hit] + d[hit] * t_full[hit][:,None]
        n = np.broadcast_to(self.normal, (hit.sum(), 3))
        alb = np.broadcast_to(self.color, (hit.sum(), 3))
        rgh = np.full(hit.sum(), self.rough, dtype=np.float32)
        met = np.zeros(hit.sum(), dtype=np.float32)
        return t_full, hit, n, pos, alb, rgh, met

# ============================================================
# 3. 局部 BVH
# ============================================================
class BVHNode:
    __slots__ = ('aabb','left','right','prims','is_leaf')
    def __init__(self, aabb, left=None, right=None, prims=None):
        self.aabb = aabb
        self.left = left
        self.right = right
        self.prims = prims
        self.is_leaf = prims is not None

class LocalBVH:
    def __init__(self, prims):
        self.prims = prims
        self.root = self._build(prims)

    def _union_aabb(self, prims):
        mins = np.array([p.aabb.min for p in prims], dtype=np.float32)
        maxs = np.array([p.aabb.max for p in prims], dtype=np.float32)
        return AABB(mins.min(axis=0), maxs.max(axis=0))

    def _build(self, prims):
        if len(prims) <= 3:
            return BVHNode(self._union_aabb(prims), prims=prims)
        centers = np.array([p.aabb.center for p in prims], dtype=np.float32)
        aabb = self._union_aabb(prims)
        axis = np.argmax(aabb.max - aabb.min)
        sorted_idx = np.argsort(centers[:, axis])
        mid = len(prims) // 2
        left  = self._build([prims[i] for i in sorted_idx[:mid]])
        right = self._build([prims[i] for i in sorted_idx[mid:]])
        return BVHNode(aabb, left=left, right=right)

    def intersect_batch(self, o, d, t_max):
        N = o.shape[0]
        best_t = np.full(N, np.inf, dtype=np.float32)
        best_n = np.zeros((N, 3), dtype=np.float32)
        best_p = np.zeros((N, 3), dtype=np.float32)
        best_a = np.zeros((N, 3), dtype=np.float32)
        best_r = np.zeros(N, dtype=np.float32)
        best_m = np.zeros(N, dtype=np.float32)

        stack = [self.root]
        while stack:
            node = stack.pop()
            t_aabb, hit = node.aabb.intersect_batch(o, d, t_max=best_t)
            if not np.any(hit): continue
            if node.is_leaf:
                idx = np.where(hit)[0]
                o_sub = o[idx]; d_sub = d[idx]
                t_sub_best = best_t[idx]
                for prim in node.prims:
                    t_sub, h_sub, n_sub, p_sub, a_sub, r_sub, m_sub = prim.intersect_batch(o_sub, d_sub, t_max=t_sub_best)
                    if not np.any(h_sub): continue
                    sub_local = np.where(h_sub)[0]
                    t_valid = t_sub[h_sub]
                    best_valid = t_sub_best[h_sub]
                    update_mask = t_valid < best_valid
                    if not np.any(update_mask): continue
                    u = idx[sub_local[update_mask]]
                    best_t[u] = t_valid[update_mask]
                    best_n[u] = n_sub[update_mask]
                    best_p[u] = p_sub[update_mask]
                    best_a[u] = a_sub[update_mask]
                    best_r[u] = r_sub[update_mask]
                    best_m[u] = m_sub[update_mask]
                    t_sub_best[sub_local[update_mask]] = t_valid[update_mask]
            else:
                stack.append(node.right)
                stack.append(node.left)
        hit = best_t < np.inf
        return best_t, hit, best_n, best_p, best_a, best_r, best_m

# ============================================================
# 4. TLAS
# ============================================================
class TLASNode:
    __slots__ = ('aabb','left','right','objs','is_leaf')
    def __init__(self, aabb, left=None, right=None, objs=None):
        self.aabb = aabb
        self.left = left
        self.right = right
        self.objs = objs
        self.is_leaf = objs is not None

class TLAS:
    def __init__(self, objects):
        self.objects = objects
        self.root = self._build(objects)

    def _build(self, objs):
        if len(objs) <= 2:
            mins = np.array([o.world_aabb.min for o in objs], dtype=np.float32)
            maxs = np.array([o.world_aabb.max for o in objs], dtype=np.float32)
            return TLASNode(AABB(mins.min(axis=0), maxs.max(axis=0)), objs=objs)
        centers = np.array([o.world_aabb.center for o in objs], dtype=np.float32)
        mins = np.array([o.world_aabb.min for o in objs], dtype=np.float32)
        maxs = np.array([o.world_aabb.max for o in objs], dtype=np.float32)
        aabb = AABB(mins.min(axis=0), maxs.max(axis=0))
        axis = np.argmax(aabb.max - aabb.min)
        sorted_idx = np.argsort(centers[:, axis])
        mid = len(objs) // 2
        left  = self._build([objs[i] for i in sorted_idx[:mid]])
        right = self._build([objs[i] for i in sorted_idx[mid:]])
        return TLASNode(aabb, left=left, right=right)

    def intersect_batch(self, o, d, t_max=None):
        N = o.shape[0]
        if t_max is None:
            t_max = np.full(N, np.inf, dtype=np.float32)
        best_t = t_max.copy()
        best_obj = np.full(N, -1, dtype=np.int32)
        best_n = np.zeros((N, 3), dtype=np.float32)
        best_p = np.zeros((N, 3), dtype=np.float32)
        best_a = np.zeros((N, 3), dtype=np.float32)
        best_r = np.zeros(N, dtype=np.float32)
        best_m = np.zeros(N, dtype=np.float32)

        stack = [self.root]
        while stack:
            node = stack.pop()
            t_aabb, hit = node.aabb.intersect_batch(o, d, t_max=best_t)
            if not np.any(hit): continue
            if node.is_leaf:
                idx = np.where(hit)[0]
                o_sub = o[idx]; d_sub = d[idx]
                t_sub_best = best_t[idx]
                for obj in node.objs:
                    t_sub, h_sub, n_sub, p_sub, a_sub, r_sub, m_sub = obj.local_bvh.intersect_batch(o_sub, d_sub, t_max=t_sub_best)
                    if not np.any(h_sub): continue
                    update = h_sub & (t_sub < t_sub_best)
                    if not np.any(update): continue
                    u = idx[update]
                    best_t[u] = t_sub[update]
                    best_obj[u] = obj.id
                    best_n[u] = n_sub[update]
                    best_p[u] = p_sub[update]
                    best_a[u] = a_sub[update]
                    best_r[u] = r_sub[update]
                    best_m[u] = m_sub[update]
                    t_sub_best[update] = t_sub[update]
            else:
                stack.append(node.right)
                stack.append(node.left)
        hit = best_t < np.inf
        return best_t, best_obj, best_n, best_p, best_a, best_r, best_m

# ============================================================
# 5. 场景对象 (含重建)
# ============================================================
class SceneObject:
    __slots__ = ('id','prims','local_bvh','world_aabb','has_plane')
    def __init__(self, obj_id, prims):
        self.id = obj_id
        self.prims = prims
        self.local_bvh = LocalBVH(prims)
        self.has_plane = any(isinstance(p, PlanePrim) for p in prims)
        self.update_world_aabb()

    def update_world_aabb(self):
        mins = np.array([p.aabb.min for p in self.prims], dtype=np.float32)
        maxs = np.array([p.aabb.max for p in self.prims], dtype=np.float32)
        self.world_aabb = AABB(mins.min(axis=0), maxs.max(axis=0))

    def rebuild(self):
        """动画后重建局部 BVH 与世界 AABB"""
        self.local_bvh = LocalBVH(self.prims)
        self.update_world_aabb()

# ============================================================
# 6. 骨骼系统
# ============================================================
class Bone:
    __slots__ = ('name','parent','local_pos','length','radius','color','rough',
                 'angle_x','angle_z','world_mat')
    def __init__(self, name, parent, local_pos, length, radius, color, rough=0.7):
        self.name = name
        self.parent = parent
        self.local_pos = np.asarray(local_pos, np.float32)
        self.length = float(length)
        self.radius = float(radius)
        self.color = np.asarray(color, np.float32)
        self.rough = float(rough)
        self.angle_x = 0.0
        self.angle_z = 0.0
        self.world_mat = np.eye(4, dtype=np.float32)

    def compute_matrix(self, parent_mat):
        cx, cz = np.cos(self.angle_x), np.cos(self.angle_z)
        sx, sz = np.sin(self.angle_x), np.sin(self.angle_z)
        Rx = np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]], dtype=np.float32)
        Rz = np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]], dtype=np.float32)
        R = Rz @ Rx
        M = np.eye(4, dtype=np.float32)
        M[:3, :3] = R
        M[:3, 3] = self.local_pos
        self.world_mat = parent_mat @ M
        return self.world_mat

    def endpoints(self):
        start = self.world_mat[:3, 3].copy()
        local_end = np.array([0, 0, self.length], dtype=np.float32)
        end = start + self.world_mat[:3, :3] @ local_end
        return start, end

class Skeleton:
    def __init__(self, root_pos):
        self.root_pos = np.asarray(root_pos, np.float32)
        self.bones = {}
        self.prims = []
        self._build()

    def _build(self):
        def add(name, parent, lp, ln, rad, col, rough=0.7):
            b = Bone(name, parent, lp, ln, rad, col, rough)
            self.bones[name] = b
            return b

        hips   = add("hips",   None, [0,0,0],       0.22, 0.24, [0.25,0.30,0.55])
        spine  = add("spine",  hips,  [0,0,0.22],   0.42, 0.20, [0.25,0.30,0.55])
        head   = add("head",   spine, [0,0,0.42],   0.22, 0.18, [0.90,0.75,0.65])
        l_arm  = add("l_arm",  spine, [0.38,0,0.35],0.32, 0.10, [0.90,0.75,0.65])
        l_fore = add("l_fore", l_arm, [0,0,0.32],   0.28, 0.09, [0.90,0.75,0.65])
        r_arm  = add("r_arm",  spine, [-0.38,0,0.35],0.32,0.10, [0.90,0.75,0.65])
        r_fore = add("r_fore", r_arm, [0,0,0.32],   0.28, 0.09, [0.90,0.75,0.65])
        l_thigh= add("l_thigh",hips, [0.16,0,0],   0.38, 0.13, [0.15,0.15,0.28])
        l_shin = add("l_shin",l_thigh,[0,0,-0.38], 0.38, 0.11, [0.15,0.15,0.28])
        r_thigh= add("r_thigh",hips, [-0.16,0,0],  0.38, 0.13, [0.15,0.15,0.28])
        r_shin = add("r_shin",r_thigh,[0,0,-0.38], 0.38, 0.11, [0.15,0.15,0.28])

        for b in self.bones.values():
            prim = CapsulePrim([0,0,0],[0,0,1], b.radius, b.color, b.rough, 0.0)
            prim.bone_name = b.name
            self.prims.append(prim)

    def animate(self, t):
        swing  = np.sin(t * 3.0) * 0.45
        swing2 = np.sin(t * 3.0 + np.pi) * 0.45

        self.bones["l_thigh"].angle_x = swing
        self.bones["l_shin"].angle_x  = -abs(swing) * 0.9
        self.bones["r_thigh"].angle_x = swing2
        self.bones["r_shin"].angle_x  = -abs(swing2) * 0.9

        self.bones["l_arm"].angle_x  = swing2 * 0.55
        self.bones["l_fore"].angle_x = -abs(swing2) * 0.35
        self.bones["r_arm"].angle_x  = swing * 0.55
        self.bones["r_fore"].angle_x = -abs(swing) * 0.35

        self.bones["spine"].angle_z = np.sin(t * 1.5) * 0.06
        self.bones["hips"].local_pos[2] = self.root_pos[2] + abs(np.sin(t * 3.0)) * 0.06

        def walk(bone, parent_mat):
            if bone.parent is None:
                bone.local_pos[:2] = self.root_pos[:2]
            bone.compute_matrix(parent_mat)
            for child in self.bones.values():
                if child.parent is bone:
                    walk(child, bone.world_mat)
        walk(self.bones["hips"], np.eye(4, dtype=np.float32))

        for prim in self.prims:
            a, b = self.bones[prim.bone_name].endpoints()
            prim.update_pose(a, b)

# ============================================================
# 7. 场景构建
# ============================================================
def build_scene(t_anim):
    objects = []
    oid = 0

    # 地面 (提亮)
    objects.append(SceneObject(oid, [PlanePrim([0,-0.5,0], [0,1,0], [0.75,0.70,0.55], 0.95)]))
    oid += 1

    # 建筑群
    buildings = [
        ([-3.0, -3.0, -0.5], [ 0.5, 0.8, 2.8], [0.82, 0.72, 0.52]),
        ([ 1.2, -2.2, -0.5], [ 3.0, 1.2, 3.2], [0.42, 0.52, 0.65]),
        ([-2.0,  2.0, -0.5], [ 1.2, 0.9, 2.0], [0.72, 0.32, 0.22]),
        ([ 3.0,  1.0, -0.5], [ 1.0, 1.8, 2.5], [0.65, 0.60, 0.55]),
        ([ 0.0,  4.0, -0.5], [ 2.2, 1.0, 2.2], [0.55, 0.70, 0.45]),
        ([-4.5, -1.0, -0.5], [ 0.8, 2.0, 1.8], [0.60, 0.58, 0.50]),
    ]
    for bmin, bmax, col in buildings:
        objects.append(SceneObject(oid, [BoxPrim(np.array(bmin), np.array(bmax), np.array(col), 0.85, 0.0)]))
        oid += 1

    # 球体
    spheres = [
        ([0.0, 0.0, 2.8], 0.75, [0.92, 0.08, 0.08], 0.75, 0.0),
        ([2.8, 1.8, 1.6], 0.50, [0.08, 0.82, 0.12], 0.65, 0.0),
        ([-1.8, -1.8, 2.2], 0.60, [0.90, 0.90, 0.95], 0.03, 1.0),
        ([3.8, -1.2, 1.2], 0.40, [0.92, 0.72, 0.08], 0.08, 1.0),
        ([-3.8, 3.2, 1.8], 0.55, [0.85, 0.50, 0.90], 0.20, 0.5),
        ([1.5, -3.5, 0.9], 0.45, [0.20, 0.60, 0.90], 0.10, 0.3),
    ]
    for c, r, col, rough, metal in spheres:
        objects.append(SceneObject(oid, [SpherePrim(np.array(c), r, np.array(col), rough, metal)]))
        oid += 1

    # 树木
    trees = [
        ([-4.2, -2.2, -0.5], 0.14, 2.6, [0.22, 0.16, 0.10], [0.18, 0.42, 0.14]),
        ([ 4.8, -2.8, -0.5], 0.12, 2.1, [0.28, 0.20, 0.14], [0.14, 0.38, 0.12]),
        ([ 2.8,  3.8, -0.5], 0.16, 3.2, [0.20, 0.16, 0.12], [0.20, 0.46, 0.16]),
        ([-1.2, -4.2, -0.5], 0.18, 2.9, [0.26, 0.18, 0.14], [0.22, 0.44, 0.15]),
        ([ 5.2,  2.2, -0.5], 0.13, 2.4, [0.24, 0.17, 0.11], [0.16, 0.40, 0.13]),
        ([-3.0,  4.5, -0.5], 0.15, 2.7, [0.21, 0.15, 0.09], [0.19, 0.43, 0.14]),
        ([ 0.5, -5.0, -0.5], 0.14, 2.3, [0.25, 0.19, 0.13], [0.17, 0.39, 0.12]),
        ([-5.0,  0.5, -0.5], 0.16, 2.8, [0.23, 0.17, 0.12], [0.18, 0.41, 0.13]),
    ]
    for center, r, h, trunk_col, leaf_col in trees:
        c = np.array(center)
        trunk = CapsulePrim(c, c + np.array([0, h, 0]), r, np.array(trunk_col), 0.95, 0.0)
        crown = SpherePrim(c + np.array([0, h + 0.4, 0]), 0.85, np.array(leaf_col), 0.9, 0.0)
        objects.append(SceneObject(oid, [trunk, crown]))
        oid += 1

    # 骨骼人物 (放在场景中心前方)
    skel = Skeleton(np.array([0.5, 0.5, -0.5]))
    skel.animate(t_anim)
    person_obj = SceneObject(oid, skel.prims)
    objects.append(person_obj)
    oid += 1

    return objects, skel, person_obj

# ============================================================
# 8. 自动相机 (更近更高)
# ============================================================
def auto_camera(objects):
    pts = []
    for obj in objects:
        if obj.has_plane:
            continue
        pts.append(obj.world_aabb.min)
        pts.append(obj.world_aabb.max)
    if len(pts) == 0:
        return np.array([5.0, -8.0, 5.0], np.float32), np.array([0.0, 0.0, 1.0], np.float32)
    V = np.array(pts, dtype=np.float32)
    vmin = V.min(axis=0)
    vmax = V.max(axis=0)
    o = (vmin + vmax) * 0.5
    diag = np.linalg.norm(vmax - vmin)
    # 更高更侧的角度
    d_vec = np.array([2.5, -2.5, 2.0], dtype=np.float32)
    d_vec = d_vec / (np.linalg.norm(d_vec) + 1e-9)
    cam = o + diag * 0.55 * d_vec
    lookat = o + np.array([0, 0, 0.3], dtype=np.float32)
    print("  Scene bbox: min=%s max=%s diag=%.2f" % (str(vmin.round(2)), str(vmax.round(2)), diag))
    return cam, lookat

# ============================================================
# 9. 渲染 (Y轴修复: y=1 在图像顶部)
# ============================================================
lights = [
    {"pos": np.array([6.0, 10.0, 8.0],  dtype=np.float32), "intensity": 400.0, "color": np.array([1.00, 0.95, 0.80], np.float32)},
    {"pos": np.array([-2.0, 4.0, 6.0], dtype=np.float32), "intensity": 45.0,  "color": np.array([0.40, 0.60, 1.00], np.float32)},
    {"pos": np.array([0.0, -8.0, 4.0], dtype=np.float32), "intensity": 40.0,  "color": np.array([1.00, 0.90, 0.55], np.float32)},
]

def render_frame(cam_pos, lookat, objects, W, H, noise_sigma=0.05):
    N = H * W
    vd = lookat - cam_pos
    vd = vd / (np.linalg.norm(vd) + 1e-9)
    rt = np.cross(vd, np.array([0.0, 0.0, 1.0]))
    rt = rt / (np.linalg.norm(rt) + 1e-9)
    up = np.cross(rt, vd)

    fov = np.pi / 2.5
    asp = W / H
    tf = np.tan(fov / 2)
    x = np.linspace(-1, 1, W)
    # 修复: y=1 对应图像顶部 (行0), y=-1 对应底部 (行H-1)
    y = np.linspace(1, -1, H)
    X, Y = np.meshgrid(x, y)

    rd = (vd[:,None,None] + X[None,:,:]*rt[:,None,None]*tf*asp + Y[None,:,:]*up[:,None,None]*tf)
    rd = rd / np.linalg.norm(rd, axis=0, keepdims=True)

    o = np.broadcast_to(cam_pos[:,None,None], (3,H,W)).reshape(3,N).T.astype(np.float32)
    d = rd.reshape(3,N).T.astype(np.float32)

    t0 = time.time()
    tlas = TLAS(objects)
    tlas_t = time.time() - t0

    t0 = time.time()
    t, obj_id, normals, pos, albedo, rough, metal = tlas.intersect_batch(o, d)
    trace_t = time.time() - t0

    hit = t < np.inf
    color = np.zeros((N, 3), dtype=np.float32)

    for lt in lights:
        lp = lt["pos"]; li = lt["intensity"]; lc = lt["color"]
        lv = lp - pos
        ld = np.linalg.norm(lv, axis=1)
        ldr = lv / (ld[:,None] + 1e-9)

        shadow_t, _, _, _, _, _, _ = tlas.intersect_batch(pos + normals*0.002, ldr, t_max=ld)
        in_shadow = shadow_t < ld

        nd = np.clip(np.sum(normals * ldr, axis=1), 0, 1)
        nd[in_shadow] = 0.03
        att = li / (ld**2 + 1.0)
        diffuse = albedo * (1.0 - metal[:,None]) * nd[:,None] / np.pi

        view = -d / (np.linalg.norm(d, axis=1, keepdims=True) + 1e-9)
        half = ldr + view
        half = half / (np.linalg.norm(half, axis=1, keepdims=True) + 1e-9)
        ndh = np.clip(np.sum(normals * half, axis=1), 0, 1)
        spec_p = np.maximum(2.0/(rough**2 + 0.001) - 2.0, 1.0)
        spec = np.power(ndh, spec_p) * (1.0 - rough) * metal
        color += (diffuse + spec[:,None]) * lc * att[:,None]

    color += albedo * 0.05  # 增强环境光
    sky_color = np.array([0.5, 0.6, 0.8], dtype=np.float32)
    color[~hit] = sky_color * (1.0 + 0.2*d[~hit, 2][:,None])
    albedo[~hit] = sky_color
    normals[~hit] = [0,0,1]

    ns = np.random.normal(0, noise_sigma, color.shape) * np.sqrt(np.abs(color) + 0.05)
    hdr = np.clip(color + ns, 0, None)

    return hdr.reshape(H,W,3), albedo.reshape(H,W,3), normals.reshape(H,W,3), tlas_t, trace_t

# ============================================================
# 10. 后处理
# ============================================================
def tm_reinhard(hdr, exposure=1.2, white=5.0):
    s = hdr * exposure
    mapped = s * (1.0 + s / (white**2 + 1e-6)) / (1.0 + s + 1e-6)
    ldr = np.where(mapped <= 0.0031308, 12.92*mapped, 1.055*np.power(mapped, 1/2.4) - 0.055)
    return np.clip(ldr, 0, 1)

def bloom(img, th=0.75, st=0.05):
    b = np.maximum(img - th, 0)
    out = np.zeros_like(img)
    for k in range(1,4):
        sig = 1.5 * (2**(k-1))
        w = 1.0 / (2**k)
        for c in range(3):
            out[:,:,c] += w * gaussian_filter(b[:,:,c], sigma=sig)
    return np.clip(img + st*out, 0, 1)

def pipe_unity(hdr, alb, nrm, bloom_th=0.75, bloom_st=0.05):
    H,W,C = hdr.shape
    ss, sr, sn, rad = 2.0, 0.25, 0.5, 3
    bil = np.zeros_like(hdr)
    for c in range(C):
        ph = np.pad(hdr[:,:,c], rad, mode="reflect")
        pa = np.pad(alb[:,:,c], rad, mode="reflect")
        pn = np.pad(nrm[:,:,c], rad, mode="reflect")
        out = np.zeros((H,W)); norm = np.zeros((H,W))
        for dy in range(-rad, rad+1):
            for dx in range(-rad, rad+1):
                nh = ph[rad+dy:rad+dy+H, rad+dx:rad+dx+W]
                na = pa[rad+dy:rad+dy+H, rad+dx:rad+dx+W]
                nn = pn[rad+dy:rad+dy+H, rad+dx:rad+dx+W]
                sw = np.exp(-(dx**2+dy**2)/(2*ss**2))
                cw = np.exp(-((na-alb[:,:,c])**2)/(2*sr**2))
                nw = np.exp(-((nn-nrm[:,:,c])**2)/(2*sn**2))
                w = sw * cw * nw
                out += w * nh; norm += w
        bil[:,:,c] = out / (norm + 1e-9)
    ldr = tm_reinhard(bil, exposure=1.2, white=5.0)
    return bloom(ldr, th=bloom_th, st=bloom_st)

# ============================================================
# 11. 评估指标
# ============================================================
def evaluate(img, gt):
    mse = np.mean((img - gt)**2)
    lum = 0.2126*img[:,:,0] + 0.7152*img[:,:,1] + 0.0722*img[:,:,2]
    lum_gt = 0.2126*gt[:,:,0] + 0.7152*gt[:,:,1] + 0.0722*gt[:,:,2]
    ssim_v = ssim_numpy(lum, lum_gt, data_range=1.0)
    grad = np.sqrt(np.gradient(lum, axis=0)**2 + np.gradient(lum, axis=1)**2)
    grad_gt = np.sqrt(np.gradient(lum_gt, axis=0)**2 + np.gradient(lum_gt, axis=1)**2)
    return mse, ssim_v, grad, grad_gt

# ============================================================
# 12. 主执行
# ============================================================
if __name__ == "__main__":
    W, H = 2048, 1536
    N_FRAMES = 8
    print("=" * 70)
    print("DyVeRT 双层 BVH + 骨骼动画 + GIF + 评估")
    print("分辨率: %dx%d (%.2fM 像素) | %d 帧 | 1 SPP | CPU" % (W, H, W*H/1e6, N_FRAMES))
    print("=" * 70)

    # 构建基础场景
    print("[1/6] 构建场景...")
    t0 = time.time()
    objects, skel, person_obj = build_scene(0.0)
    scene_t = time.time() - t0
    print("  场景构建: %.2fs | 对象数: %d" % (scene_t, len(objects)))

    print("[2/6] 自动相机...")
    cam_pos, lookat = auto_camera(objects)
    print("  Camera: %s -> LookAt: %s" % (str(cam_pos.round(2)), str(lookat.round(2))))

    # 渲染动画帧
    print("[3/6] 渲染 %d 帧动画..." % N_FRAMES)
    anim_frames = []
    anim_raw = []
    for fi in range(N_FRAMES):
        print("  帧 %d/%d..." % (fi+1, N_FRAMES), end="", flush=True)
        t_anim = fi * 0.25
        # 动画人物
        skel.animate(t_anim)
        person_obj.rebuild()
        # 微扰相机模拟手持
        np.random.seed(42 + fi)
        j_cam = cam_pos + np.random.normal(0, 0.03, 3)
        j_look = lookat + np.random.normal(0, 0.015, 3)
        hdr, alb, nrm, tlas_t, trace_t = render_frame(j_cam, j_look, objects, W, H, noise_sigma=0.05)
        final = pipe_unity(hdr, alb, nrm)
        anim_frames.append(final)
        anim_raw.append(hdr)
        print(" 渲染%.2fs" % (tlas_t + trace_t))

    # 保存 GIF (用 Pillow, 降低分辨率到 512 宽度以控制文件大小)
    print("[4/6] 生成 GIF...")
    gif_path = "dyvert_anim.gif"
    gif_frames = []
    for f in anim_frames:
        # 下采样到 512x384 用于 GIF
        small = f[::3, ::4, :]  # 近似降采样
        small = np.clip(small, 0, 1)
        gif_frames.append(Image.fromarray((small * 255).astype(np.uint8)))
    gif_frames[0].save(gif_path, save_all=True, append_images=gif_frames[1:], duration=400, loop=0)
    print("  GIF: %s (%d frames)" % (gif_path, N_FRAMES))

    # GT: 第4帧 16x 累积
    print("[5/6] GT (第4帧 16x 累积)...")
    gt_acc = np.zeros((H, W, 3), dtype=np.float32)
    for s in range(16):
        np.random.seed(100 + s)
        skel.animate(3 * 0.25)
        person_obj.rebuild()
        j_cam = cam_pos + np.random.normal(0, 0.02, 3)
        j_look = lookat + np.random.normal(0, 0.01, 3)
        hdr, _, _, _, _ = render_frame(j_cam, j_look, objects, W, H, noise_sigma=0.02)
        gt_acc += hdr
    gt_hdr = gt_acc / 16
    gt_final = pipe_unity(gt_hdr, alb, nrm)  # 复用最后一帧的 alb/nrm 近似

    # 评估第4帧
    print("[6/6] 评估 (第4帧 vs GT)...")
    frame4 = anim_frames[3]
    mse, ssim_v, grad, grad_gt = evaluate(frame4, gt_final)
    gdiff = abs(grad.mean() - grad_gt.mean())

    print("")
    print("=" * 70)
    print("评估结果 [第4帧]")
    print("=" * 70)
    print("%-15s %12s" % ("指标", "数值"))
    print("%-15s %12.6f" % ("MSE", mse))
    print("%-15s %12.4f" % ("SSIM", ssim_v))
    print("%-15s %12.5f" % ("梯度均值", grad.mean()))
    print("%-15s %12.5f" % ("GT梯度均值", grad_gt.mean()))
    print("%-15s %12.5f" % ("梯度差", gdiff))
    print("=" * 70)

    # 保存对比图
    print("保存对比图...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for ax, img, title in zip(axes, [gt_final, frame4, anim_raw[3]], ["GT (16x)", "Frame 4 (1spp)", "Raw HDR"]):
        disp = np.clip(img, 0, 1)
        if title == "Raw HDR":
            disp = np.clip(disp / (disp.max() + 1e-6), 0, 1) ** 0.45
        ax.imshow(disp)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.axis("off")
    plt.tight_layout()
    plt.savefig("dyvert_compare.png", dpi=150, bbox_inches="tight")
    plt.close()

    # 保存最终帧
    plt.imsave("dyvert_v95_final.png", np.clip(anim_frames[-1], 0, 1))

    print("全部完成!")
    print("  - dyvert_anim.gif")
    print("  - dyvert_compare.png")
    print("  - dyvert_final.png")
