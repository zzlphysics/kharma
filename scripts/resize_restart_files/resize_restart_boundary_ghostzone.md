# resize_restart 的边界和 Ghost Zone 处理

## 1. resize_restart 初始化时的边界处理

### 1.1 数据读取阶段 (`ReadIharmRestart`)

**读取范围**：`IndexDomain::entire`
```cpp
const auto& kb = pmb->cellbounds.GetBoundsK(IndexDomain::entire);
const auto& jb = pmb->cellbounds.GetBoundsJ(IndexDomain::entire);
const auto& ib = pmb->cellbounds.GetBoundsI(IndexDomain::entire);
```

这意味着读取：
- **物理区域（interior）**：实际模拟的网格点
- **Ghost zones**：边界外的额外层（通常 `nghost=4`）

### 1.2 边界插值处理

对于球坐标系统 (`is_spherical=true`)，边界有特殊处理：

```cpp
// resize_restart.cpp line 490-506
const bool repeat_x1i = is_spherical;  // 内径边界
const bool repeat_x1o = is_spherical;  // 外径边界
const bool repeat_x2i = is_spherical;  // theta 下边界（极点）
const bool repeat_x2o = is_spherical;  // theta 上边界（极点）

// 边界修正：重复最后一个 zone（等效于最近邻插值）
if (repeat_x1i && gi < 0) { gi = 0; del[1] = 0; }
if (repeat_x1o && gi > n1tot-2) { gi = n1tot - 2; del[1] = 1; }
if (repeat_x2i && gj < 0) { gj = 0; del[2] = 0; }
if (repeat_x2o && gj > n2tot-2) { gj = n2tot - 2; del[2] = 1; }
```

**这意味着**：
- 在 r 方向（径向）：内外边界重复最内/最外层的值
- 在 theta 方向（极角）：极点重复边界值
- 在 phi 方向（方位角）：周期性边界（自动处理）

### 1.3 Dirichlet 边界冻结

在 `ProblemGenerator` 结束时：
```cpp
// problem.cpp line 169
KBoundaries::FreezeDirichletBlock(rc.get());
```

**作用**：
1. 检查每个边界面是否设置为 `dirichlet` 类型
2. 如果是，将当前 ghost zone 的值保存到边界缓存 `Boundaries.[face_name]`
3. 这些值在整个模拟过程中保持不变

**实现细节** (`dirichlet.cpp` line 131-146):
```cpp
void FreezeDirichletBlock(MeshBlockData<Real> *rc) {
    for (int i=0; i < BOUNDARY_NFACES; i++) {
        BoundaryFace bface = (BoundaryFace) i;
        if (is_dirichlet && is_physical_boundary) {
            // 将 ghost zone 的值保存到缓存
            SetDomainDirichlet(rc, domain, false);
        }
    }
}
```

## 2. 正常模拟过程中的边界处理

### 2.1 时间步进中的边界更新流程

**典型的时间步任务顺序** (`kharma_step.cpp`):

```
1. StartReceive / StartSend        ← 启动 MPI 通信
2. FluxDivergence                  ← 计算通量散度
3. FinishReceive                   ← 完成 MPI 通信（同步保守变量）
4. UtoP (entire domain)            ← U → P 转换（包括 ghost zones）
5. ApplyFloors                     ← 应用物理下限
6. FixUtoP                         ← 修复失败的反演
7. ApplyBoundaryConditions         ← 应用域边界条件
8. PtoU (entire domain)            ← P → U 转换
```

### 2.2 边界条件应用 (`ApplyBoundaryConditions`)

**Parthenon 调用**：
```cpp
// kharma_step.cpp line 244-246
auto t_set_bc = tl.AddTask(t_fix_p, parthenon::ApplyBoundaryConditionsOnCoarseOrFineMD, 
                           md_sub_step_final, false);
```

这会调用 KHARMA 的边界包装函数：
```cpp
// boundaries.cpp: ApplyBoundary()
1. 检查边界类型（dirichlet/outflow/reflecting/periodic）
2. 对于 dirichlet：从缓存恢复值到 ghost zones
3. 对于 outflow/reflecting：调用相应的处理函数
4. 特殊处理：极点平均、通量修正等
```

### 2.3 Dirichlet 边界应用

**每个时间步** (`dirichlet.cpp` line 96-102):
```cpp
if (set) {
    // 保存当前值到缓存（初始化时）
    bound(...) = q(...);
} else {
    // 从缓存恢复到 ghost zones（每个时间步）
    q(...) = bound(...);
}
```

## 3. 关键区别对比

| 方面 | resize_restart 初始化 | 正常模拟 |
|------|----------------------|----------|
| **Ghost zone 来源** | 从 HDF5 文件插值 | MPI 通信 + 边界条件 |
| **边界处理时机** | `ProblemGenerator` 结束时 | 每个时间步 |
| **Dirichlet 边界** | 一次性冻结 | 每步从缓存恢复 |
| **插值方法** | 线性插值（或最近邻） | 无插值，直接复制/计算 |
| **域** | `IndexDomain::entire` | `IndexDomain::entire` |

## 4. 一致性分析

### 4.1 相同之处 ✓

1. **域范围**：两者都处理 `entire` domain（物理 + ghost zones）
2. **Dirichlet 机制**：都使用相同的边界缓存机制
3. **冻结操作**：`FreezeDirichletBlock` 在两种情况下行为一致
4. **最终状态**：初始化结束时，ghost zones 都已正确填充

### 4.2 不同之处 ⚠️

1. **初始填充方式**：
   - **resize_restart**：从 HDF5 插值（可能有微小差异）
   - **正常模拟**：从邻居块精确复制

2. **边界类型处理**：
   - **resize_restart**：球坐标边界使用重复/最近邻
   - **正常模拟**：根据边界类型（outflow/reflecting 等）计算

3. **磁场处理**：
   - **resize_restart**：cell-centered B 需要转换为 face-centered
   - **正常模拟**：face-centered B 直接同步

### 4.3 潜在问题点

1. **插值误差**：
   - 如果源网格和目标网格分辨率不同，边界插值可能引入误差
   - 特别是在磁场的情况下，`DangerousPtoU` 会再次插值

2. **边界一致性**：
   - 如果 HDF5 文件中的边界值不是正确的边界条件结果
   - Ghost zones 可能包含不一致的数据

3. **MPI 边界**：
   - `resize_restart` 在单个 block 内插值，不涉及 MPI 通信
   - 如果 block 分解改变，可能需要 `PostInitialize` 中的同步

## 5. PostInitialize 中的同步

### 5.1 边界同步 (`post_initialize.cpp`)

```cpp
// line 223-224
KBoundaries::FreezeDirichlet(md);
KHARMADriver::SyncAllBounds(md);
```

**作用**：
1. **FreezeDirichlet**：确保所有 Dirichlet 边界已冻结
2. **SyncAllBounds**：
   - 同步所有 MPI 边界
   - 应用所有物理边界条件
   - 填充所有 ghost zones

### 5.2 这解决了什么问题？

- **MPI 一致性**：确保跨进程边界的数据一致
- **B 场同步**：在 `DangerousPtoU` 之后，`cons.fB` 需要在边界上同步
- **完整性检查**：在开始模拟前确保所有边界正确

## 6. 结论

### resize_restart 的边界处理是**基本一致的**：

✓ **域范围一致**：都处理整个 domain 包括 ghost zones
✓ **Dirichlet 机制一致**：使用相同的冻结/恢复机制  
✓ **最终状态一致**：`PostInitialize` 中的同步确保一致性

### 但有一些**初始化阶段的差异**：

⚠️ **插值方法**：resize_restart 使用插值，可能引入小误差
⚠️ **边界来源**：从文件读取 vs 从物理计算  
⚠️ **时序**：一次性设置 vs 每步更新

### 建议：

1. **验证边界值**：检查 HDF5 文件中 ghost zones 是否合理
2. **B 场清理**：resize_restart 后建议运行 `b_cleanup` 来减小 divB
3. **测试一致性**：比较 resize_restart 和正常 restart 的前几步演化

## 7. 相关代码位置

- `kharma/prob/resize_restart.cpp`: HDF5 读取和插值
- `kharma/prob/problem.cpp`: `FreezeDirichletBlock` 调用
- `kharma/prob/post_initialize.cpp`: `FreezeDirichlet` 和边界同步
- `kharma/boundaries/dirichlet.cpp`: Dirichlet 边界实现
- `kharma/driver/kharma_step.cpp`: 时间步中的边界更新流程

