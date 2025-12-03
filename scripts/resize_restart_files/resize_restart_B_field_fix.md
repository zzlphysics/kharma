# resize_restart 磁场数据丢失问题的诊断与修复

## 问题描述

在使用 `resize_restart` 问题生成器从 iharm3d HDF5 文件重启时，所有流体变量（`rho`, `u`, `uvec`）都正确读取，但磁场分量（`B1`, `B2`, `B3`）在输出文件中全部为零。

## 根本原因

### 1. Parthenon 的变量分配机制

KHARMA 使用 Parthenon 的变量管理系统，其中有两类关键的元数据标志：

- **`Derived`**：派生变量，由其他变量计算得到，Parthenon 不为其分配持久存储
- **`Independent`**：独立变量，Parthenon 为其分配持久存储并管理生命周期

### 2. prims.B 的元数据配置

在 `b_ct.cpp` 中，`prims.B`（cell-centered 原始磁场）被定义为：
```cpp
Metadata::Derived + Metadata::Restart
```

这个配置的含义：
- `Derived`：表示这是从 face-centered `cons.fB` 派生的
- `Restart`：表示在 restart 时需要分配空间

### 3. 初始化流程导致的数据丢失

KHARMA 的初始化流程：

```
ProblemGenerator (problem.cpp)
  └─> ReadIharmRestart (resize_restart.cpp)
      └─> 读取 HDF5 数据到 prims.B  ✓ 数据正确
      └─> 返回

Mesh::Initialize (Parthenon 内部)
  └─> SetAllVariablesToInitialized()
  └─> PreCommFillDerived()
      └─> 重新分配所有 Derived 字段  ✗ prims.B 被清零！

PostInitialize (post_initialize.cpp)
  └─> DangerousPtoU: prims.B → cons.fB
      └─> 但 prims.B 已经是零了  ✗
```

**关键问题**：`prims.B` 标记为 `Derived`，在 `ProblemGenerator` 和 `PostInitialize` 之间，Parthenon 会重新分配 `Derived` 变量，导致我们在 `ReadIharmRestart` 中填充的数据被清零。

## 解决方案

### 核心思路

使用一个临时的 `Independent + OneCopy` 缓存字段来保存磁场数据，该字段不会被 Parthenon 重新分配：

```
ReadIharmRestart:
  prims.B (Derived) ←─ 从 HDF5 读取
  prims.B_cache (Independent) ←─ 同时保存副本

[Parthenon 重新分配 Derived 字段]
  prims.B 被清零
  prims.B_cache 保持不变 ✓

PostInitialize:
  prims.B ←─ 从 prims.B_cache 恢复
  cons.fB ←─ DangerousPtoU(prims.B)
```

### 具体实现

#### 1. 添加 B_cache 字段 (`b_ct.cpp`)

```cpp
if (pin->GetOrAddBoolean("b_field", "restart_from_prims", false)) {
    flags_prim.push_back(Metadata::Restart);
    
    // 添加临时缓存，防止 resize_restart 期间数据丢失
    std::vector<MetadataFlag> flags_cache = {
        Metadata::Real, 
        Metadata::Cell, 
        Metadata::Independent,  // 关键：不会被重新分配
        Metadata::OneCopy,      // 关键：只存在于 host，节省内存
        Metadata::Vector
    };
    std::vector<int> s_vector_cache({NVEC});
    Metadata m_cache = Metadata(flags_cache, s_vector_cache);
    pkg->AddField("prims.B_cache", m_cache);
}
```

#### 2. 保存数据到缓存 (`resize_restart.cpp`)

```cpp
GridVector B_P = rc->Get("prims.B").data;
GridVector B_cache = rc->Get("prims.B_cache").data;

auto B_host = B_P.GetHostMirror();
auto B_cache_host = B_cache.GetHostMirror();

// 插值循环中同时填充两个字段
VLOOP B_host(v, k, j, i) = /* 从 HDF5 读取或插值 */;
VLOOP B_cache_host(v, k, j, i) = B_host(v, k, j, i);

// 拷贝到设备
B_P.DeepCopy(B_host);
B_cache.DeepCopy(B_cache_host);  // OneCopy，实际只在 host
```

#### 3. 从缓存恢复数据 (`post_initialize.cpp`)

```cpp
// 检查是否有缓存（resize_restart 问题）
bool has_cache = false;
if (prob_name == "resize_restart") {
    for (int i = 0; i < pmesh->GetNumMeshBlocksThisRank(); i++) {
        auto rc = pmesh->block_list[i]->meshblock_data.Get("base");
        if (rc->IsAllocated("prims.B_cache")) {
            has_cache = true;
            break;
        }
    }
}

if (has_cache) {
    // 从缓存恢复到 prims.B
    for (int i = 0; i < pmesh->GetNumMeshBlocksThisRank(); i++) {
        auto pmb = pmesh->block_list[i];
        auto rc = pmb->meshblock_data.Get("base");
        auto B_P = rc->Get("prims.B").data;
        auto B_cache = rc->Get("prims.B_cache").data;
        
        auto B_cache_host = B_cache.GetHostMirrorAndCopy();
        B_P.DeepCopy(B_cache_host);  // 直接拷贝到设备
    }
    Kokkos::fence();
}

// 然后调用 DangerousPtoU 将 prims.B 转换为 cons.fB
B_CT::DangerousPtoU(md.get(), IndexDomain::interior, false);
```

#### 4. 跳过不必要的 BlockPtoU (`problem.cpp`)

```cpp
// resize_restart 已经手动填充了 prims.B，不需要从 cons 计算
if (prob != "resize_restart") {
    Flux::BlockPtoU(rc.get(), IndexDomain::entire);
}
```

## 调试过程中遇到的问题

### 问题 1：段错误访问 cons.fB

**错误**：尝试直接对 `cons.fB` 调用 `GetHostMirrorAndCopy()`
```cpp
auto B_U_host = rc->Get("cons.fB").data.GetHostMirrorAndCopy();  // 段错误
```

**原因**：`cons.fB` 是 `Face` 类型（`Metadata::Face`），不是普通的 4D 数组，需要通过 `PackVariables` 访问并用 F1、F2、F3 索引。

**解决**：删除这段调试代码。

### 问题 2：restart_from_prims = 0 时找不到 prims.B_cache

**错误**：`Couldn't find variable 'prims.B_cache'`

**原因**：`prims.B_cache` 只在 `restart_from_prims = 1` 时创建，但 `resize_restart.cpp` 无条件访问它。

**解决**：在参数文件中强制要求 `restart_from_prims = 1`，并在代码中添加文档说明。

## 使用要求

**重要**：使用 `resize_restart` 必须在参数文件中设置：

```
<b_field>
restart_from_prims = 1
```

否则会因为找不到 `prims.B_cache` 而失败。

## 内存开销

`prims.B_cache` 使用 `OneCopy` 标志，只存在于 host 端，不占用 GPU 内存。对于典型的 288×128×128 网格：
- 每个 block 的 B_cache: 3 × 296 × 40 × 136 × 8 bytes ≈ 10 MB
- 4 个 blocks: ≈ 40 MB

这是可接受的临时开销。

## 未来改进方向

### 选项 1：使 prims.B 为 Independent

将 `prims.B` 改为 `Independent + Restart`，这样就不需要缓存。但这会改变 KHARMA 的基本设计（`prims.B` 应该从 `cons.fB` 派生）。

### 选项 2：修改 Parthenon

在 Parthenon 中添加 `Derived + Persistent` 标志，允许派生变量在重新分配时保留数据。

### 选项 3：改用 KHARMA 原生 restart

从 KHARMA 的 `.phdf` 文件重启，而不是 iharm3d HDF5 文件，这样就可以使用标准的 restart 机制。

## 相关文件

- `kharma/b_ct/b_ct.cpp`: 定义 `prims.B` 和 `prims.B_cache` 字段
- `kharma/prob/resize_restart.cpp`: 读取 HDF5 并保存到缓存
- `kharma/prob/post_initialize.cpp`: 从缓存恢复数据
- `kharma/prob/problem.cpp`: 跳过 resize_restart 的 BlockPtoU

## 参考

- Parthenon 文档: https://github.com/parthenon-hpc-lab/parthenon
- KHARMA 仓库: https://github.com/AFD-Illinois/kharma

