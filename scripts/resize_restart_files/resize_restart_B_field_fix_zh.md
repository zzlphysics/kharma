# resize_restart 磁场数据丢失问题修复说明（简要版）

## 问题

使用 `resize_restart` 从 iharm3d HDF5 文件读取数据时，流体变量正常但磁场 B 全为零。

## 原因

`prims.B` 被标记为 `Derived`（派生变量），Parthenon 会在 `ProblemGenerator` 和 `PostInitialize` 之间重新分配它，导致在 `ReadIharmRestart` 中填充的数据被清零。

## 解决方案

使用临时缓存 `prims.B_cache`（标记为 `Independent + OneCopy`）保存磁场数据：

1. **读取阶段** (`resize_restart.cpp`)：同时填充 `prims.B` 和 `prims.B_cache`
2. **Parthenon 重新分配**：`prims.B` 被清零，但 `prims.B_cache` 保持不变
3. **恢复阶段** (`post_initialize.cpp`)：从 `prims.B_cache` 恢复到 `prims.B`
4. **转换阶段**：`DangerousPtoU` 将 `prims.B` 转换为 `cons.fB`

## 修改文件

- `kharma/b_ct/b_ct.cpp`：添加 `prims.B_cache` 字段定义
- `kharma/prob/resize_restart.cpp`：保存数据到缓存
- `kharma/prob/post_initialize.cpp`：从缓存恢复数据
- `kharma/prob/problem.cpp`：跳过 `resize_restart` 的 `BlockPtoU` 调用

## 使用要求

**必须在参数文件中设置**：
```
<b_field>
restart_from_prims = 1
```

## 内存开销

`prims.B_cache` 使用 `OneCopy` 标志，只在 CPU 端，典型 288×128×128 网格约占 40 MB 临时内存。

## 详细文档

参见 `resize_restart_B_field_fix.md` 了解完整技术细节和调试过程。

