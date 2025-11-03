#!/usr/bin/env python3
"""
修复HDF5文件：修复KHARMA resize_restart所需的字段位置和命名
"""
import h5py
import sys
import shutil
import os
import numpy as np

def fix_h5_file(filename):
    """
    修复HDF5文件以符合resize_restart.cpp的要求：
    1. 将 /header/version 复制到 /version
    2. 将 prims 重命名为 p，调整维度顺序并转换为双精度
       - 维度: iharm3d格式 (n1,n2,n3,nprim) → KHARMA格式 (nprim,n3,n2,n1)
       - 例如: (288,128,128,8) → (8,128,128,288)
       - 精度: float32 → float64
    3. 将 /header/ 下的必需字段复制到根目录
    4. 计算 Rin = exp(startx1) 和 Rout = r_out
    
    Args:
        filename: HDF5文件路径
    """
    # 先备份原文件
    backup_file = filename + '.backup'
    if not os.path.exists(backup_file):
        print(f"备份原文件到: {backup_file}")
        shutil.copy2(filename, backup_file)
    else:
        print(f"备份文件已存在: {backup_file}")
    
    # resize_restart.cpp需要的根目录字段
    required_root_fields = ['n1', 'n2', 'n3', 'gam', 't', 'dt', 'tf', 'a', 'hslope']
    
    fixed_count = 0
    
    try:
        # 打开文件（读写模式）
        with h5py.File(filename, 'r+') as f:
            print(f"\n打开文件: {filename}")
            print("="*70)
            
            # 1. 修复 version 字段
            print("\n1. 检查 version 字段...")
            if '/version' in f:
                print("  ✓ /version 已存在")
            elif '/header/version' in f:
                version_data = f['/header/version'][()]
                f.create_dataset('/version', data=version_data)
                if isinstance(version_data, bytes):
                    version_data = version_data.decode('utf-8')
                print(f"  ✓ 已从 /header/version 复制到 /version: {version_data}")
                fixed_count += 1
            else:
                print("  ✗ 警告: version 字段不存在")
            
            # 2. 修复 prims -> p (并转换为双精度和调整维度顺序)
            print("\n2. 检查主数据数组...")
            if '/p' in f:
                # 检查现有的 /p 是否需要转换精度和维度
                p_dtype = f['/p'].dtype
                p_shape = f['/p'].shape
                print(f"  ✓ /p 已存在，shape={p_shape}, dtype={p_dtype}")
                
                needs_fixing = False
                p_data = f['/p'][...]
                
                # 检查维度顺序：如果是 (n1,n2,n3,nprim) 需要转为 (nprim,n3,n2,n1)
                if len(p_shape) == 4 and p_shape[3] < p_shape[0]:
                    # 假设最后一维较小的是 nprim
                    print(f"    ⚠ /p 维度顺序错误: {p_shape}")
                    print(f"    ⚠ iharm3d格式 (n1,n2,n3,nprim) → KHARMA格式 (nprim,n3,n2,n1)")
                    # 转置: (n1,n2,n3,nprim) -> (nprim,n3,n2,n1)
                    # 即: 维度 [0,1,2,3] -> [3,2,1,0]
                    p_data = np.transpose(p_data, (3, 2, 1, 0))
                    # 确保数组是C-contiguous（内存连续），这对HDF5读取很重要
                    p_data = np.ascontiguousarray(p_data)
                    print(f"    ✓ 已调整维度顺序: {p_shape} → {p_data.shape}")
                    print(f"    ✓ 内存布局: C-contiguous")
                    needs_fixing = True
                
                # 检查精度
                if p_dtype == np.float32:
                    print(f"    ⚠ /p 是单精度(float32)，转换为双精度(float64)")
                    p_data = p_data.astype(np.float64)
                    needs_fixing = True
                
                if needs_fixing:
                    del f['/p']
                    # 重要：显式指定chunks=None以禁用chunking，确保数据连续存储
                    f.create_dataset('/p', data=p_data, dtype=np.float64, chunks=None)
                    print(f"    ✓ 已更新 /p: shape={p_data.shape}, dtype=float64, chunks=None")
                    fixed_count += 1
                    
            elif '/prims' in f:
                print("  ⚠ 发现 /prims，需要重命名为 /p")
                # HDF5不支持直接重命名，需要复制然后删除
                prims_dtype = f['/prims'].dtype
                prims_data = f['/prims'][...]
                prims_shape = prims_data.shape
                print(f"    读取 prims 数据: shape={prims_shape}, dtype={prims_dtype}")
                
                # 调整维度顺序：iharm3d (n1,n2,n3,nprim) → KHARMA (nprim,n3,n2,n1)
                if len(prims_shape) == 4:
                    print(f"    ⚠ 调整维度顺序: iharm3d格式 (n1,n2,n3,nprim) → KHARMA格式 (nprim,n3,n2,n1)")
                    # 转置: (n1,n2,n3,nprim) -> (nprim,n3,n2,n1)
                    # 即: 维度 [0,1,2,3] -> [3,2,1,0]
                    prims_data = np.transpose(prims_data, (3, 2, 1, 0))
                    # 确保数组是C-contiguous（内存连续），这对HDF5读取很重要
                    prims_data = np.ascontiguousarray(prims_data)
                    print(f"    ✓ 维度已调整: {prims_shape} → {prims_data.shape}")
                    print(f"    ✓ 内存布局: C-contiguous")
                
                # 检查是否需要转换精度
                if prims_dtype == np.float32:
                    print(f"    ⚠ prims 是单精度(float32)，转换为双精度(float64)")
                    prims_data = prims_data.astype(np.float64)
                    target_dtype = np.float64
                else:
                    print(f"    ✓ prims 已经是双精度(float64)")
                    target_dtype = prims_dtype
                
                # 创建新数据集 /p（禁用chunking以确保连续存储）
                f.create_dataset('/p', data=prims_data, dtype=target_dtype, chunks=None)
                print(f"    ✓ 已创建 /p: shape={prims_data.shape}, dtype={target_dtype}, chunks=None")
                
                # 删除旧的 /prims
                del f['/prims']
                print(f"    ✓ 已删除 /prims")
                fixed_count += 1
            else:
                print("  ✗ 错误: 主数据数组 p 或 prims 都不存在")
            
            # 3. 复制其他必需字段从 /header/ 到根目录
            print("\n3. 检查其他必需字段...")
            for field in required_root_fields:
                if f'/{field}' in f:
                    print(f"  ✓ /{field} 已存在")
                elif f'/header/{field}' in f:
                    field_data = f[f'/header/{field}'][()]
                    f.create_dataset(f'/{field}', data=field_data)
                    print(f"  ✓ 已从 /header/{field} 复制到 /{field}: {field_data}")
                    fixed_count += 1
                else:
                    print(f"  ⚠ {field} 在根目录和header中都不存在")
            
            # 4. 计算并创建 Rin 和 Rout
            print("\n4. 检查并计算 Rin 和 Rout...")
            
            # 处理 Rin: Rin = exp(startx1)
            if '/Rin' in f:
                print(f"  ✓ /Rin 已存在")
            else:
                # 尝试从 startx1 计算
                startx1 = None
                if '/header/geom/startx1' in f:
                    startx1 = f['/header/geom/startx1'][()]
                    source = '/header/geom/startx1'
                elif '/header/startx1' in f:
                    startx1 = f['/header/startx1'][()]
                    source = '/header/startx1'
                elif '/startx1' in f:
                    startx1 = f['/startx1'][()]
                    source = '/startx1'
                
                if startx1 is not None:
                    Rin = np.exp(startx1)
                    f.create_dataset('/Rin', data=Rin)
                    print(f"  ✓ 从 {source} 计算 Rin = exp({startx1}) = {Rin}")
                    fixed_count += 1
                else:
                    print(f"  ✗ 无法找到 startx1 来计算 Rin")
            
            # 处理 Rout: 直接从 r_out 复制
            if '/Rout' in f:
                print(f"  ✓ /Rout 已存在")
            else:
                # 尝试从不同位置找 r_out 或 rout
                rout_value = None
                rout_source = None
                
                if '/r_out' in f:
                    rout_value = f['/r_out'][()]
                    rout_source = '/r_out'
                elif '/header/r_out' in f:
                    rout_value = f['/header/r_out'][()]
                    rout_source = '/header/r_out'
                elif '/header/geom/fmks/r_out' in f:
                    rout_value = f['/header/geom/fmks/r_out'][()]
                    rout_source = '/header/geom/fmks/r_out'
                elif '/rout' in f:
                    rout_value = f['/rout'][()]
                    rout_source = '/rout'
                elif '/header/rout' in f:
                    rout_value = f['/header/rout'][()]
                    rout_source = '/header/rout'
                
                if rout_value is not None:
                    f.create_dataset('/Rout', data=rout_value)
                    print(f"  ✓ 从 {rout_source} 复制到 /Rout: {rout_value}")
                    fixed_count += 1
                else:
                    print(f"  ✗ 无法找到 r_out 或 rout")
            
            print("\n" + "="*70)
            print(f"修复完成！共修复了 {fixed_count} 个字段")
            
            return True
            
    except Exception as e:
        print(f"\n✗ 错误: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python fix_h5.py <hdf5文件>")
        print("示例: python fix_h5.py torus.out0.05000.h5")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    if not os.path.exists(filename):
        print(f"错误: 文件不存在: {filename}")
        sys.exit(1)
    
    print(f"开始修复文件: {filename}")
    print("="*70)
    
    success = fix_h5_file(filename)
    
    print("\n" + "="*70)
    if success:
        print("✓✓✓ 文件修复成功！")
        print("\n建议：运行 read_h5.py 验证修复结果")
        print(f"命令: python3 read_h5.py {os.path.basename(filename)}")
        sys.exit(0)
    else:
        print("✗✗✗ 文件修复失败")
        print(f"\n可以从备份恢复: {filename}.backup")
        sys.exit(1)

