#!/usr/bin/env python3
"""
读取HDF5文件并验证其结构是否符合 resize_restart.cpp 的要求
"""
import h5py
import sys
import os

def check_h5_structure(filename):
    """
    检查HDF5文件结构，对比resize_restart.cpp中需要读取的字段
    """
    print(f"{'='*70}")
    print(f"检查文件: {filename}")
    print(f"{'='*70}\n")
    
    try:
        with h5py.File(filename, 'r') as f:
            # 根据resize_restart.cpp第61-130行，列出所有需要的字段
            required_fields = {
                'root': {
                    # 第71行：version (关键!)
                    'version': {'location': '/', 'type': 'string', 'critical': True},
                    # 第78-80行：网格尺寸
                    'n1': {'location': '/', 'type': 'int', 'critical': True},
                    'n2': {'location': '/', 'type': 'int', 'critical': True},
                    'n3': {'location': '/', 'type': 'int', 'critical': True},
                    # 第90-95行：可选的精确边界
                    'x1Min': {'location': '/', 'type': 'float', 'critical': False},
                    'x1Max': {'location': '/', 'type': 'float', 'critical': False},
                    'x2Min': {'location': '/', 'type': 'float', 'critical': False},
                    'x2Max': {'location': '/', 'type': 'float', 'critical': False},
                    'x3Min': {'location': '/', 'type': 'float', 'critical': False},
                    'x3Max': {'location': '/', 'type': 'float', 'critical': False},
                    # 第99-100行：球坐标参数（如果x1Min不存在）
                    'Rin': {'location': '/', 'type': 'float', 'critical': False},
                    'Rout': {'location': '/', 'type': 'float', 'critical': False},
                    # 第109-110行：黑洞参数
                    'a': {'location': '/', 'type': 'float', 'critical': False},
                    'hslope': {'location': '/', 'type': 'float', 'critical': False},
                    # 第116-119行：其他参数
                    'gam': {'location': '/', 'type': 'float', 'critical': True},
                    't': {'location': '/', 'type': 'float', 'critical': True},
                    'dt': {'location': '/', 'type': 'float', 'critical': True},
                    'tf': {'location': '/', 'type': 'float', 'critical': True},
                    # 第123行：电子相关（可选）
                    'game': {'location': '/', 'type': 'float', 'critical': False},
                    # 第367行：主数据数组
                    'p': {'location': '/', 'type': 'array', 'critical': True},
                }
            }
            
            print("1. 检查根目录 (/) 下的必需字段:")
            print("-" * 70)
            
            missing_critical = []
            missing_optional = []
            found_fields = []
            
            for field_name, field_info in required_fields['root'].items():
                exists_in_root = field_name in f
                exists_in_header = f'/header/{field_name}' in f
                
                # 特殊处理：p 数组可能叫 prims
                if field_name == 'p':
                    exists_as_prims = 'prims' in f
                else:
                    exists_as_prims = False
                
                status = "✗ 缺失"
                detail = ""
                
                if exists_in_root:
                    status = "✓ 存在"
                    detail = f"在根目录 /"
                    found_fields.append(field_name)
                    
                    # 显示值（如果不是数组）
                    if field_info['type'] != 'array':
                        try:
                            value = f[field_name][()]
                            if isinstance(value, bytes):
                                value = value.decode('utf-8')
                            detail += f", 值: {value}"
                        except:
                            pass
                    else:
                        # 对于数组，显示形状
                        try:
                            shape = f[field_name].shape
                            dtype = f[field_name].dtype
                            detail += f", shape: {shape}, dtype: {dtype}"
                        except:
                            pass
                
                elif exists_as_prims:
                    status = "⚠ 警告"
                    detail = f"数组名为 'prims' 而非 'p'（程序期望 'p'）"
                    try:
                        shape = f['prims'].shape
                        dtype = f['prims'].dtype
                        detail += f", shape: {shape}, dtype: {dtype}"
                    except:
                        pass
                            
                elif exists_in_header:
                    if field_info['critical']:
                        status = "⚠ 警告"
                        detail = f"仅在 /header/ 下找到（程序期望在根目录）"
                    else:
                        status = "✓ 存在"
                        detail = f"在 /header/ 下"
                    
                    # 显示值
                    if field_info['type'] != 'array':
                        try:
                            value = f[f'/header/{field_name}'][()]
                            if isinstance(value, bytes):
                                value = value.decode('utf-8')
                            detail += f", 值: {value}"
                        except:
                            pass
                else:
                    if field_info['critical']:
                        status = "✗ 缺失"
                        detail = "关键字段缺失！"
                        missing_critical.append(field_name)
                    else:
                        status = "○ 可选"
                        detail = "可选字段，不影响运行"
                        missing_optional.append(field_name)
                
                critical_mark = "[必需]" if field_info['critical'] else "[可选]"
                print(f"{status} {critical_mark:8} {field_name:15} {detail}")
            
            print("\n" + "="*70)
            print("2. 文件结构摘要:")
            print("-" * 70)
            
            def print_structure(group, prefix="", max_depth=3, current_depth=0):
                """递归打印HDF5结构"""
                if current_depth >= max_depth:
                    return
                    
                items = []
                try:
                    for key in group.keys():
                        items.append(key)
                except:
                    pass
                    
                for key in sorted(items):
                    item = group[key]
                    if isinstance(item, h5py.Dataset):
                        shape_info = f"shape: {item.shape}, dtype: {item.dtype}"
                        print(f"{prefix}├── {key} (dataset, {shape_info})")
                    elif isinstance(item, h5py.Group):
                        print(f"{prefix}├── {key}/ (group)")
                        print_structure(item, prefix + "│   ", max_depth, current_depth + 1)
            
            print("文件层次结构:")
            print_structure(f, "", max_depth=3)
            
            print("\n" + "="*70)
            print("3. 诊断结果:")
            print("-" * 70)
            
            if missing_critical:
                print(f"✗ 发现 {len(missing_critical)} 个缺失的关键字段:")
                for field in missing_critical:
                    print(f"  - {field}")
                    # 检查是否在header中
                    if f'/header/{field}' in f:
                        print(f"    提示: 该字段存在于 /header/{field}，需要复制到根目录 /{field}")
            else:
                print("✓ 所有关键字段都存在（或在可接受位置）")
            
            if missing_optional:
                print(f"\n○ {len(missing_optional)} 个可选字段缺失（不影响运行）:")
                for field in missing_optional:
                    print(f"  - {field}")
            
            # 特别检查version字段
            print("\n" + "="*70)
            print("4. 特别检查 - version 字段 (最常见的错误原因):")
            print("-" * 70)
            
            version_in_root = 'version' in f
            version_in_header = '/header/version' in f
            
            if version_in_root:
                version_value = f['version'][()]
                if isinstance(version_value, bytes):
                    version_value = version_value.decode('utf-8')
                print(f"✓ /version 存在，值: {version_value}")
                print("  文件结构正确！")
            elif version_in_header:
                version_value = f['/header/version'][()]
                if isinstance(version_value, bytes):
                    version_value = version_value.decode('utf-8')
                print(f"✗ /version 不存在于根目录")
                print(f"⚠ /header/version 存在，值: {version_value}")
                print(f"\n  这就是错误的原因！")
                print(f"  resize_restart.cpp 第71行期望在根目录读取 /version")
                print(f"  但文件中 version 位于 /header/version")
                print(f"\n  解决方案: 运行 fix_h5.py 修复此问题")
                print(f"  命令: python3 fix_h5.py {os.path.basename(filename)}")
            else:
                print("✗ version 字段完全不存在！")
                print("  文件可能已损坏或格式不正确")
            
            print("\n" + "="*70)
            
            # 返回是否有关键问题
            return len(missing_critical) == 0 and version_in_root
            
    except Exception as e:
        print(f"\n✗ 错误: 无法读取文件")
        print(f"  详细信息: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python reda_h5.py <hdf5文件>")
        print("示例: python reda_h5.py torus.out0.05000.h5")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    if not os.path.exists(filename):
        print(f"错误: 文件不存在: {filename}")
        sys.exit(1)
    
    is_valid = check_h5_structure(filename)
    
    if is_valid:
        print("\n✓✓✓ 文件结构验证通过！可以正常使用。")
        sys.exit(0)
    else:
        print("\n✗✗✗ 文件结构存在问题，需要修复后才能使用。")
        sys.exit(1)

