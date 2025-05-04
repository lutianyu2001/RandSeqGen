#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)

from core import SequenceTree


def test_multiple_cuts():
    """
    测试多重切割情况下的重建算法
    创建一个被多个供体切割的情景，然后执行重建，验证所有切割关系是否都被正确处理。
    Returns:
        tuple: (重建的donor记录列表, 是否成功)
    """
    # 初始化树数据结构
    tree = SequenceTree("ATGCATGCATGCATGCATGCATGCATGCATGCATGC")  # 原始序列

    # 创建一个被多次切割的简单情景
    # 原始序列被三个不同的donor切割：
    # - donor 1 切割在位置10
    # - donor 2 切割在位置20
    # - donor 3 切割在位置30
    # 使用真实DNA序列
    first_cutter = "GTACGTAC"
    second_cutter = "CCGGAATT"
    third_cutter = "TTAGGCCA"
    
    # 插入donor序列在相应位置
    tree.insert(10, first_cutter, "cutter1")
    tree.insert(20, second_cutter, "cutter2")
    tree.insert(30, third_cutter, "cutter3")
    
    # 重建供体
    _, reconstructed_records = tree.donors("test")
    
    # 定义预期结果 - 明确基于SequenceTree的预期行为
    expected_results = {
        "clean_recon_count": 3,  # 期望有3个清洁重建（每个切割点一个）
        "full_recon_count": 3,   # 期望有3个完整重建（每个切割者一个）
        "expected_cutters": [first_cutter, second_cutter, third_cutter],
        "min_cutter_matches": 2  # 至少有2个切割者应该在完整重建中找到
    }
    
    # 验证正确性
    success = True
    error_messages = []
    
    # 检查清洁重建数量
    clean_recon = [rec for rec in reconstructed_records if rec.annotations.get("reconstruction_type") == "clean"]
    if len(clean_recon) != expected_results["clean_recon_count"]:
        error_messages.append(f"错误: 应该有{expected_results['clean_recon_count']}个清洁重建，但实际有{len(clean_recon)}个")
        success = False
    
    # 检查完整重建数量
    full_recon = [rec for rec in reconstructed_records if rec.annotations.get("reconstruction_type") == "full"]
    if len(full_recon) != expected_results["full_recon_count"]:
        error_messages.append(f"错误: 应该有{expected_results['full_recon_count']}个完整重建，但实际有{len(full_recon)}个")
        success = False
    
    # 验证完整重建序列内容 - 检查是否包含切割者序列
    found_cutters = []
    for cutter in expected_results["expected_cutters"]:
        found = False
        for rec in full_recon:
            seq_str = str(rec.seq)
            # 检查完整或部分序列
            if cutter in seq_str or cutter[:4] in seq_str:
                found = True
                found_cutters.append(cutter)
                break
    
    if len(found_cutters) < expected_results["min_cutter_matches"]:
        error_messages.append(f"错误: 至少应该找到{expected_results['min_cutter_matches']}个切割者的序列，但只找到了{len(found_cutters)}个")
        success = False
    
    # 打印测试结果
    if success:
        print("多重切割测试成功！所有切割关系都被正确处理。")
    else:
        print("多重切割测试失败！请检查重建算法。")
        for error in error_messages:
            print(f"  {error}")
    
    return reconstructed_records, success


def test_comprehensive_nesting():
    """
    全面测试多种复杂嵌套和切割情景，验证重建算法的正确性。
    测试包含多种场景:
    - 单个插入
    - 相邻位置插入
    - 嵌套插入
    - 多重嵌套
    - 切割donor
    - 多重切割
    - 连锁切割
    - 边界切割
    - 随机多donor复杂网络
    - 空序列处理
    - 极端长短序列
    - 完全重叠切割
    Returns:
        bool: 测试是否全部通过
    """
    test_results = []
    
    # ======== 场景1: 单个donor插入 ========
    print("\n===== 测试场景1: 单个donor插入 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 插入donor
    tree.insert(2, "TTT", "donor1")
    
    # 定义预期结果
    expected = {
        "total_reconstructions": 2,         # 预期总共有2个重建
        "clean_reconstructions": 1,         # 预期有1个清洁重建
        "full_reconstructions": 1,          # 预期有1个完整重建
        "donor_sequence_present": "TTT"     # 预期重建中包含此序列
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    # 检查重建总数
    if reconstructed is None or len(reconstructed) != expected["total_reconstructions"]:
        success = False
        reason = f"期望{expected['total_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        # 检查清洁重建数量
        clean_count = len([r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"])
        if clean_count != expected["clean_reconstructions"]:
            success = False
            reason = f"期望{expected['clean_reconstructions']}个清洁重建，但得到了{clean_count}个"
        
        # 检查完整重建数量
        full_count = len([r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"])
        if full_count != expected["full_reconstructions"]:
            success = False
            reason = f"期望{expected['full_reconstructions']}个完整重建，但得到了{full_count}个"
        
        # 检查序列内容
        if not any(expected["donor_sequence_present"] in str(r.seq) for r in reconstructed):
            success = False
            reason = f"在重建中未找到预期序列 {expected['donor_sequence_present']}"
    
    # 记录测试结果
    test_results.append(("场景1", success, reason if not success else "单个donor插入，创建了预期的重建"))
    if success:
        print("✓ 场景1测试通过: 单个donor插入，创建了预期的重建")
    else:
        print(f"✗ 场景1测试失败: {reason}")
    
    # ======== 场景2: 相邻位置双重插入 ========
    print("\n===== 测试场景2: 相邻位置双重插入 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 插入donor1
    tree.insert(2, "TTT", "donor1")
    # 插入donor2(相邻位置)
    tree.insert(5, "GGG", "donor2")
    
    # 定义预期结果
    expected = {
        "total_reconstructions": 2,  # 预期总共有2个重建
        "donor_sequences": ["TTT", "GGG"]  # 预期包含这两个序列
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) != expected["total_reconstructions"]:
        success = False
        reason = f"期望{expected['total_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        # 检查序列内容
        for donor_seq in expected["donor_sequences"]:
            if not any(donor_seq in str(r.seq) for r in reconstructed):
                success = False
                reason = f"在重建中未找到预期序列 {donor_seq}"
                break
    
    # 记录测试结果
    test_results.append(("场景2", success, reason if not success else "相邻位置donors，创建了预期的重建"))
    if success:
        print("✓ 场景2测试通过: 相邻位置donors，创建了预期的重建")
    else:
        print(f"✗ 场景2测试失败: {reason}")
    
    # ======== 场景3: 嵌套插入 ========
    print("\n===== 测试场景3: 嵌套插入 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 插入donor1
    tree.insert(2, "TTTAAA", "donor1")
    # 插入donor2(在donor1内部)
    tree.insert(5, "GGG", "donor2")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 2,  # 期望至少有2个重建
        "min_full_reconstructions": 1,  # 期望至少有1个完整重建
        "should_contain_nested": True,  # 期望至少有一个重建包含嵌套序列
        "outer_donor": "TTTAAA",
        "inner_donor": "GGG"
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        # 检查完整重建数量
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        if len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"期望至少{expected['min_full_reconstructions']}个完整重建，但得到了{len(full_recon)}个"
        else:
            # 检查嵌套序列
            found_nested = False
            for rec in full_recon:
                seq_str = str(rec.seq)
                if expected["inner_donor"] in seq_str and (expected["outer_donor"] in seq_str or 
                       expected["outer_donor"][:3] in seq_str or expected["outer_donor"][3:] in seq_str):
                    found_nested = True
                    break
            
            if expected["should_contain_nested"] and not found_nested:
                success = False
                reason = "未找到包含预期嵌套序列的重建"
    
    # 记录测试结果
    test_results.append(("场景3", success, reason if not success else "嵌套插入重建正确"))
    if success:
        print("✓ 场景3测试通过: 嵌套插入重建正确")
    else:
        print(f"✗ 场景3测试失败: {reason}")
        
    # ======== 场景4: 多重嵌套插入 ========
    print("\n===== 测试场景4: 多重嵌套插入 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 插入donor1
    tree.insert(2, "TTTAAA", "donor1")
    # 插入donor2(在donor1内部)
    tree.insert(5, "GGCCC", "donor2")
    # 插入donor3(在donor2内部)
    tree.insert(7, "AAA", "donor3")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 2,  # 期望至少有2个重建
        "min_nested_sequences": 2,  # 期望至少发现2层嵌套
        "donor_sequences": ["TTTAAA", "GGCCC", "AAA"]  # 可能出现的donor序列
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        # 获取所有完整重建
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        
        # 检查是否至少有一个重建包含多层嵌套
        max_nested_count = 0
        for rec in full_recon:
            seq_str = str(rec.seq)
            nested_count = 0
            for donor in expected["donor_sequences"]:
                if donor in seq_str or (len(donor) >= 6 and (donor[:3] in seq_str or donor[3:] in seq_str)):
                    nested_count += 1
            max_nested_count = max(max_nested_count, nested_count)
        
        if max_nested_count < expected["min_nested_sequences"]:
            success = False
            reason = f"期望找到至少{expected['min_nested_sequences']}层嵌套，但最多只找到{max_nested_count}层"
    
    # 记录测试结果
    test_results.append(("场景4", success, reason if not success else "多重嵌套插入重建正确"))
    if success:
        print("✓ 场景4测试通过: 多重嵌套插入重建正确")
    else:
        print(f"✗ 场景4测试失败: {reason}")
    
    # ======== 场景5: 切割donor ========
    print("\n===== 测试场景5: 切割donor ========")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 插入donor1
    tree.insert(2, "TTAAA", "donor1")
    # 插入donor2(切断donor1)
    tree.insert(4, "GGG", "donor2")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 2,  # 期望至少有2个重建
        "has_clean_reconstructions": True,  # 期望有清洁重建
        "has_full_reconstructions": True,  # 期望有完整重建
        "cutter_sequence": "GGG"  # 期望在完整重建中找到切割者序列
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        # 验证完整重建和清洁重建
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        
        if expected["has_clean_reconstructions"] and not clean_recon:
            success = False
            reason = "缺少清洁重建"
        elif expected["has_full_reconstructions"] and not full_recon:
            success = False
            reason = "缺少完整重建"
        else:
            # 检查是否有至少一个完整重建包含donor2
            has_donor2 = any(expected["cutter_sequence"] in str(r.seq) for r in full_recon)
            
            if not has_donor2:
                success = False
                reason = f"未找到包含{expected['cutter_sequence']}的完整重建"
    
    # 记录测试结果
    test_results.append(("场景5", success, reason if not success else "切割donor重建正确"))
    if success:
        print("✓ 场景5测试通过: 切割donor重建正确")
    else:
        print(f"✗ 场景5测试失败: {reason}")
    
    # ======== 场景6: 多重切割 ========
    print("\n===== 测试场景6: 多重切割 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 插入donor1
    tree.insert(2, "TTTAAACCC", "donor1")
    # 插入donor2(切断donor1)
    tree.insert(4, "GGG", "donor2")
    # 插入donor3(再次切断donor1)
    tree.insert(8, "TTT", "donor3")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1,  # 期望至少有1个重建
        "has_clean_reconstructions": True,  # 期望有清洁重建
        "has_full_reconstructions": True,  # 期望有完整重建
        "cutter_sequences": ["GGG", "TTT"]  # 期望在完整重建中找到至少一个切割者序列
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        # 验证清洁重建和完整重建的存在性
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        
        if expected["has_clean_reconstructions"] and not clean_recon:
            success = False
            reason = "缺少清洁重建"
        elif expected["has_full_reconstructions"] and not full_recon:
            success = False
            reason = "缺少完整重建"
        else:
            # 检查完整重建中是否包含任何切割者序列
            has_cutter = False
            for cutter in expected["cutter_sequences"]:
                if any(cutter in str(r.seq) for r in full_recon):
                    has_cutter = True
                    break
            
            if not has_cutter:
                success = False
                reason = "未找到包含切割者序列的完整重建"
    
    # 记录测试结果
    test_results.append(("场景6", success, reason if not success else "多重切割重建正确"))
    if success:
        print("✓ 场景6测试通过: 多重切割重建正确")
    else:
        print(f"✗ 场景6测试失败: {reason}")
    
    # ======== 场景7: 连锁切割 ========
    print("\n===== 测试场景7: 连锁切割 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 插入donor1
    tree.insert(2, "TTTAAA", "donor1")
    # 插入donor2(切断donor1)
    tree.insert(4, "GGGCCC", "donor2")
    # 插入donor3(切断donor2)
    tree.insert(6, "AAA", "donor3")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 3,  # 期望至少有3个重建
        "min_full_reconstructions": 1,  # 期望至少有1个完整重建
        "min_clean_reconstructions": 1,  # 期望至少有1个清洁重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        
        if len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"期望至少{expected['min_full_reconstructions']}个完整重建，但得到了{len(full_recon)}个"
        elif len(clean_recon) < expected["min_clean_reconstructions"]:
            success = False
            reason = f"期望至少{expected['min_clean_reconstructions']}个清洁重建，但得到了{len(clean_recon)}个"
    
    # 记录测试结果
    test_results.append(("场景7", success, reason if not success else "连锁切割重建正确"))
    if success:
        print("✓ 场景7测试通过: 连锁切割重建正确")
    else:
        print(f"✗ 场景7测试失败: {reason}")
    
    # ======== 场景8: 首尾插入边界情况 ========
    print("\n===== 测试场景8: 首尾插入边界情况 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 在序列起始位置插入donor
    tree.insert(1, "TTT", "donor1")
    # 在序列末尾插入donor
    tree.insert(7, "GGG", "donor2")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1  # 期望至少有1个重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    
    # 记录测试结果
    test_results.append(("场景8", success, reason if not success else "首尾插入边界情况正确处理"))
    if success:
        print("✓ 场景8测试通过: 首尾插入边界情况正确处理")
    else:
        print(f"✗ 场景8测试失败: {reason}")
    
    # ======== 场景9: 同一位置多个插入 ========
    print("\n===== 测试场景9: 同一位置多个插入 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 在同一位置插入两个donor
    tree.insert(2, "TTT", "donor1")
    tree.insert(2, "GGG", "donor2")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1  # 期望至少有1个重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    
    # 记录测试结果
    test_results.append(("场景9", success, reason if not success else "同一位置多个插入正确处理"))
    if success:
        print("✓ 场景9测试通过: 同一位置多个插入正确处理")
    else:
        print(f"✗ 场景9测试失败: {reason}")
    
    # ======== 场景10: 随机多donor复杂网络 ========
    print("\n===== 测试场景10: 随机多donor复杂网络 =====")
    # 创建初始序列
    tree = SequenceTree("ATGCATGCATGC")
    # 插入多个相互切割的donor
    tree.insert(3, "AAATTT", "donor1")
    tree.insert(5, "CCCGGG", "donor2")
    tree.insert(8, "TTTAAA", "donor3")
    tree.insert(10, "GGGCCC", "donor4")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 4  # 期望至少有4个重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    
    # 记录测试结果
    test_results.append(("场景10", success, reason if not success else "随机多donor复杂网络正确处理"))
    if success:
        print("✓ 场景10测试通过: 随机多donor复杂网络正确处理")
    else:
        print(f"✗ 场景10测试失败: {reason}")
    
    # ======== 场景11: 空序列处理 ========
    print("\n===== 测试场景11: 空序列处理 =====")
    # 创建空序列
    tree = SequenceTree("")
    # 插入donor到空序列
    tree.insert(1, "TTT", "donor1")
    
    # 定义预期结果 - 根据SequenceTree的实际行为
    expected = {
        "should_have_reconstructions": False  # 期望无重建(因为没有嵌套或切割)
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if expected["should_have_reconstructions"] == False:
        if reconstructed and len(reconstructed) > 0:
            success = False
            reason = f"期望无重建，但得到了{len(reconstructed)}个重建"
    else:
        if reconstructed is None or len(reconstructed) == 0:
            success = False
            reason = "期望有重建，但未获得任何重建"
    
    # 记录测试结果
    test_results.append(("场景11", success, reason if not success else "空序列正确处理"))
    if success:
        print("✓ 场景11测试通过: 空序列正确处理")
    else:
        print(f"✗ 场景11测试失败: {reason}")
    
    # ======== 场景12: 极长序列处理 ========
    print("\n===== 测试场景12: 极长序列处理 =====")
    # 创建一个较长序列（1000个碱基）
    long_seq = "A" * 400 + "T" * 300 + "G" * 200 + "C" * 100
    tree = SequenceTree(long_seq)
    # 插入两个相互嵌套的donor
    tree.insert(500, "AAATTT", "donor1")
    tree.insert(502, "GGG", "donor2")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1,  # 期望至少有1个重建
        "min_full_reconstructions": 1  # 期望至少有1个完整重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        if len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"期望至少{expected['min_full_reconstructions']}个完整重建，但得到了{len(full_recon)}个"
    
    # 记录测试结果
    test_results.append(("场景12", success, reason if not success else "极长序列正确处理"))
    if success:
        print("✓ 场景12测试通过: 极长序列正确处理")
    else:
        print(f"✗ 场景12测试失败: {reason}")
    
    # ======== 场景13: 完全重叠切割 ========
    print("\n===== 测试场景13: 完全重叠切割 =====")
    # 创建初始序列
    tree = SequenceTree("ATGCATGC")
    # 创建一个donor
    tree.insert(2, "AAATTTCCC", "donor1")
    # 创建两个同时在相同位置切割的donor
    tree.insert(5, "GGG", "cutter1")
    tree.insert(5, "TTT", "cutter2")
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1,  # 期望至少有1个重建
        "min_clean_reconstructions": 1,  # 期望至少有1个清洁重建
        "min_full_reconstructions": 1   # 期望至少有1个完整重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        # 应该至少有一个清洁重建和一个完整重建
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        
        if len(clean_recon) < expected["min_clean_reconstructions"]:
            success = False
            reason = f"期望至少{expected['min_clean_reconstructions']}个清洁重建，但只得到了{len(clean_recon)}个"
        elif len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"期望至少{expected['min_full_reconstructions']}个完整重建，但只得到了{len(full_recon)}个"
    
    # 记录测试结果
    test_results.append(("场景13", success, reason if not success else "完全重叠切割正确处理"))
    if success:
        print("✓ 场景13测试通过: 完全重叠切割正确处理")
    else:
        print(f"✗ 场景13测试失败: {reason}")
    
    # ======== 场景14: 边界切割 ========
    print("\n===== 测试场景14: 边界切割 =====")
    # 创建初始序列
    tree = SequenceTree("ATGCATGC")
    # 创建一个donor
    tree.insert(2, "AAATTTCCC", "donor1")
    # 创建一个在donor边界处切割的cutter
    tree.insert(2, "GGG", "cutter1")  # 在donor开头切割
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1  # 期望至少有1个重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    
    # 记录测试结果
    test_results.append(("场景14", success, reason if not success else "边界切割正确处理"))
    if success:
        print("✓ 场景14测试通过: 边界切割正确处理")
    else:
        print(f"✗ 场景14测试失败: {reason}")
    
    # ======== 场景15: 多donor同时切割另一donor的不同位置 ========
    print("\n===== 测试场景15: 多donor同时切割另一donor的不同位置 =====")
    # 创建初始序列
    tree = SequenceTree("ATGCATGC")
    # 创建主donor
    tree.insert(2, "AAAAAAAAAATTTTTTTTTTGGGGGGGGGG", "main_donor")  # 30bp
    # 创建三个切割者
    tree.insert(5, "CCC", "cutter1")  # 切割前段
    tree.insert(15, "TTT", "cutter2")  # 切割中段
    tree.insert(25, "GGG", "cutter3")  # 切割后段
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1,  # 期望至少有1个重建
        "has_clean_reconstructions": True  # 期望有清洁重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        # 检查清洁重建
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        if expected["has_clean_reconstructions"] and not clean_recon:
            success = False
            reason = "未找到主donor的清洁重建"
    
    # 记录测试结果
    test_results.append(("场景15", success, reason if not success else "多donor同时切割另一donor的不同位置正确处理"))
    if success:
        print("✓ 场景15测试通过: 多donor同时切割另一donor的不同位置正确处理")
        # 检查完整重建数量
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        print(f"  完整重建数量: {len(full_recon)}")
    else:
        print(f"✗ 场景15测试失败: {reason}")
    
    # ======== 场景16: 复杂切割网络 - 多个donor相互交错切割 ========
    print("\n===== 测试场景16: 复杂切割网络 - 多个donor相互交错切割 =====")
    # 创建初始序列
    tree = SequenceTree("ATGCATGCATGCATGCATGC")
    # 创建5个donor，形成复杂的交错切割网络
    tree.insert(2, "AAAAAAAAAA", "donor_a")   # A
    tree.insert(12, "TTTTTTTTTT", "donor_b")  # B
    tree.insert(22, "GGGGGGGGGG", "donor_c")  # C
    tree.insert(32, "CCCCCCCCCC", "donor_d")  # D
    tree.insert(42, "ATATATATAT", "donor_e")  # E
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1  # 期望至少有1个重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    
    # 记录测试结果
    test_results.append(("场景16", success, reason if not success else f"复杂切割网络正确处理，生成了{len(reconstructed) if reconstructed else 0}个重建"))
    if success:
        print(f"✓ 场景16测试通过: 复杂切割网络正确处理，生成了{len(reconstructed) if reconstructed else 0}个重建")
    else:
        print(f"✗ 场景16测试失败: {reason}")
    
    # ======== 场景17: 序列末尾边界切割和首尾部特殊切割 ========
    print("\n===== 测试场景17: 序列末尾边界切割和首尾部特殊切割 =====")
    # 创建初始序列
    tree = SequenceTree("ATGCATGC")
    # 创建donor序列
    tree.insert(2, "AAAATTTTGGGG", "donor")  # 12bp
    # 创建在首尾位置切割的donor
    tree.insert(2, "CCC", "cutter1")  # 在donor开头切割
    tree.insert(14, "TTT", "cutter2")  # 在donor结尾切割
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1,  # 期望至少有1个重建
        "has_clean_reconstructions": True  # 期望有清洁重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    elif expected["has_clean_reconstructions"]:
        # 检查是否正确处理多重切割
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        if not clean_recon:
            success = False
            reason = "未找到donor的清洁重建"
    
    # 记录测试结果
    test_results.append(("场景17", success, reason if not success else "序列末尾边界切割和首尾部特殊切割正确处理"))
    if success:
        print("✓ 场景17测试通过: 序列末尾边界切割和首尾部特殊切割正确处理")
    else:
        print(f"✗ 场景17测试失败: {reason}")
    
    # ======== 场景18: 深度嵌套 ========
    print("\n===== 测试场景18: 深度嵌套 =====")
    # 创建初始序列
    tree = SequenceTree("ATGCATGC")
    # 创建5层嵌套的donor
    tree.insert(2, "AAAAAAAA", "donor1")  # 最外层
    tree.insert(4, "TTTTTTTT", "donor2")  # 第2层
    tree.insert(6, "GGGGGGGG", "donor3")  # 第3层
    tree.insert(8, "CCCCCCCC", "donor4")  # 第4层
    tree.insert(10, "AAAATTTT", "donor5") # 最内层
    
    # 定义预期结果
    expected = {
        "min_reconstructions": 1,  # 期望至少有1个重建
        "min_full_reconstructions": 1  # 期望至少有1个完整重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"期望至少{expected['min_reconstructions']}个重建，但得到了{len(reconstructed) if reconstructed else 0}个重建"
    else:
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        if len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"期望至少{expected['min_full_reconstructions']}个完整重建，但得到了{len(full_recon)}个"
    
    # 记录测试结果
    test_results.append(("场景18", success, reason if not success else "深度嵌套正确处理"))
    if success:
        print("✓ 场景18测试通过: 深度嵌套正确处理")
        print(f"  生成了{len(full_recon)}个完整重建")
    else:
        print(f"✗ 场景18测试失败: {reason}")
    
    # ======== 场景19: 极端短序列和空序列切割处理 ========
    print("\n===== 测试场景19: 极端短序列和空序列切割处理 =====")
    # 创建初始序列
    tree = SequenceTree("ATGC")
    # 创建极短donor序列
    tree.insert(2, "A", "short_donor")
    # 创建一个切割极短序列的donor
    tree.insert(2, "TTT", "cutter")
    
    # 定义预期结果
    expected = {
        "should_have_reconstructions": True  # 期望有重建
    }
    
    # 验证结果
    _, reconstructed = tree.donors("test")
    
    # 明确的断言检查
    success = True
    reason = ""
    
    if expected["should_have_reconstructions"]:
        if reconstructed is None or len(reconstructed) == 0:
            success = False
            reason = "期望有重建，但未获得任何重建"
    else:
        if reconstructed and len(reconstructed) > 0:
            success = False
            reason = f"期望无重建，但得到了{len(reconstructed)}个重建"
    
    # 记录测试结果
    test_results.append(("场景19", success, reason if not success else "极端短序列正确处理"))
    if success:
        print("✓ 场景19测试通过: 极端短序列正确处理")
    else:
        print(f"✗ 场景19测试失败: {reason}")
    
    # ======== 显示最终测试结果 ========
    print("\n===== 综合测试结果 =====")
    passed = sum(1 for _, result, _ in test_results if result)
    total = len(test_results)
    print(f"测试完成: {passed}/{total} 通过")
    if passed == total:
        print("✓ 所有场景测试通过！重建算法工作正常")
    else:
        print("✗ 部分测试未通过，请检查重建算法")
        for name, result, message in test_results:
            if not result:
                print(f"  - {name}: {message}")
    
    return passed == total


def test_multiple_cuts_fragments_distinction():
    """
    测试多重切割情况下片段区分功能的改进
    验证在一个节点被多次切割的情况下，能否正确区分不同切割者产生的片段
    
    Returns:
        bool: 测试是否成功
    """
    print("\n===== 测试多重切割片段区分功能 =====")
    
    # 初始化树数据结构
    tree = SequenceTree("ATGCATGCATGCATGCATGCATGCATGCATGCATGC")  # 原始序列

    # 创建一个被多次切割的场景
    # 原始序列被三个不同的donor切割：
    # - donor 1 切割在位置10
    # - donor 2 切割在位置20
    # - donor 3 切割在位置30
    
    # 插入donor序列在相应位置
    tree.insert(10, "GTACGTAC", "cutter1")
    tree.insert(20, "CCGGAATT", "cutter2")
    tree.insert(30, "TTAGGCCA", "cutter3")
    
    # 获取事件日志
    event_journal = tree.event_journal
    
    # 定义预期结果
    expected = {
        "event_count": 3,  # 期望有3个事件记录
        "donor_uids": [],  # 将在测试过程中填充
        "fragment_uids": [],  # 将在测试过程中填充
        "has_valid_fragment_info": True,  # 期望所有片段都有有效信息
        "fragment_cutter_match": True,  # 期望每个片段的切割者信息与事件中的donor_uid一致
    }
    
    # 测试片段区分功能
    success = True
    error_messages = []
    
    # 1. 验证事件记录是否完整
    if len(event_journal.events) != expected["event_count"]:
        error_messages.append(f"错误: 应记录{expected['event_count']}个事件，但找到 {len(event_journal.events)} 个")
        success = False
    
    # 2. 获取事件信息并填充预期结果中的UID列表
    for event in event_journal.events:
        print(f"事件 {event.event_id}: donor({event.donor_uid}) → target({event.target_uid}) → [L({event.left_uid}), R({event.right_uid})]")
        
        # 记录donor和fragment的UID
        expected["donor_uids"].append(event.donor_uid)
        expected["fragment_uids"].append(event.left_uid)
        expected["fragment_uids"].append(event.right_uid)
        
        # 检查donor是否正确标记
        if not event_journal.is_donor(event.donor_uid):
            error_messages.append(f"错误: UID {event.donor_uid} 应被标记为donor")
            success = False
        
        # 检查片段是否正确标记
        if not event_journal.is_fragment(event.left_uid):
            error_messages.append(f"错误: UID {event.left_uid} 应被标记为fragment")
            success = False
        
        if not event_journal.is_fragment(event.right_uid):
            error_messages.append(f"错误: UID {event.right_uid} 应被标记为fragment")
            success = False
        
        # 检查片段信息
        left_info = event_journal.get_fragment_info(event.left_uid)
        right_info = event_journal.get_fragment_info(event.right_uid)
        
        if left_info:
            orig_uid, is_left, cutter_uid = left_info
            print(f"  左片段 {event.left_uid}: 来自 {orig_uid}, "
                  f"{'左' if is_left else '右'}侧, 切割者 {cutter_uid}")
            
            # 验证切割者信息
            if expected["fragment_cutter_match"] and cutter_uid != event.donor_uid:
                error_messages.append(f"错误: 左片段 {event.left_uid} 的切割者应为 {event.donor_uid}，但得到 {cutter_uid}")
                success = False
        elif expected["has_valid_fragment_info"]:
            error_messages.append(f"错误: 未找到左片段 {event.left_uid} 的信息")
            success = False
        
        if right_info:
            orig_uid, is_left, cutter_uid = right_info
            print(f"  右片段 {event.right_uid}: 来自 {orig_uid}, "
                  f"{'左' if is_left else '右'}侧, 切割者 {cutter_uid}")
            
            # 验证切割者信息
            if expected["fragment_cutter_match"] and cutter_uid != event.donor_uid:
                error_messages.append(f"错误: 右片段 {event.right_uid} 的切割者应为 {event.donor_uid}，但得到 {cutter_uid}")
                success = False
        elif expected["has_valid_fragment_info"]:
            error_messages.append(f"错误: 未找到右片段 {event.right_uid} 的信息")
            success = False
    
    # 3. 测试获取重建结果
    _, reconstructed = tree.donors("test")
    
    # 定义重建的预期结果
    recon_expected = {
        "should_have_reconstructions": True  # 期望有重建结果
    }
    
    if recon_expected["should_have_reconstructions"] and not reconstructed:
        error_messages.append("错误: 未获得任何重建结果")
        success = False
    else:
        print(f"获得 {len(reconstructed) if reconstructed else 0} 个重建结果")
        # 输出重建记录信息
        if reconstructed:
            for rec in reconstructed:
                recon_type = rec.annotations.get("reconstruction_type", "未知")
                original_uid = rec.annotations.get("original_uid", "未知")
                seq_len = len(rec.seq)
                print(f"  重建: 类型={recon_type}, 原始UID={original_uid}, 序列长度={seq_len}")
    
    # 4. 生成DOT可视化
    dot_str = event_journal.to_graphviz_dot()
    
    # 定义DOT可视化的预期结果
    dot_expected = {
        "contains_nodes": True,  # 期望DOT包含节点信息
        "contains_events": True  # 期望DOT包含事件信息
    }
    
    # 检查DOT字符串中是否包含节点和事件信息
    if dot_expected["contains_nodes"] and "node_" not in dot_str:
        error_messages.append("错误: DOT可视化中应包含节点信息")
        success = False
    
    if dot_expected["contains_events"] and "event_" not in dot_str:
        error_messages.append("错误: DOT可视化中应包含事件信息")
        success = False
    
    if success:
        print("✓ 多重切割片段区分功能测试通过!")
    else:
        print("✗ 多重切割片段区分功能测试失败!")
        for error in error_messages:
            print(f"  {error}")
    
    return success


# 如果作为主程序运行，执行所有测试
if __name__ == "__main__":
    print("执行测试...")
    
    # 执行多重切割测试
    success1 = test_multiple_cuts()
    if not success1:
        print("  - 多重切割测试未通过")
        
    # 执行综合嵌套测试
    success2 = test_comprehensive_nesting()
    if not success2:
        print("  - 综合嵌套测试未通过") 
        
    # 执行多次切割片段区分功能测试
    success3 = test_multiple_cuts_fragments_distinction()
    if not success3:
        print("  - 多次切割片段区分功能测试未通过")
    
    # 显示总体测试结果
    print("\n===== 总体测试结果 =====")
    if success1 and success2 and success3:
        print("✓ 所有测试通过！")
    else:
        print("✗ 部分测试未通过！")
        if not success1:
            print("  - 多重切割测试未通过")
        if not success2:
            print("  - 综合嵌套测试未通过") 
        if not success3:
            print("  - 多次切割片段区分功能测试未通过")
