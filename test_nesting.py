#!/usr/bin/env python3

from RandSeqInsert import DonorNestingGraph

# 创建测试图实例
test_graph = DonorNestingGraph()

# 测试comprehensive_nesting方法
print("======= 测试comprehensive_nesting方法 =======")
test_graph.test_comprehensive_nesting()

print("\n测试完成！") 
