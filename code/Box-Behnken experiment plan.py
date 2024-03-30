# Box-Behnken设计
import pyDOE2
# 创建一个包含4个因素的Box-Behnken设计
design = pyDOE2.bbdesign(4)
print(design)
