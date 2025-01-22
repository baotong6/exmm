'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2025-01-06 17:52:06
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2025-01-07 15:42:25
FilePath: /code/match_eSASS/Bayesian_prob.py
Description: 

Copyright (c) 2025 by baotong, All Rights Reserved. 
'''
import math

def calculate_posterior_probabilities(rho, r, p_true_match):
    # Step 1: 计算 P(Match | False Match)
    pi_r2 = math.pi * r**2
    p_match_given_false_match = 1 - math.exp(-rho * pi_r2)
    
    # Step 2: 计算 P(Match)
    p_match_given_true_match = 1  # 真实匹配的概率接近 1
    p_false_match = 1 - p_true_match
    p_match = (p_match_given_true_match * p_true_match) + (p_match_given_false_match * p_false_match)
    
    # Step 3: 使用贝叶斯公式计算 P(True Match | Match)
    p_true_match_given_match = (p_match_given_true_match * p_true_match) / p_match
    
    # Step 4: 计算 P(False Match | Match)
    p_false_match_given_match = 1 - p_true_match_given_match
    
    return p_true_match_given_match, p_false_match_given_match


# 定义输入参数
rho = 0.035967052  # sources per square arcsecond
r = 4  # arcsecond
p_true_match = 0.7  # 假设 P(True Match) 为 0.5

# 计算后验概率
p_true_match_given_match, p_false_match_given_match = calculate_posterior_probabilities(rho, r, p_true_match)

# 输出结果
print(f"P(True Match | Match) = {p_true_match_given_match:.4f}")
print(f"P(False Match | Match) = {p_false_match_given_match:.4f}")
