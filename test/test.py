#!/usr/bin/env python3

import argparse
import logging
from Bio import SeqIO
import re
from SingleSequenceMatcher import SingleSequenceMatcher

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def extract_chr_id(seq_id):
    """
    从序列ID中提取染色体ID
    
    Args:
        seq_id (str): 序列ID (如 "ChrC_ins10_tsd9")
        
    Returns:
        str: 提取的染色体ID (如 "ChrC")
    """
    parts = seq_id.split('_')
    if not parts:
        return seq_id  # 如果没有下划线，返回整个ID
    return parts[0]

def modify_sequence_ids(input_file):
    """
    读取输入文件中的序列并修改ID只保留染色体部分
    
    Args:
        input_file (str): 输入FASTA文件路径
        
    Returns:
        dict: 染色体ID到序列的映射字典
    """
    # 使用helper function和列表推导式创建(chr_id, sequence)元组列表
    records = [(extract_chr_id(record.id), record.seq) for record in SeqIO.parse(input_file, "fasta")]
    
    # 转换为字典
    sequences = {}
    for chr_id, seq in records:
        if chr_id in sequences:
            logger.warning(f"Found multiple sequences for chromosome {chr_id}. Overwriting.")
        sequences[chr_id] = seq
    
    logger.info(f"Loaded {len(sequences)} unique chromosome sequences")
    return sequences

def parse_donor_id(donor_id):
    """
    解析donor ID提取染色体、起始和结束位置
    
    Args:
        donor_id (str): Donor ID (如 "ChrC_9177_9472-+-295")
        
    Returns:
        tuple: (chromosome, start, end) 或解析失败时返回None
    """
    match = re.match(r'([^_]+)_(\d+)_(\d+).*', donor_id)
    if not match:
        return None
    
    chr_id, start, end = match.groups()
    return chr_id, int(start), int(end)

def compare_donor_sequences(donor_file, sequences, verbose=False):
    """
    比较donor文件中的序列与main序列的切片
    
    Args:
        donor_file (str): Donor FASTA文件路径
        sequences (dict): 染色体ID到序列的映射字典
        verbose (bool): 是否输出详细信息
        
    Returns:
        float: 成功匹配的百分比
    """
    total_donors = 0
    successful_matches = 0
    
    for record in SeqIO.parse(donor_file, "fasta"):
        total_donors += 1
        
        # 解析donor ID
        parse_result = parse_donor_id(record.id)
        if parse_result is None:
            logger.warning(f"Could not parse donor ID format: {record.id}")
            continue
            
        chr_id, start, end = parse_result
        
        # 检查染色体是否存在于main序列中
        if chr_id not in sequences:
            logger.warning(f"Chromosome {chr_id} not found in main sequences")
            continue
            
        # 检查坐标是否有效
        if start <= 0 or end > len(sequences[chr_id]):
            logger.warning(f"Invalid coordinates for {record.id} (start={start}, end={end}, sequence length={len(sequences[chr_id])})")
            continue
            
        # 提取main序列的切片 (转换为0-based索引)
        main_seq = sequences[chr_id][start-1:end]  
        
        # 比较序列(转为大写进行比较，避免大小写差异)
        if str(main_seq).upper() == str(record.seq).upper():
            successful_matches += 1
            if verbose:
                logger.info(f"Match found for {record.id}")
        elif verbose:
            logger.info(f"No match for {record.id}")
            logger.debug(f"Main sequence slice: {main_seq}")
            logger.debug(f"Donor sequence: {record.seq}")
    
    if total_donors == 0:
        logger.warning("No donor sequences processed")
        return 0.0
    
    return successful_matches / total_donors

def find_donors_in_library(donor_file, library_file, min_identity=95.0, verbose=False):
    """
    在donor库中查找donors，要求最小一致性百分比
    
    Args:
        donor_file (str): Donor FASTA文件路径
        library_file (str): Donor库FASTA文件路径
        min_identity (float): 成功匹配的最小一致性百分比
        verbose (bool): 是否输出详细信息
        
    Returns:
        float: 成功匹配的百分比
    """
    # 初始化matcher并设置最小一致性阈值为用户指定的值
    matcher = SingleSequenceMatcher(min_identity=min_identity)
    
    total_donors = 0
    successful_matches = 0
    
    for record in SeqIO.parse(donor_file, "fasta"):
        total_donors += 1
        query_sequence = str(record.seq)
        
        # 跳过空序列
        if not query_sequence:
            logger.warning(f"Empty sequence found for {record.id}")
            continue
        
        # 在库中查找最佳匹配
        best_ref_id, best_identity, best_length, _ = matcher.find_best_match(
            query_sequence=query_sequence,
            reference_file=library_file
        )
        
        if best_ref_id is not None:
            successful_matches += 1
            if verbose:
                logger.info(f"Match found for {record.id}: {best_ref_id} (Identity: {best_identity:.2f}%, Length: {best_length})")
        elif verbose:
            logger.info(f"No match found for {record.id}")
    
    if total_donors == 0:
        logger.warning("No donor sequences processed")
        return 0.0
    
    return successful_matches / total_donors

def main():
    parser = argparse.ArgumentParser(description='Compare sequences and donors')
    parser.add_argument('-s', '--sequences', required=True, help='Path to main sequences FASTA file')
    parser.add_argument('-d', '--donors', required=True, help='Path to donor sequences FASTA file')
    parser.add_argument('-l', '--library', required=True, help='Path to donor library FASTA file')
    parser.add_argument('--min-identity', type=float, default=95.0, 
                        help='Minimum identity percentage for library matching (default: 95.0)')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
    # 设置日志级别
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    
    # 步骤1: 读取并修改main序列
    sequences = modify_sequence_ids(args.sequences)
    logger.info(f"Read {len(sequences)} main sequences")
    
    # 步骤2: 比较donor序列与main序列
    donor_match_perc = compare_donor_sequences(args.donors, sequences, verbose=args.verbose)
    print(f"Donor sequence matching: {donor_match_perc:.2%}")
    
    # 步骤3: 在库中查找donors
    library_match_perc = find_donors_in_library(args.donors, args.library, args.min_identity, verbose=args.verbose)
    print(f"Donor library matching: {library_match_perc:.2%}")
    
    return 0

if __name__ == "__main__":
    exit(main())
