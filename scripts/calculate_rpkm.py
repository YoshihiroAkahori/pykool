#!/usr/bin/env python3
"""
RPKM計算スクリプト
SAMファイルから転写産物の発現量（RPKM）を計算
"""

import argparse
import pysam
import pandas as pd
from collections import defaultdict
import sys

def calculate_rpkm(sam_file, output_file):
    """
    SAMファイルからRPKM値を計算
    
    Args:
        sam_file (str): 入力SAMファイルパス
        output_file (str): 出力RPKMファイルパス
    """
    print(f"SAMファイル読み込み: {sam_file}")
    
    # 転写産物ごとのリード数とリード長をカウント
    transcript_counts = defaultdict(int)
    transcript_lengths = defaultdict(int)
    total_reads = 0
    
    try:
        with pysam.AlignmentFile(sam_file, "r") as samfile:
            for read in samfile:
                if not read.is_unmapped:
                    transcript_name = read.reference_name
                    transcript_counts[transcript_name] += 1
                    if transcript_lengths[transcript_name] == 0:
                        transcript_lengths[transcript_name] = read.reference_length
                    total_reads += 1
    
        print(f"総リード数: {total_reads}")
        print(f"マップされた転写産物数: {len(transcript_counts)}")
        
        # RPKM計算
        rpkm_data = []
        for transcript in transcript_counts:
            read_count = transcript_counts[transcript]
            length_kb = transcript_lengths[transcript] / 1000.0
            reads_per_million = total_reads / 1000000.0
            
            if length_kb > 0 and reads_per_million > 0:
                rpkm = read_count / (length_kb * reads_per_million)
            else:
                rpkm = 0.0
            
            rpkm_data.append({
                'transcript_id': transcript,
                'read_count': read_count,
                'length_bp': transcript_lengths[transcript],
                'rpkm': rpkm
            })
        
        # DataFrameに変換してソート
        rpkm_df = pd.DataFrame(rpkm_data)
        rpkm_df = rpkm_df.sort_values('rpkm', ascending=False)
        
        # ファイルに保存
        rpkm_df.to_csv(output_file, sep='\t', index=False)
        print(f"RPKM結果保存: {output_file}")
        print(f"上位5転写産物のRPKM:")
        print(rpkm_df.head())
        
    except Exception as e:
        print(f"エラー: RPKM計算中にエラーが発生しました - {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="SAMファイルからRPKM値を計算")
    parser.add_argument("-i", "--input", required=True, help="入力SAMファイル")
    parser.add_argument("-o", "--output", required=True, help="出力RPKMファイル")
    
    args = parser.parse_args()
    
    calculate_rpkm(args.input, args.output)

if __name__ == "__main__":
    main()