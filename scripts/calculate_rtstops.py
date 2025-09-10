#!/usr/bin/env python3
"""
RTstop計算スクリプト
逆転写停止サイトをカウント
"""

import argparse
import pysam
import pandas as pd
from collections import defaultdict
import sys

def calculate_rtstops(sam_file, rpkm_file, output_file, min_rpkm=1):
    """
    SAMファイルから逆転写停止サイト（RTstops）を計算
    
    Args:
        sam_file (str): 入力SAMファイルパス
        rpkm_file (str): RPKMファイルパス
        output_file (str): 出力RTstopファイルパス
        min_rpkm (float): 最小RPKM閾値
    """
    print(f"RTstop計算開始")
    print(f"SAMファイル: {sam_file}")
    print(f"RPKMファイル: {rpkm_file}")
    print(f"最小RPKM閾値: {min_rpkm}")
    
    # RPKM情報読み込み
    try:
        rpkm_df = pd.read_csv(rpkm_file, sep='\t')
        high_expr_transcripts = set(rpkm_df[rpkm_df['rpkm'] >= min_rpkm]['transcript_id'])
        print(f"高発現転写産物数 (RPKM >= {min_rpkm}): {len(high_expr_transcripts)}")
    except Exception as e:
        print(f"警告: RPKMファイル読み込みエラー - {e}")
        high_expr_transcripts = set()
    
    # RTstopカウント
    rtstop_counts = defaultdict(lambda: defaultdict(int))
    base_coverage = defaultdict(lambda: defaultdict(int))
    
    try:
        with pysam.AlignmentFile(sam_file, "r") as samfile:
            for read in samfile:
                if not read.is_unmapped:
                    transcript_name = read.reference_name
                    
                    # 高発現転写産物のみ処理（RPKMファイルがある場合）
                    if high_expr_transcripts and transcript_name not in high_expr_transcripts:
                        continue
                    
                    # リードの5'端をRTstopサイトとして記録
                    if not read.is_reverse:
                        rtstop_position = read.reference_start
                    else:
                        rtstop_position = read.reference_end - 1
                    
                    rtstop_counts[transcript_name][rtstop_position] += 1
                    
                    # カバレッジ計算
                    for pos in range(read.reference_start, read.reference_end):
                        base_coverage[transcript_name][pos] += 1
        
        print(f"RTstop解析完了: {len(rtstop_counts)} 転写産物")
        
        # RTstop結果をファイルに保存
        with open(output_file, 'w') as f:
            f.write("transcript_id\\tposition\\trtstop_count\\tbase_coverage\\trtstop_ratio\\n")
            
            for transcript_id in rtstop_counts:
                for position in rtstop_counts[transcript_id]:
                    rtstop_count = rtstop_counts[transcript_id][position]
                    coverage = base_coverage[transcript_id].get(position, 0)
                    
                    if coverage > 0:
                        rtstop_ratio = rtstop_count / coverage
                    else:
                        rtstop_ratio = 0.0
                    
                    f.write(f"{transcript_id}\\t{position}\\t{rtstop_count}\\t{coverage}\\t{rtstop_ratio:.6f}\\n")
        
        print(f"RTstop結果保存: {output_file}")
        
        # サマリー統計
        total_rtstops = sum(sum(counts.values()) for counts in rtstop_counts.values())
        print(f"総RTstop数: {total_rtstops}")
        
    except Exception as e:
        print(f"エラー: RTstop計算中にエラーが発生しました - {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="SAMファイルからRTstopを計算")
    parser.add_argument("-i", "--input", required=True, help="入力SAMファイル")
    parser.add_argument("-r", "--rpkm", required=True, help="RPKMファイル")
    parser.add_argument("-o", "--output", required=True, help="出力RTstopファイル")
    parser.add_argument("-c", "--min_rpkm", type=float, default=1, help="最小RPKM閾値")
    
    args = parser.parse_args()
    
    calculate_rtstops(args.input, args.rpkm, args.output, args.min_rpkm)

if __name__ == "__main__":
    main()