#!/usr/bin/env python3
"""
RTファイルマージスクリプト
複数のRTファイルを統合してバックグラウンドファイルを作成
"""

import argparse
import pandas as pd
import numpy as np
import sys
from collections import defaultdict

def merge_rt_files(input_files, output_file):
    """
    複数のRTファイルをマージ
    
    Args:
        input_files (list): 入力RTファイルのリスト
        output_file (str): 出力マージ済みファイル
    """
    print(f"RTファイルマージ開始")
    print(f"入力ファイル数: {len(input_files)}")
    
    merged_data = defaultdict(lambda: defaultdict(list))
    
    # 各ファイルからデータ読み込み
    for i, rt_file in enumerate(input_files):
        print(f"読み込み中 ({i+1}/{len(input_files)}): {rt_file}")
        
        try:
            with open(rt_file, 'r') as f:
                header = f.readline().strip()
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        transcript_id = parts[0]
                        position = int(parts[1])
                        rtstop_count = int(parts[2])
                        base_coverage = int(parts[3])
                        rtstop_ratio = float(parts[4])
                        
                        merged_data[transcript_id][position].append({
                            'rtstop_count': rtstop_count,
                            'base_coverage': base_coverage,
                            'rtstop_ratio': rtstop_ratio
                        })
        
        except Exception as e:
            print(f"警告: ファイル読み込みエラー - {rt_file}: {e}")
            continue
    
    print(f"マージ対象転写産物数: {len(merged_data)}")
    
    # マージ済みデータをファイルに保存
    try:
        with open(output_file, 'w') as f:
            f.write("transcript_id\\tposition\\trtstop_count\\tbase_coverage\\trtstop_ratio\\n")
            
            for transcript_id in merged_data:
                for position in merged_data[transcript_id]:
                    data_points = merged_data[transcript_id][position]
                    
                    # 平均値計算
                    if data_points:
                        avg_rtstop_count = int(np.mean([d['rtstop_count'] for d in data_points]))
                        avg_base_coverage = int(np.mean([d['base_coverage'] for d in data_points]))
                        avg_rtstop_ratio = np.mean([d['rtstop_ratio'] for d in data_points])
                        
                        f.write(f"{transcript_id}\\t{position}\\t{avg_rtstop_count}\\t"
                               f"{avg_base_coverage}\\t{avg_rtstop_ratio:.6f}\\n")
        
        print(f"マージ完了: {output_file}")
        
        # 統計情報
        total_positions = sum(len(positions) for positions in merged_data.values())
        print(f"総ポジション数: {total_positions}")
        
    except Exception as e:
        print(f"エラー: RTファイルマージ中にエラーが発生しました - {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="複数のRTファイルをマージ")
    parser.add_argument("-i", "--input", required=True, 
                       help="入力RTファイル（:区切り）")
    parser.add_argument("-o", "--output", required=True, help="出力マージ済みファイル")
    
    args = parser.parse_args()
    
    # 入力ファイルリストを解析
    input_files = args.input.split(':')
    input_files = [f.strip() for f in input_files if f.strip()]
    
    if len(input_files) < 2:
        print("エラー: 2つ以上の入力ファイルが必要です")
        sys.exit(1)
    
    merge_rt_files(input_files, args.output)

if __name__ == "__main__":
    main()