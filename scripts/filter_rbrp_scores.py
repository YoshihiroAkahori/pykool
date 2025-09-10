#!/usr/bin/env python3
"""
RBRPスコアフィルタリングスクリプト
低品質なRBRPスコアを除去
"""

import argparse
import pandas as pd
import numpy as np
import sys

def filter_rbrp_scores(input_file, output_file, min_depth=200, skip_start=5, skip_end=30):
    """
    RBRPスコアをフィルタリング
    
    Args:
        input_file (str): 入力RBRPスコアファイル
        output_file (str): 出力フィルタリング済みファイル
        min_depth (int): 最小シーケンス深度
        skip_start (int): 転写産物開始部分をスキップ
        skip_end (int): 転写産物終端部分をスキップ
    """
    print(f"RBRPスコアフィルタリング開始")
    print(f"入力ファイル: {input_file}")
    print(f"最小シーケンス深度: {min_depth}")
    print(f"開始スキップ: {skip_start} nt")
    print(f"終端スキップ: {skip_end} nt")
    
    try:
        # RBRPスコアファイル読み込み
        df = pd.read_csv(input_file, sep='\t')
        original_count = len(df)
        print(f"元のデータポイント数: {original_count}")
        
        # フィルタリング条件
        # 1. 最小シーケンス深度
        df_filtered = df[df['base_coverage'] >= min_depth].copy()
        depth_filtered = len(df_filtered)
        print(f"深度フィルタリング後: {depth_filtered} ({depth_filtered/original_count*100:.1f}%)")
        
        # 2. 転写産物の端部除去
        transcript_lengths = df_filtered.groupby('transcript_id')['position'].max()
        
        def is_valid_position(row):
            transcript_id = row['transcript_id']
            position = row['position']
            max_pos = transcript_lengths[transcript_id]
            
            # 開始部分と終端部分をスキップ
            return (position >= skip_start) and (position <= max_pos - skip_end)
        
        df_filtered = df_filtered[df_filtered.apply(is_valid_position, axis=1)].copy()
        position_filtered = len(df_filtered)
        print(f"位置フィルタリング後: {position_filtered} ({position_filtered/original_count*100:.1f}%)")
        
        # 3. 統計情報計算
        if len(df_filtered) > 0:
            print(f"\\nフィルタリング済みRBRPスコア統計:")
            print(f"  平均: {df_filtered['rbrp_score'].mean():.4f}")
            print(f"  中央値: {df_filtered['rbrp_score'].median():.4f}")
            print(f"  標準偏差: {df_filtered['rbrp_score'].std():.4f}")
            print(f"  最大値: {df_filtered['rbrp_score'].max():.4f}")
            print(f"  最小値: {df_filtered['rbrp_score'].min():.4f}")
            
            # 高スコア領域の統計
            high_score_threshold = 0.12
            high_scores = df_filtered[df_filtered['rbrp_score'] >= high_score_threshold]
            print(f"  高スコア領域 (>= {high_score_threshold}): {len(high_scores)} ポイント")
            
            if len(high_scores) > 0:
                unique_transcripts = high_scores['transcript_id'].nunique()
                print(f"  高スコア転写産物数: {unique_transcripts}")
        
        # フィルタリング済みデータ保存
        df_filtered.to_csv(output_file, sep='\t', index=False)
        print(f"\\nフィルタリング済みRBRPスコア保存: {output_file}")
        print(f"最終データポイント数: {len(df_filtered)}")
        
    except Exception as e:
        print(f"エラー: RBRPスコアフィルタリング中にエラーが発生しました - {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="RBRPスコアをフィルタリング")
    parser.add_argument("-i", "--input", required=True, help="入力RBRPスコアファイル")
    parser.add_argument("-o", "--output", required=True, help="出力フィルタリング済みファイル")
    parser.add_argument("-t", "--min_depth", type=int, default=200, help="最小シーケンス深度")
    parser.add_argument("-s", "--skip_start", type=int, default=5, help="開始部分スキップ (nt)")
    parser.add_argument("-e", "--skip_end", type=int, default=30, help="終端部分スキップ (nt)")
    
    args = parser.parse_args()
    
    filter_rbrp_scores(args.input, args.output, args.min_depth, args.skip_start, args.skip_end)

if __name__ == "__main__":
    main()