#!/usr/bin/env python3
"""
RTファイル正規化スクリプト
RTファイルのフォーマットを正規化
"""

import argparse
import pandas as pd
import sys

def normalize_rt_file(input_file, output_file, skip_start=32, skip_end=32):
    """
    RTファイルを正規化
    
    Args:
        input_file (str): 入力RTファイル
        output_file (str): 出力正規化ファイル
        skip_start (int): 開始部分スキップ
        skip_end (int): 終端部分スキップ
    """
    print(f"RTファイル正規化開始")
    print(f"入力ファイル: {input_file}")
    print(f"開始スキップ: {skip_start} nt")
    print(f"終端スキップ: {skip_end} nt")
    
    try:
        # RTファイル読み込み
        df = pd.read_csv(input_file, sep='\t')
        original_count = len(df)
        print(f"元のデータポイント数: {original_count}")
        
        # 転写産物ごとの長さ計算
        transcript_lengths = df.groupby('transcript_id')['position'].max()
        
        # 正規化フィルタリング
        def is_valid_position(row):
            transcript_id = row['transcript_id']
            position = row['position']
            max_pos = transcript_lengths[transcript_id]
            
            # 開始部分と終端部分をスキップ
            return (position >= skip_start) and (position <= max_pos - skip_end)
        
        df_normalized = df[df.apply(is_valid_position, axis=1)].copy()
        normalized_count = len(df_normalized)
        
        print(f"正規化後データポイント数: {normalized_count} ({normalized_count/original_count*100:.1f}%)")
        
        # 正規化済みファイル保存
        df_normalized.to_csv(output_file, sep='\t', index=False)
        print(f"正規化完了: {output_file}")
        
    except Exception as e:
        print(f"エラー: RTファイル正規化中にエラーが発生しました - {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="RTファイルを正規化")
    parser.add_argument("-i", "--input", required=True, help="入力RTファイル")
    parser.add_argument("-o", "--output", required=True, help="出力正規化ファイル")
    parser.add_argument("-d", "--skip_start", type=int, default=32, help="開始部分スキップ")
    parser.add_argument("-l", "--skip_end", type=int, default=32, help="終端部分スキップ")
    
    args = parser.parse_args()
    
    normalize_rt_file(args.input, args.output, args.skip_start, args.skip_end)

if __name__ == "__main__":
    main()