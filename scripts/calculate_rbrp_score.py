#!/usr/bin/env python3
"""
RBRPスコア計算スクリプト
プローブ処理サンプルとコントロールサンプル間でRBRPスコアを計算
"""

import argparse
import pandas as pd
import numpy as np
import sys
from collections import defaultdict

def load_rt_file(rt_file):
    """
    RTファイルを読み込み
    
    Args:
        rt_file (str): RTファイルパス
    
    Returns:
        dict: 転写産物ごとのRTstopデータ
    """
    rt_data = defaultdict(dict)
    
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
                    
                    rt_data[transcript_id][position] = {
                        'rtstop_count': rtstop_count,
                        'base_coverage': base_coverage,
                        'rtstop_ratio': rtstop_ratio
                    }
        
        print(f"RTファイル読み込み完了: {rt_file}")
        print(f"転写産物数: {len(rt_data)}")
        return rt_data
        
    except Exception as e:
        print(f"エラー: RTファイル読み込み失敗 - {e}")
        return {}

def calculate_rbrp_scores(probe_rt_data, background_rt_data=None, method='dividing', factor=0.5):
    """
    RBRPスコア計算
    
    Args:
        probe_rt_data (dict): プローブ処理サンプルのRTデータ
        background_rt_data (dict): バックグラウンドサンプルのRTデータ
        method (str): 計算方法
        factor (float): バックグラウンド補正係数
    
    Returns:
        dict: RBRPスコア
    """
    rbrp_scores = defaultdict(dict)
    
    for transcript_id in probe_rt_data:
        for position in probe_rt_data[transcript_id]:
            probe_ratio = probe_rt_data[transcript_id][position]['rtstop_ratio']
            
            # バックグラウンド補正
            if background_rt_data and transcript_id in background_rt_data and position in background_rt_data[transcript_id]:
                bg_ratio = background_rt_data[transcript_id][position]['rtstop_ratio']
            else:
                bg_ratio = 0.0
            
            # RBRPスコア計算
            if method == 'dividing':
                rbrp_score = probe_ratio - (factor * bg_ratio)
            else:
                rbrp_score = probe_ratio
            
            rbrp_scores[transcript_id][position] = {
                'probe_ratio': probe_ratio,
                'background_ratio': bg_ratio,
                'rbrp_score': rbrp_score,
                'base_coverage': probe_rt_data[transcript_id][position]['base_coverage']
            }
    
    return rbrp_scores

def save_rbrp_scores(rbrp_scores, output_file):
    """
    RBRPスコアをファイルに保存
    
    Args:
        rbrp_scores (dict): RBRPスコアデータ
        output_file (str): 出力ファイルパス
    """
    try:
        with open(output_file, 'w') as f:
            f.write("transcript_id\\tposition\\tprobe_ratio\\tbackground_ratio\\trbrp_score\\tbase_coverage\\n")
            
            for transcript_id in rbrp_scores:
                for position in rbrp_scores[transcript_id]:
                    data = rbrp_scores[transcript_id][position]
                    f.write(f"{transcript_id}\\t{position}\\t{data['probe_ratio']:.6f}\\t"
                           f"{data['background_ratio']:.6f}\\t{data['rbrp_score']:.6f}\\t"
                           f"{data['base_coverage']}\\n")
        
        print(f"RBRPスコア保存: {output_file}")
        
        # 統計情報
        all_scores = []
        for transcript_data in rbrp_scores.values():
            for pos_data in transcript_data.values():
                all_scores.append(pos_data['rbrp_score'])
        
        if all_scores:
            print(f"RBRPスコア統計:")
            print(f"  平均: {np.mean(all_scores):.4f}")
            print(f"  中央値: {np.median(all_scores):.4f}")
            print(f"  標準偏差: {np.std(all_scores):.4f}")
            print(f"  最大値: {np.max(all_scores):.4f}")
            print(f"  最小値: {np.min(all_scores):.4f}")
        
    except Exception as e:
        print(f"エラー: RBRPスコア保存失敗 - {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="RBRPスコアを計算")
    parser.add_argument("-f", "--foreground", required=True, help="プローブ処理サンプルのRTファイル")
    parser.add_argument("-b", "--background", help="バックグラウンドサンプルのRTファイル")
    parser.add_argument("-o", "--output", required=True, help="出力RBRPスコアファイル")
    parser.add_argument("-e", "--method", default="dividing", help="計算方法")
    parser.add_argument("-y", "--factor", type=float, default=0.5, help="バックグラウンド補正係数")
    
    args = parser.parse_args()
    
    # プローブデータ読み込み
    probe_rt_data = load_rt_file(args.foreground)
    if not probe_rt_data:
        print("エラー: プローブデータの読み込みに失敗しました")
        sys.exit(1)
    
    # バックグラウンドデータ読み込み（オプション）
    background_rt_data = None
    if args.background:
        background_rt_data = load_rt_file(args.background)
    
    # RBRPスコア計算
    print("RBRPスコア計算中...")
    rbrp_scores = calculate_rbrp_scores(probe_rt_data, background_rt_data, args.method, args.factor)
    
    # 結果保存
    save_rbrp_scores(rbrp_scores, args.output)
    print("RBRPスコア計算完了")

if __name__ == "__main__":
    main()