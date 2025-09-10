#!/usr/bin/env python3
"""
Bedgraphファイル生成スクリプト
RBRPスコアからIGV用のbedgraphファイルを生成
"""

import argparse
import pandas as pd
import sys
from collections import defaultdict

def load_gtf_info(gtf_file):
    """
    GTFファイルから転写産物情報を読み込み
    
    Args:
        gtf_file (str): GTFファイルパス
    
    Returns:
        dict: 転写産物の染色体位置情報
    """
    transcript_info = {}
    
    if not gtf_file:
        print("警告: GTFファイルが指定されていません")
        return transcript_info
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == 'transcript':
                    chrom = parts[0]
                    start = int(parts[3]) - 1  # 0-based
                    end = int(parts[4])
                    strand = parts[6]
                    
                    # transcript_idを抽出
                    attributes = parts[8]
                    transcript_id = None
                    for attr in attributes.split(';'):
                        if 'transcript_id' in attr:
                            transcript_id = attr.split('"')[1]
                            break
                    
                    if transcript_id:
                        transcript_info[transcript_id] = {
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'strand': strand
                        }
        
        print(f"GTF情報読み込み完了: {len(transcript_info)} 転写産物")
        return transcript_info
        
    except Exception as e:
        print(f"警告: GTFファイル読み込みエラー - {e}")
        return {}

def generate_bedgraph(rbrp_file, output_file, gtf_file=None):
    """
    RBRPスコアからbedgraphファイルを生成
    
    Args:
        rbrp_file (str): 入力RBRPスコアファイル
        output_file (str): 出力bedgraphファイル
        gtf_file (str): GTFアノテーションファイル（オプション）
    """
    print(f"Bedgraph生成開始")
    print(f"入力RBRPファイル: {rbrp_file}")
    print(f"出力bedgraphファイル: {output_file}")
    
    # GTF情報読み込み
    transcript_info = load_gtf_info(gtf_file) if gtf_file else {}
    
    try:
        # RBRPスコア読み込み
        df = pd.read_csv(rbrp_file, sep='\t')
        print(f"RBRPデータポイント数: {len(df)}")
        
        # Bedgraphエントリ生成
        bedgraph_entries = []
        
        for _, row in df.iterrows():
            transcript_id = row['transcript_id']
            position = int(row['position'])
            rbrp_score = float(row['rbrp_score'])
            
            # 染色体位置情報があるかチェック
            if transcript_id in transcript_info:
                chrom = transcript_info[transcript_id]['chrom']
                transcript_start = transcript_info[transcript_id]['start']
                
                # ゲノム座標に変換
                genomic_pos = transcript_start + position
                
                bedgraph_entries.append({
                    'chrom': chrom,
                    'start': genomic_pos,
                    'end': genomic_pos + 1,
                    'score': rbrp_score
                })
            else:
                # GTF情報がない場合は転写産物IDをそのまま使用
                bedgraph_entries.append({
                    'chrom': transcript_id,
                    'start': position,
                    'end': position + 1,
                    'score': rbrp_score
                })
        
        # Bedgraphファイル書き込み
        if bedgraph_entries:
            bedgraph_df = pd.DataFrame(bedgraph_entries)
            
            # 染色体とポジションでソート
            bedgraph_df = bedgraph_df.sort_values(['chrom', 'start'])
            
            with open(output_file, 'w') as f:
                # Bedgraphヘッダー
                f.write("track type=bedGraph name=\\\"RBRP_Score\\\" description=\\\"RBRP Score Track\\\"\\n")
                
                # データエントリ
                for _, row in bedgraph_df.iterrows():
                    f.write(f"{row['chrom']}\\t{row['start']}\\t{row['end']}\\t{row['score']:.6f}\\n")
            
            print(f"Bedgraph生成完了: {output_file}")
            print(f"エントリ数: {len(bedgraph_entries)}")
            
            # 統計情報
            print(f"\\nBedgraphスコア統計:")
            print(f"  平均: {bedgraph_df['score'].mean():.4f}")
            print(f"  中央値: {bedgraph_df['score'].median():.4f}")
            print(f"  最大値: {bedgraph_df['score'].max():.4f}")
            print(f"  最小値: {bedgraph_df['score'].min():.4f}")
            
        else:
            print("警告: Bedgraphエントリが生成されませんでした")
            
    except Exception as e:
        print(f"エラー: Bedgraph生成中にエラーが発生しました - {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="RBRPスコアからbedgraphファイルを生成")
    parser.add_argument("-i", "--input", required=True, help="入力RBRPスコアファイル")
    parser.add_argument("-o", "--output", required=True, help="出力bedgraphファイル")
    parser.add_argument("-g", "--gtf", help="GTFアノテーションファイル")
    parser.add_argument("-a", "--annotation", help="追加アノテーションファイル（使用しない）")
    
    args = parser.parse_args()
    
    generate_bedgraph(args.input, args.output, args.gtf)

if __name__ == "__main__":
    main()