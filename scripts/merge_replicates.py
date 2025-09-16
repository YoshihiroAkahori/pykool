#!/usr/bin/env python3
"""
Replicate統合スクリプト
同一実験条件のreplicateファイルを統計的に統合
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys
from scipy import stats

def setup_logging():
    """ログ設定"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def detect_file_format(file_path):
    """
    ファイル形式を自動検出

    Args:
        file_path (str): ファイルパス

    Returns:
        str: ファイル形式 ('rt', 'rpkm', 'rbrp')
    """
    file_path = Path(file_path)

    if '.rt' in file_path.name:
        return 'rt'
    elif '.rpkm' in file_path.name:
        return 'rpkm'
    elif '.rbrp' in file_path.name or '.out' in file_path.name:
        return 'rbrp'
    else:
        # ファイル内容から推定
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                if 'gene_id' in first_line or 'transcript_id' in first_line:
                    return 'rpkm'
                elif 'position' in first_line or 'coord' in first_line:
                    return 'rt'
                else:
                    return 'rbrp'
        except:
            return 'unknown'

def load_data_file(file_path, file_format):
    """
    データファイルを読み込み

    Args:
        file_path (str): ファイルパス
        file_format (str): ファイル形式

    Returns:
        pd.DataFrame: データフレーム
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"ファイルが見つかりません: {file_path}")

    try:
        if file_format == 'rt':
            # RTファイル形式: transcript_id, position, rt_count
            df = pd.read_csv(file_path, sep='\t', header=0)
            if 'transcript_id' not in df.columns:
                df.columns = ['transcript_id', 'position', 'rt_count']

        elif file_format == 'rpkm':
            # RPKMファイル形式: gene_id, rpkm_value
            df = pd.read_csv(file_path, sep='\t', header=0)
            if 'gene_id' not in df.columns:
                df.columns = ['gene_id', 'rpkm_value']

        elif file_format == 'rbrp':
            # RBRPファイル形式: transcript_id, position, rbrp_score
            df = pd.read_csv(file_path, sep='\t', header=0)
            if len(df.columns) >= 3:
                df.columns = df.columns[:3].tolist() + ['rbrp_score'] if len(df.columns) > 3 else df.columns.tolist()

        else:
            # 汎用的な読み込み
            df = pd.read_csv(file_path, sep='\t', header=0)

        return df

    except Exception as e:
        raise ValueError(f"ファイル読み込みエラー {file_path}: {e}")

def merge_replicates(input_files, method='mean', outlier_detection=True):
    """
    Replicateファイルを統合

    Args:
        input_files (list): 入力ファイルパスのリスト
        method (str): 統計手法 ('mean', 'median', 'sum')
        outlier_detection (bool): 外れ値検出を実行

    Returns:
        pd.DataFrame: 統合されたデータフレーム
        dict: 統計情報
    """
    logger = logging.getLogger(__name__)

    if len(input_files) < 2:
        raise ValueError("最低2つのreplicateファイルが必要です")

    # ファイル形式を検出
    file_format = detect_file_format(input_files[0])
    logger.info(f"ファイル形式検出: {file_format}")

    # データ読み込み
    dfs = []
    for i, file_path in enumerate(input_files):
        logger.info(f"読み込み中: {file_path}")
        df = load_data_file(file_path, file_format)
        df = df.add_suffix(f'_rep{i+1}' if i > 0 else '')
        dfs.append(df)

    # データフレーム結合
    if file_format == 'rt':
        # transcript_id と position でマージ
        merged = dfs[0]
        for i, df in enumerate(dfs[1:], 1):
            key_cols = ['transcript_id', 'position']
            merged = merged.merge(df, on=key_cols, how='outer', suffixes=('', f'_rep{i+1}'))

        value_cols = [col for col in merged.columns if 'rt_count' in col]

    elif file_format == 'rpkm':
        # gene_id でマージ
        merged = dfs[0]
        for i, df in enumerate(dfs[1:], 1):
            merged = merged.merge(df, on='gene_id', how='outer', suffixes=('', f'_rep{i+1}'))

        value_cols = [col for col in merged.columns if 'rpkm' in col]

    elif file_format == 'rbrp':
        # transcript_id と position でマージ
        merged = dfs[0]
        for i, df in enumerate(dfs[1:], 1):
            key_cols = [merged.columns[0], merged.columns[1]]  # 最初の2列をキーとして使用
            df_renamed = df.copy()
            df_renamed.columns = merged.columns[:2].tolist() + [f'{merged.columns[2]}_rep{i+1}'] + df.columns[3:].tolist()
            merged = merged.merge(df_renamed, on=key_cols, how='outer')

        value_cols = [col for col in merged.columns if 'rbrp' in col or col == merged.columns[2]]

    else:
        raise ValueError(f"不明なファイル形式: {file_format}")

    # 欠損値の処理
    merged[value_cols] = merged[value_cols].fillna(0)

    # 外れ値検出
    outliers_info = {}
    if outlier_detection and len(value_cols) >= 3:
        logger.info("外れ値検出実行中...")

        # Z-scoreによる外れ値検出
        z_scores = np.abs(stats.zscore(merged[value_cols], axis=1, nan_policy='omit'))
        outlier_threshold = 2.5

        outlier_mask = z_scores > outlier_threshold
        outliers_info = {
            'outlier_count': np.sum(outlier_mask),
            'outlier_threshold': outlier_threshold,
            'outlier_positions': merged[outlier_mask.any(axis=1)].index.tolist()
        }

        logger.info(f"外れ値検出結果: {outliers_info['outlier_count']} 箇所")

    # 統計的統合
    if method == 'mean':
        merged['merged_value'] = merged[value_cols].mean(axis=1)
        merged['std_dev'] = merged[value_cols].std(axis=1)
    elif method == 'median':
        merged['merged_value'] = merged[value_cols].median(axis=1)
        merged['mad'] = merged[value_cols].mad(axis=1)  # Median Absolute Deviation
    elif method == 'sum':
        merged['merged_value'] = merged[value_cols].sum(axis=1)
        merged['std_dev'] = merged[value_cols].std(axis=1)
    else:
        raise ValueError(f"不明な統計手法: {method}")

    # 統計情報
    stats_info = {
        'method': method,
        'input_files': len(input_files),
        'total_features': len(merged),
        'missing_values': merged[value_cols].isnull().sum().sum(),
        'correlation_matrix': merged[value_cols].corr().to_dict() if len(value_cols) > 1 else {},
        'outliers': outliers_info
    }

    logger.info(f"統合完了: {method}手法, {len(merged)} features")

    return merged, stats_info

def save_merged_data(merged_df, output_file, file_format, stats_info):
    """
    統合データを保存

    Args:
        merged_df (pd.DataFrame): 統合データフレーム
        output_file (str): 出力ファイルパス
        file_format (str): ファイル形式
        stats_info (dict): 統計情報
    """
    logger = logging.getLogger(__name__)

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # メインデータ保存
    merged_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"統合データ保存: {output_path}")

    # 統計情報保存
    stats_file = output_path.with_suffix('.stats.txt')
    with open(stats_file, 'w') as f:
        f.write("# Replicate Merge Statistics\n")
        f.write(f"Method: {stats_info['method']}\n")
        f.write(f"Input files: {stats_info['input_files']}\n")
        f.write(f"Total features: {stats_info['total_features']}\n")
        f.write(f"Missing values: {stats_info['missing_values']}\n")

        if stats_info['outliers']:
            f.write(f"Outliers detected: {stats_info['outliers']['outlier_count']}\n")
            f.write(f"Outlier threshold: {stats_info['outliers']['outlier_threshold']}\n")

        if stats_info['correlation_matrix']:
            f.write("\n# Correlation Matrix\n")
            for col1, correlations in stats_info['correlation_matrix'].items():
                for col2, corr in correlations.items():
                    if col1 != col2:
                        f.write(f"{col1} vs {col2}: {corr:.3f}\n")

    logger.info(f"統計情報保存: {stats_file}")

def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(description="Replicate統合スクリプト")
    parser.add_argument('-i', '--input', required=True,
                       help='入力ファイル（:区切りで複数指定）')
    parser.add_argument('-o', '--output', required=True,
                       help='出力ファイルパス')
    parser.add_argument('-m', '--method', default='mean',
                       choices=['mean', 'median', 'sum'],
                       help='統計手法 (default: mean)')
    parser.add_argument('-c', '--condition',
                       help='実験条件名（ログ用）')
    parser.add_argument('--no-outlier-detection', action='store_true',
                       help='外れ値検出を無効化')

    args = parser.parse_args()

    # ログ設定
    logger = setup_logging()

    try:
        # 入力ファイル解析
        input_files = args.input.split(':')
        logger.info(f"統合開始: {len(input_files)} replicates")
        logger.info(f"条件: {args.condition or 'Unknown'}")
        logger.info(f"手法: {args.method}")

        # ファイル存在確認
        for file_path in input_files:
            if not Path(file_path).exists():
                raise FileNotFoundError(f"入力ファイルが見つかりません: {file_path}")

        # Replicate統合実行
        merged_df, stats_info = merge_replicates(
            input_files,
            method=args.method,
            outlier_detection=not args.no_outlier_detection
        )

        # 結果保存
        file_format = detect_file_format(input_files[0])
        save_merged_data(merged_df, args.output, file_format, stats_info)

        logger.info("✅ Replicate統合完了")

    except Exception as e:
        logger.error(f"❌ エラー: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()