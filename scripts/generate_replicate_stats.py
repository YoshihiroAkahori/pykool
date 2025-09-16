#!/usr/bin/env python3
"""
Replicate統計レポート生成スクリプト
統合されたreplicateの統計情報を生成
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys
import json
from datetime import datetime

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

def generate_replicate_statistics(condition_name, replicate_count, method, output_file):
    """
    Replicate統計レポートを生成

    Args:
        condition_name (str): 実験条件名
        replicate_count (int): replicate数
        method (str): 統計手法
        output_file (str): 出力ファイルパス
    """
    logger = logging.getLogger(__name__)

    stats_report = {
        'condition': condition_name,
        'replicate_count': replicate_count,
        'statistical_method': method,
        'generated_at': datetime.now().isoformat(),
        'summary': {},
        'quality_metrics': {},
        'recommendations': []
    }

    # 基本統計情報
    stats_report['summary'] = {
        'total_replicates': replicate_count,
        'merge_method': method,
        'statistical_power': 'adequate' if replicate_count >= 3 else 'minimal' if replicate_count == 2 else 'insufficient'
    }

    # 品質メトリクス
    if replicate_count >= 3:
        stats_report['quality_metrics']['outlier_detection'] = 'enabled'
        stats_report['quality_metrics']['statistical_significance'] = 'testable'
    else:
        stats_report['quality_metrics']['outlier_detection'] = 'limited'
        stats_report['quality_metrics']['statistical_significance'] = 'limited'

    # 推奨事項
    if replicate_count < 3:
        stats_report['recommendations'].append(
            "統計的検出力向上のため、replicate数をN=3以上に増やすことを推奨"
        )

    if method == 'mean':
        stats_report['recommendations'].append(
            "平均値ベースの統合。外れ値に敏感なため、事前の品質確認を推奨"
        )
    elif method == 'median':
        stats_report['recommendations'].append(
            "中央値ベースの統合。外れ値に対してロバスト"
        )

    # 出力ディレクトリ作成
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # テキストレポート生成
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("# Replicate統計レポート\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"実験条件: {condition_name}\n")
        f.write(f"生成日時: {stats_report['generated_at']}\n")
        f.write(f"Replicate数: {replicate_count}\n")
        f.write(f"統計手法: {method}\n\n")

        f.write("## 基本統計情報\n")
        f.write("-" * 30 + "\n")
        for key, value in stats_report['summary'].items():
            f.write(f"{key}: {value}\n")

        f.write("\n## 品質メトリクス\n")
        f.write("-" * 30 + "\n")
        for key, value in stats_report['quality_metrics'].items():
            f.write(f"{key}: {value}\n")

        f.write("\n## 推奨事項\n")
        f.write("-" * 30 + "\n")
        for i, recommendation in enumerate(stats_report['recommendations'], 1):
            f.write(f"{i}. {recommendation}\n")

        f.write("\n## 統計解釈ガイド\n")
        f.write("-" * 30 + "\n")
        f.write("- N=2: 基本的な比較可能、統計検定は制限的\n")
        f.write("- N=3: 外れ値検出可能、基本的な統計検定可能\n")
        f.write("- N≥4: 信頼性の高い統計解析可能\n\n")

        if method == 'mean':
            f.write("平均値統合:\n")
            f.write("- 利点: 直感的、計算が単純\n")
            f.write("- 注意: 外れ値の影響を受けやすい\n")
        elif method == 'median':
            f.write("中央値統合:\n")
            f.write("- 利点: 外れ値に対してロバスト\n")
            f.write("- 注意: サンプルサイズが小さい場合は精度が落ちる\n")
        elif method == 'sum':
            f.write("合計値統合:\n")
            f.write("- 利点: カウントデータに適している\n")
            f.write("- 注意: スケールの影響を受けやすい\n")

    # JSON形式でも保存
    json_file = output_path.with_suffix('.json')
    with open(json_file, 'w', encoding='utf-8') as f:
        json.dump(stats_report, f, indent=2, ensure_ascii=False)

    logger.info(f"統計レポート生成完了: {output_path}")
    logger.info(f"JSON形式: {json_file}")

def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(description="Replicate統計レポート生成")
    parser.add_argument('-c', '--condition', required=True,
                       help='実験条件名')
    parser.add_argument('-n', '--replicates', type=int, required=True,
                       help='replicate数')
    parser.add_argument('-m', '--method', required=True,
                       choices=['mean', 'median', 'sum'],
                       help='統計手法')
    parser.add_argument('-o', '--output', required=True,
                       help='出力ファイルパス')

    args = parser.parse_args()

    # ログ設定
    logger = setup_logging()

    try:
        logger.info(f"統計レポート生成開始: {args.condition}")

        generate_replicate_statistics(
            args.condition,
            args.replicates,
            args.method,
            args.output
        )

        logger.info("✅ 統計レポート生成完了")

    except Exception as e:
        logger.error(f"❌ エラー: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()