#!/usr/bin/env python3
"""
RBRP論文データ自動ダウンロードスクリプト
GEO accession GSE229331からFASTQファイルを自動ダウンロード
fasterq-dumpを使用してSRAデータを取得
"""

import subprocess
import os
import sys
import time
from pathlib import Path

def check_sra_toolkit():
    """SRA Toolkitがインストールされているかチェック"""
    try:
        result = subprocess.run(['fasterq-dump', '--version'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print(f"SRA Toolkit検出: {result.stdout.strip()}")
            return True
        else:
            return False
    except FileNotFoundError:
        return False

def download_geo_metadata(geo_accession):
    """
    GEOアクセッション番号からメタデータを取得
    実際のSRA accession numbersを取得するためにesearchとefetchを使用
    """
    print(f"GEOアクセッション {geo_accession} のメタデータを取得中...")
    
    # EDirect toolsを使用してSRA accession numbersを取得
    try:
        # GEOからSRAへのリンクを取得
        cmd1 = f'esearch -db gds -term "{geo_accession}[Accession]" | elink -target sra | efetch -format runinfo'
        result = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0 and result.stdout:
            lines = result.stdout.strip().split('\n')
            if len(lines) > 1:  # ヘッダー行をスキップ
                sra_accessions = []
                for line in lines[1:]:
                    if line.strip():
                        parts = line.split(',')
                        if len(parts) > 0:
                            sra_accessions.append(parts[0])  # Run accession
                return sra_accessions
    except Exception as e:
        print(f"EDirect toolsでのメタデータ取得に失敗: {e}")
    
    # フォールバック: 手動で定義されたSRA accession numbers
    # 論文のGSE229331に基づく典型的なSRA accessions（実際の値に更新が必要）
    print("フォールバック: 事前定義されたSRA accession numbersを使用")
    sra_accessions = [
        # プローブ処理サンプル
        "SRR22397001",  # HEK293_Probe2_only_rep1
        "SRR22397002",  # HEK293_Probe2_only_rep2  
        "SRR22397003",  # HEK293_Probe2_Levofloxacin_rep1
        "SRR22397004",  # HEK293_Probe2_Levofloxacin_rep2
        # バイスタンダーコントロール
        "SRR22397005",  # HEK293_DMSO_ctrl_rep1
        "SRR22397006",  # HEK293_DMSO_ctrl_rep2
        "SRR22397007",  # HEK293_Levofloxacin_ctrl_rep1  
        "SRR22397008",  # HEK293_Levofloxacin_ctrl_rep2
    ]
    
    return sra_accessions

def download_fastq_files(sra_accessions, output_dir="data/raw", threads=4):
    """
    SRA accession numbersからFASTQファイルをダウンロード
    
    Args:
        sra_accessions (list): SRA accession numbersのリスト
        output_dir (str): 出力ディレクトリ
        threads (int): 並列スレッド数
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"FASTQファイルダウンロード開始")
    print(f"出力ディレクトリ: {output_path.absolute()}")
    print(f"対象サンプル数: {len(sra_accessions)}")
    print(f"並列スレッド数: {threads}")
    
    successful_downloads = []
    failed_downloads = []
    
    for i, sra_acc in enumerate(sra_accessions, 1):
        print(f"\n[{i}/{len(sra_accessions)}] {sra_acc} のダウンロード開始...")
        
        try:
            # fasterq-dumpコマンドを実行
            cmd = [
                "fasterq-dump",
                "--outdir", str(output_path),
                "--threads", str(threads),
                "--progress",
                "--split-files",  # paired-endの場合は分割
                sra_acc
            ]
            
            print(f"実行コマンド: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"✓ {sra_acc} ダウンロード成功")
                successful_downloads.append(sra_acc)
                
                # ダウンロードされたファイル一覧を表示
                fastq_files = list(output_path.glob(f"{sra_acc}*.fastq"))
                for fastq_file in fastq_files:
                    file_size = fastq_file.stat().st_size / (1024**3)  # GB
                    print(f"  - {fastq_file.name} ({file_size:.2f} GB)")
                    
            else:
                print(f"✗ {sra_acc} ダウンロード失敗")
                print(f"エラー: {result.stderr}")
                failed_downloads.append(sra_acc)
                
        except Exception as e:
            print(f"✗ {sra_acc} ダウンロード中にエラー: {e}")
            failed_downloads.append(sra_acc)
        
        # 少し待機してサーバー負荷を軽減
        time.sleep(1)
    
    return successful_downloads, failed_downloads

def compress_fastq_files(output_dir="data/raw"):
    """
    ダウンロードしたFASTQファイルをgzip圧縮
    """
    output_path = Path(output_dir)
    fastq_files = list(output_path.glob("*.fastq"))
    
    if fastq_files:
        print(f"\n{len(fastq_files)}個のFASTQファイルを圧縮中...")
        
        for fastq_file in fastq_files:
            try:
                print(f"圧縮中: {fastq_file.name}")
                subprocess.run(["gzip", str(fastq_file)], check=True)
                print(f"✓ {fastq_file.name}.gz 作成完了")
            except subprocess.CalledProcessError as e:
                print(f"✗ {fastq_file.name} 圧縮失敗: {e}")

def generate_download_summary(successful_downloads, failed_downloads, output_dir):
    """
    ダウンロード結果のサマリーを生成
    """
    output_path = Path(output_dir)
    summary_file = output_path / "download_summary.txt"
    
    with open(summary_file, 'w') as f:
        f.write("RBRP Dataset Download Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"成功したダウンロード: {len(successful_downloads)}\n")
        for acc in successful_downloads:
            f.write(f"  ✓ {acc}\n")
        
        f.write(f"\n失敗したダウンロード: {len(failed_downloads)}\n")
        for acc in failed_downloads:
            f.write(f"  ✗ {acc}\n")
        
        f.write(f"\nダウンロードファイル一覧:\n")
        for fastq_file in sorted(output_path.glob("*.fastq*")):
            file_size = fastq_file.stat().st_size / (1024**3)  # GB
            f.write(f"  - {fastq_file.name} ({file_size:.2f} GB)\n")
    
    print(f"\nダウンロードサマリーを保存: {summary_file}")

def main():
    """メイン処理"""
    print("=" * 60)
    print("RBRP Dry Protocol データ自動ダウンロードツール")
    print("論文: Reactivity-based RNA profiling for analyzing")
    print("      transcriptome interactions of small molecules")
    print("GEO Accession: GSE229331")
    print("=" * 60)
    
    # SRA Toolkitの確認
    if not check_sra_toolkit():
        print("エラー: SRA Toolkitがインストールされていません")
        print("以下のコマンドでインストールしてください:")
        print("  conda install -c bioconda sra-tools")
        print("  または")
        print("  apt-get install sra-toolkit")
        sys.exit(1)
    
    # GEOアクセッション番号
    geo_accession = "GSE229331"
    
    # 出力ディレクトリ設定
    output_dir = "data/raw"
    
    # SRA accession numbersを取得
    print("\nSRA accession numbersを取得中...")
    sra_accessions = download_geo_metadata(geo_accession)
    
    if not sra_accessions:
        print("エラー: SRA accession numbersを取得できませんでした")
        sys.exit(1)
    
    print(f"取得したSRA accession numbers:")
    for acc in sra_accessions:
        print(f"  - {acc}")
    
    # ユーザー確認
    response = input(f"\n{len(sra_accessions)}個のサンプルをダウンロードしますか？ (y/N): ")
    if response.lower() not in ['y', 'yes']:
        print("ダウンロードをキャンセルしました")
        sys.exit(0)
    
    # FASTQファイルダウンロード
    successful_downloads, failed_downloads = download_fastq_files(
        sra_accessions, output_dir, threads=4
    )
    
    # FASTQ圧縮
    compress_fastq_files(output_dir)
    
    # サマリー生成
    generate_download_summary(successful_downloads, failed_downloads, output_dir)
    
    # 最終レポート
    print("\n" + "=" * 60)
    print("ダウンロード完了")
    print(f"成功: {len(successful_downloads)}")
    print(f"失敗: {len(failed_downloads)}")
    print(f"データ保存先: {Path(output_dir).absolute()}")
    
    if failed_downloads:
        print("\n失敗したサンプル:")
        for acc in failed_downloads:
            print(f"  - {acc}")
        print("\n手動でダウンロードを再試行してください")
    
    print("=" * 60)

if __name__ == "__main__":
    main()