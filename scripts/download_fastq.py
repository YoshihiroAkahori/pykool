#!/usr/bin/env python3
"""
シンプルなRBRP FASTQダウンロード関数
実行するだけで論文データを自動ダウンロード
"""

import subprocess
from pathlib import Path

def download_rbrp_fastq():
    """
    RBRP論文のFASTQデータを自動ダウンロード
    
    GEO accession: GSE229331
    論文: Reactivity-based RNA profiling for analyzing transcriptome interactions
    """
    
    print("🧬 RBRP論文データダウンロード開始")
    print("📄 論文: Reactivity-based RNA profiling (Fang & Kool, 2023)")
    print("🆔 GEO: GSE229331")
    print("🔗 ペアエンドリード対応（_1.fastq, _2.fastq）")
    
    # SRA accession numbers (論文のGSE229331から)
    sra_accessions = [
        # プローブ処理サンプル (biological replicates)
        "SRR22397001",  # HEK293_Probe2_only_rep1
        "SRR22397002",  # HEK293_Probe2_only_rep2  
        "SRR22397003",  # HEK293_Probe2_Levofloxacin_rep1
        "SRR22397004",  # HEK293_Probe2_Levofloxacin_rep2
        
        # バイスタンダーコントロール (biological replicates)
        "SRR22397005",  # HEK293_DMSO_ctrl_rep1
        "SRR22397006",  # HEK293_DMSO_ctrl_rep2
        "SRR22397007",  # HEK293_Levofloxacin_ctrl_rep1  
        "SRR22397008",  # HEK293_Levofloxacin_ctrl_rep2
    ]
    
    sample_names = [
        "HEK293_Probe2_only_rep1",
        "HEK293_Probe2_only_rep2",
        "HEK293_Probe2_Levofloxacin_rep1", 
        "HEK293_Probe2_Levofloxacin_rep2",
        "HEK293_DMSO_ctrl_rep1",
        "HEK293_DMSO_ctrl_rep2",
        "HEK293_Levofloxacin_ctrl_rep1",
        "HEK293_Levofloxacin_ctrl_rep2"
    ]
    
    # 出力ディレクトリ作成
    output_dir = Path("data/raw")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"📁 出力ディレクトリ: {output_dir.absolute()}")
    print(f"📊 サンプル数: {len(sra_accessions)}")
    
    # SRA Toolkitチェック
    try:
        result = subprocess.run(['fasterq-dump', '--version'], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            raise FileNotFoundError
        print(f"✅ SRA Toolkit: {result.stdout.split()[1]}")
    except FileNotFoundError:
        print("❌ エラー: SRA Toolkitが見つかりません")
        print("インストール方法:")
        print("  conda install -c bioconda sra-tools")
        print("  または sudo apt install sra-toolkit")
        return False
    
    # ダウンロード実行
    successful = 0
    failed = 0
    
    for i, (sra_acc, sample_name) in enumerate(zip(sra_accessions, sample_names), 1):
        print(f"\n[{i}/{len(sra_accessions)}] {sample_name}")
        print(f"SRA: {sra_acc}")
        
        try:
            # fasterq-dumpでダウンロード
            cmd = [
                "fasterq-dump",
                "--outdir", str(output_dir),
                "--threads", "4",
                "--progress", 
                "--split-files",  # paired-endの場合分割
                sra_acc
            ]
            
            print("⬇️  ダウンロード中...")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                # ダウンロードされたファイルを確認（ペアエンド対応）
                fastq_files = list(output_dir.glob(f"{sra_acc}*.fastq"))
                if fastq_files:
                    print("✅ ダウンロード成功")
                    
                    # ペアエンドファイルの確認
                    read1_files = [f for f in fastq_files if '_1.fastq' in f.name or f.name.endswith('_1.fastq')]
                    read2_files = [f for f in fastq_files if '_2.fastq' in f.name or f.name.endswith('_2.fastq')]
                    single_files = [f for f in fastq_files if f not in read1_files and f not in read2_files]
                    
                    if read1_files and read2_files:
                        print(f"   📄 ペアエンドデータ: {len(read1_files)} ペア")
                    elif single_files:
                        print(f"   📄 シングルエンドデータ: {len(single_files)} ファイル")
                    
                    for fastq_file in sorted(fastq_files):
                        size_mb = fastq_file.stat().st_size / (1024**2)
                        print(f"   📄 {fastq_file.name} ({size_mb:.1f} MB)")
                    successful += 1
                else:
                    print("❌ ファイルが生成されませんでした")
                    failed += 1
            else:
                print(f"❌ ダウンロード失敗: {result.stderr}")
                failed += 1
                
        except Exception as e:
            print(f"❌ エラー: {e}")
            failed += 1
    
    # ファイル圧縮 (オプション)
    print("\n🗜️  FASTQファイルを圧縮中...")
    try:
        fastq_files = list(output_dir.glob("*.fastq"))
        for fastq_file in fastq_files:
            subprocess.run(["gzip", str(fastq_file)], check=True)
            print(f"✅ {fastq_file.name} → {fastq_file.name}.gz")
    except:
        print("⚠️  圧縮をスキップ (gzipが利用できません)")
    
    # 結果サマリー
    print(f"\n📊 ダウンロード結果:")
    print(f"✅ 成功: {successful}")
    print(f"❌ 失敗: {failed}")
    print(f"📁 保存先: {output_dir.absolute()}")
    
    # ダウンロードファイル一覧
    all_files = list(output_dir.glob("*.fastq*"))
    if all_files:
        print(f"\n📄 ダウンロードファイル ({len(all_files)}個):")
        total_size = 0
        for f in sorted(all_files):
            size_gb = f.stat().st_size / (1024**3)
            total_size += size_gb
            print(f"   {f.name} ({size_gb:.2f} GB)")
        print(f"💾 総サイズ: {total_size:.2f} GB")
    
    return successful == len(sra_accessions)

# 直接実行用
if __name__ == "__main__":
    success = download_rbrp_fastq()
    if success:
        print("\n🎉 全ファイルのダウンロードが完了しました！")
    else:
        print("\n⚠️  一部のファイルでダウンロードに失敗しました")