# RBRP Dry Protocol 自動化パイプライン

Eric Koolラボの「Reactivity-based RNA profiling for analyzing transcriptome interactions of small molecules in human cells」論文のドライプロトコル（バイオインフォマティクス解析）を自動化するJupyterノートブックです。

## 概要

このプロジェクトは、RBRP（Reactivity-Based RNA Profiling）プロトコルのシーケンシングデータ処理部分を完全自動化し、バイオインフォマティクス初学者でも簡単に使用できるようにします。

**対象論文**: [Reactivity-based RNA profiling for analyzing transcriptome interactions of small molecules in human cells](https://www.sciencedirect.com/science/article/pii/S2666166723006378)

## 主な機能

- 🔄 **完全自動化**: FASTQファイルを指定するだけで全プロセスが自動実行
- 📊 **品質管理**: FastQCによる品質確認とレポート生成
- 🧬 **シーケンス解析**: デマルチプレックス、トリミング、アライメント
- 📈 **RBRPスコア計算**: RNA-薬物相互作用サイトの同定
- 📋 **可視化**: IGV用bigwigファイル生成
- 📝 **詳細ログ**: 各ステップの実行状況とエラーレポート

## プロジェクト構造

```
kool/
├── notebooks/
│   └── RBRP_Dry_Protocol_Pipeline.ipynb  # メインパイプライン
├── scripts/                              # ヘルパースクリプト
│   ├── calculate_rpkm.py                 # RPKM計算
│   ├── calculate_rtstops.py              # RTstop計算
│   ├── calculate_rbrp_score.py           # RBRPスコア計算
│   ├── filter_rbrp_scores.py             # スコアフィルタリング
│   ├── generate_bedgraph.py              # Bedgraph生成
│   ├── merge_rt_files.py                 # RTファイルマージ
│   └── normalize_rt_file.py              # RTファイル正規化
├── data/
│   ├── raw/                              # 元データ
│   ├── processed/                        # 処理済みデータ
│   └── results/                          # 最終結果
├── logs/                                 # ログファイル
├── requirements.txt                      # Python依存関係
└── README.md                            # このファイル
```

## 必要なソフトウェア

### Python環境
- Python 3.6.6以上
- Jupyter Lab/Notebook

### 外部ツール（ソースからインストール）
以下のバイオインフォマティクスツールが必要です：

#### 1. FastQC
```bash
# 依存関係インストール
sudo apt-get update
sudo apt-get install openjdk-8-jdk

# FastQCダウンロード・インストール
cd /opt
sudo wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
sudo unzip fastqc_v0.12.1.zip
sudo chmod +x FastQC/fastqc
sudo ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc
```

#### 2. Bowtie2
```bash
# 依存関係インストール
sudo apt-get install build-essential libtbb-dev zlib1g-dev

# Bowtie2ソースダウンロード・コンパイル
cd /tmp
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-source.zip
unzip bowtie2-2.5.1-source.zip
cd bowtie2-2.5.1
make
sudo cp bowtie2* /usr/local/bin/
```

#### 3. Cutadapt
```bash
# Pythonパッケージマネージャーでインストール
pip install cutadapt
```

#### 4. UCSC Tools（BigWig生成用）
```bash
# バイナリダウンロード・インストール
cd /tmp
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
chmod +x bedGraphToBigWig wigToBigWig
sudo mv bedGraphToBigWig wigToBigWig /usr/local/bin/
```

#### 5. gffread
```bash
# 依存関係インストール
sudo apt-get install build-essential

# gffreadeソースダウンロード・コンパイル
cd /tmp
git clone https://github.com/gpertea/gffread.git
cd gffread
make
sudo cp gffread /usr/local/bin/
```

#### 6. SAMtools
```bash
# 依存関係インストール
sudo apt-get install build-essential libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev

# SAMtoolsソースダウンロード・コンパイル
cd /tmp
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xjf samtools-1.18.tar.bz2
cd samtools-1.18
./configure --prefix=/usr/local
make
sudo make install
```

#### 7. WiggleTools（オプション）
```bash
# 依存関係インストール
sudo apt-get install libgsl-dev libcurl4-openssl-dev

# WiggleToolsソースダウンロード・コンパイル
cd /tmp
git clone https://github.com/Ensembl/WiggleTools.git
cd WiggleTools
make
sudo cp bin/wiggletools /usr/local/bin/
```

#### 8. icSHAPE/RBRPパイプライン
```bash
# icSHAPEパイプラインダウンロード
cd /opt
sudo git clone https://github.com/qczhang/icSHAPE.git
sudo chmod +x icSHAPE/scripts/*.pl

# RBRPパイプライン追加スクリプト
cd /tmp
git clone https://github.com/linglanfang/RBRP.git
sudo cp RBRP/scripts/*.pl /opt/icSHAPE/scripts/

# パス設定
echo 'export PATH="/opt/icSHAPE/scripts:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

## セットアップ手順

### 1. リポジトリクローン
```bash
git clone <repository-url>
cd kool
```

### 2. Python環境セットアップ
```bash
pip install -r requirements.txt
```

### 3. 参照ゲノム準備
```bash
# ヒトゲノム参照配列ダウンロード
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

# 解凍
gunzip *.gz

# 転写産物インデックス作成
gffread -g Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -w transcriptome.fa \
        Homo_sapiens.GRCh38.110.gtf

bowtie2-build transcriptome.fa transcriptome_index
```

### 4. ゲノムサイズファイル作成
```bash
# 染色体サイズファイル生成
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > genome.size
```

## 使用方法

### 1. Jupyter Labの起動
```bash
jupyter lab
```

### 2. メインパイプラインノートブックを開く
`notebooks/RBRP_Dry_Protocol_Pipeline.ipynb`を開きます。

### 3. 設定の編集
ノートブックの設定セクションで以下を指定：

```python
# 入力FASTQファイル
INPUT_FASTQ_FILES = {
    'sample1_probe_only': '/path/to/sample1_probe_only.fastq',
    'sample1_probe_drug': '/path/to/sample1_probe_drug.fastq',
    'sample2_DMSO_ctrl': '/path/to/sample2_DMSO_ctrl.fastq',
    'sample2_drug_ctrl': '/path/to/sample2_drug_ctrl.fastq'
}

# 参照ファイル
REFERENCE_FILES = {
    'genome_fasta': '/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
    'genome_gtf': '/path/to/Homo_sapiens.GRCh38.110.gtf',
    'transcriptome_index': '/path/to/transcriptome_index',
    'adapter_fasta': '/path/to/adapter.fa'
}
```

### 4. パイプライン実行
「Run All Cells」で全自動実行、または各セルを順次実行。

## 出力ファイル

### 処理済みデータ (`data/processed/`)
- `*_demux.fastq`: デマルチプレックス済みFASTQ
- `*_trimmed.fastq`: トリミング済みFASTQ
- `*.sam`: アライメント結果
- `*.rpkm`: 転写産物発現量
- `*.rt`: RTstop解析結果
- `*_rbrp.out`: RBRPスコア

### 最終結果 (`data/results/`)
- `*.bedgraph`: IGV用bedgraphファイル
- `*.bw`: IGV用bigwigファイル
- `processing_summary.csv`: 処理サマリー
- `figures/`: 可視化プロット

## パラメータ設定

### 解析パラメータ
```python
ANALYSIS_PARAMS = {
    'min_read_length': 25,           # 最小リード長
    'min_rpkm': 1,                   # 最小RPKM閾値
    'min_sequencing_depth': 200,     # 最小シーケンス深度
    'rbrp_score_threshold': 0.12,    # RBRPスコア閾値
    'p_value_threshold': 0.05        # p値閾値
}
```

### バーコード設定
```python
BARCODES = {
    'sample1_probe_only': 'GGTT',
    'sample1_probe_drug': 'TTGT', 
    'sample2_DMSO_ctrl': 'ACCT',
    'sample2_drug_ctrl': 'CAAT'
}
```

## 結果の解釈

### RBRPスコア
- **範囲**: -1 ～ 1
- **高スコア (≥ 0.12)**: RNA-薬物相互作用の可能性が高い部位
- **低スコア (< 0.12)**: 相互作用の可能性が低い部位

### 統計的有意性
- **p値 ≤ 0.05**: 統計的に有意な相互作用
- **シーケンス深度 ≥ 200**: 信頼性の高い結果

### IGVでの可視化
1. IGVを起動
2. ヒトゲノム (GRCh38) を選択
3. 生成された `.bw` ファイルを読み込み
4. 興味のある遺伝子領域にズーム

## トラブルシューティング

### よくある問題

**1. 外部ツールが見つからない**
```bash
# ツールの存在確認
which fastqc bowtie2 cutadapt
```

**2. メモリ不足エラー**
- 大きなファイルの場合、計算リソースを増やす
- サンプル数を減らして分割処理

**3. アライメント率が低い**
- 参照転写産物インデックスの確認
- FASTQ品質の確認
- アダプター配列の確認

**4. RBRPスコアが計算されない**
- 最小RPKM閾値を下げる
- 最小シーケンス深度を下げる

### ログファイル
詳細なエラー情報は `logs/` ディレクトリのログファイルを確認してください。

## システム要件

### 最小要件
- **CPU**: 4コア以上
- **メモリ**: 16GB以上
- **ストレージ**: サンプルあたり200GB以上

### 推奨要件
- **CPU**: 8コア以上
- **メモリ**: 32GB以上
- **ストレージ**: SSD、サンプルあたり500GB以上

## 処理時間の目安

| ステップ | 処理時間（1サンプル） |
|---------|-------------------|
| デマルチプレックス | 5-10分 |
| トリミング | 10-20分 |
| アライメント | 30-60分 |
| RBRPスコア計算 | 20-30分 |
| 可視化ファイル生成 | 5-10分 |
| **合計** | **約1-2時間** |

*処理時間はファイルサイズと計算リソースに依存

## 引用

このパイプラインを使用した場合は、以下の論文を引用してください：

```bibtex
@article{fang2023rbrp,
  title={Reactivity-based RNA profiling for analyzing transcriptome interactions of small molecules in human cells},
  author={Fang, Linglan and Kool, Eric T},
  journal={STAR Protocols},
  volume={4},
  number={4},
  pages={102670},
  year={2023},
  publisher={Elsevier}
}
```

## ライセンス

このプロジェクトはMITライセンスの下で提供されます。

## サポート

- **Issue報告**: GitHub Issues
- **質問**: [連絡先メール]
- **Wiki**: [プロジェクトWiki]

## 貢献

プロジェクトへの貢献を歓迎します：

1. Forkしてください
2. 新機能ブランチを作成 (`git checkout -b feature/AmazingFeature`)
3. 変更をコミット (`git commit -m 'Add some AmazingFeature'`)
4. ブランチにプッシュ (`git push origin feature/AmazingFeature`)
5. Pull Requestを作成

---

**注意**: このパイプラインはバイオインフォマティクス初学者向けに設計されていますが、基本的なLinuxコマンドラインとPythonの知識があることを前提としています。