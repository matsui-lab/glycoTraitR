# glycoTraitR 機能概要

**glycoTraitR** は、pGlyco3 や Glyco-Decipher などの検索エンジンから出力される GPSM (Glycopeptide-Spectrum Match) データを用いて、**グリカン（糖鎖）の構造的特徴（Traits）を抽出、集計、比較分析**するための R パッケージです。

## 1. 設計思想 (Design)

本パッケージは、質量分析データから得られる複雑な糖鎖構造情報を、定量分析可能な「特徴量（Trait）」に変換し、統計解析を容易にすることを目的としています。

*   **入力の柔軟性**: 主要な糖鎖同定ツール（pGlyco3, Glyco-Decipher）の出力をサポートし、統一的なフォーマット（Proteins, Peptide, GlycanStructure, File, Count）に変換して扱います。
*   **構造のグラフ化**: 糖鎖文字列（WURCS 2.0 や pGlyco3 形式）を解析し、`igraph` ベースのツリー構造に変換することで、トポロジーに基づいた高度な特徴抽出を可能にしています。
*   **標準データ形式の採用**: 解析の中核となるデータ構造に Bioconductor の `SummarizedExperiment` クラスを採用しており、既存のオミクス解析エコシステムとの親和性を高めています。
*   **拡張性**: 組み込みの特徴量だけでなく、ユーザーが定義した任意の糖鎖モチーフ（部分構造）を検出・定量する機能を備えています。

## 2. 主な機能 (Features)

### A. データインポートとパース
*   **GPSM データの読み込み**:
    *   `read_pGlyco3_gpsm()`: pGlyco3 の結果ファイルを読み込み、Protein/Peptide/Glycan/File ごとのスペクトル数（Count）を集計します。
    *   `read_decipher_gpsm()`: Glyco-Decipher の結果ファイルを読み込み、WURCS ID を構造文字列に変換して集計します。
*   **糖鎖構造のパース**:
    *   WURCS 2.0 文字列および pGlyco3 形式の文字列を、ノード（単糖）とエッジ（結合）を持つツリー構造に変換します。

### B. 特徴量（Trait）の計算
グリカン構造から以下の3種類のカテゴリの特徴量を計算します。

1.  **組成特徴 (Composition Traits)**:
    *   単糖の数（Hexose, HexNAc, Neu5Ac, Neu5Gc, Fucose）および糖鎖全体のサイズ。
2.  **構造特徴 (Structural Traits) [組み込み]**:
    *   **Antennas**: アンテナ数（枝分かれの数）。
    *   **Types**: High-Mannose, Hybrid, Complex などの型判定。
    *   **Fucosylation**: コア・フコース（Core-Fuc）およびアンテナ・フコース（Ant-Fuc）の有無。
    *   **Bisecting**: バイセクティング GlcNAc の有無。
3.  **ユーザー定義モチーフ (User-defined Motifs)**:
    *   ユーザーが定義した部分構造（例: 特定の分岐パターンや末端構造）をサブグラフ同型性判定（Subgraph Isomorphism）により検索し、出現回数をカウントします。

### C. データ集計と行列構築 (`SummarizedExperiment`)
*   **`build_trait_se()`**: 計算された特徴量を PSM (Peptide-Spectrum Match) ごとに付与し、**Protein** または **Site (Peptide)** レベルで集計します。
*   結果は `SummarizedExperiment` オブジェクトとして出力され、各アッセイ（Assay）が「特徴量 × PSM」の行列となります。

### D. 統計解析と可視化
*   **変動解析 (`analyze_trait_changes()`)**:
    *   2群間（例: Normal vs Symptomatic）で特徴量の変動を検定します。
    *   **Welch's t-test**: 平均値の差の検定。
    *   **Levene's test**: 分散の差（不均一性）の検定。
*   **可視化**:
    *   **`plot_trait_distribution()`**: 特定の特徴量について、群間の分布をヒストグラムとボックスプロットで表示します。
    *   **`plot_glycan_tree()`**: 解析された糖鎖のトポロジー（構造図）を描画します。

## 3. ユースケース (Use Cases)

### ケース1: 疾患バイオマーカー探索
*   **シナリオ**: 健常者と患者の血清糖タンパク質データ（pGlyco3 出力）がある。
*   **アクション**:
    1.  データを読み込み、グリカン特徴量（フコース化、シアル酸数など）を計算。
    2.  `analyze_trait_changes` を用いて、患者群で有意に上昇・減少している構造的特徴を特定する。
    3.  特定の糖鎖修飾（例: コア・フコースの増加）が見られるタンパク質を同定する。

### ケース2: 特定の糖鎖モチーフのスクリーニング
*   **シナリオ**: 特定の抗原エピトープ（例: ルイスX構造など）を持つ糖ペプチドを探したい。
*   **アクション**:
    1.  目的のモチーフ構造（ノードとエッジ）を定義する。
    2.  `build_trait_se` の `motifs` 引数に定義したモチーフを渡し、全データに対してスクリーニングを行う。
    3.  そのモチーフを含む糖ペプチドの分布や量を可視化する。

### ケース3: 糖鎖不均一性（Microheterogeneity）の評価
*   **シナリオ**: 異なる細胞株で作られたバイオ医薬品の糖鎖プロファイルを比較したい。
*   **アクション**:
    1.  サイトレベル（Peptide）でデータを集計。
    2.  各サイトにおける糖鎖の多様性（High-Mannose 型の比率やアンテナ数の分布）を比較し、品質管理や特性解析に役立てる。
