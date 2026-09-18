# 使用自己的群体数据快速开始

该入口适用于已经拥有群体 PLINK 文件、对应参考基因组以及一个或多个关注区域的用户。

## 必需输入

| 输入 | 格式要求 |
|---|---|
| 基因型 | 同一前缀的 `PREFIX.bed`、`PREFIX.bim`、`PREFIX.fam` |
| 参考基因组 | 未压缩 FASTA；序列名称必须与 PLINK 染色体名称一致 |
| 加权区域 | BED3/BED4，0-based、half-open；每个文件设置一个权重 |
| 设计窗口 | BED3/BED4，0-based、half-open；可由 `make_windows.py` 生成 |

PLINK `.bed` 是二进制基因型文件，与窗口和加权区域所用的文本 BED3/BED4
区间文件不是同一种格式。

原始 BIM 可以使用 rsID，也可以尚未确定等位基因方向。准备工具会让 PLINK2
根据 FASTA 确定 REF/ALT，仅保留二等位 A/C/G/T SNP，排序后将 ID 统一写成
`CHROM:POS:REF:ALT`。工具会生成一套新的 PLINK 文件，不会原位修改源文件。

FASTA、BIM、设计窗口和所有区域 BED 必须来自同一组装版本，并使用相同染色体名称。

## 1. 规范化 PLINK 文件

```bash
python tools/prepare_plink_input.py \
  --bfile data/my_population \
  --reference data/assembly.fa \
  --out data/prepared/my_population \
  --chr-set 31 --threads 8
```

如果 PLINK 与 FASTA 的染色体名称不同，可用 `--chrom-map` 提供无表头的两列表
`OLD_CHROM NEW_CHROM`；未列出的名称保持不变。

`--chr-set` 应设置为该物种在 PLINK 中使用的常染色体数，例如人 22、马 31。
输出包括新的 BED/BIM/FAM 和准备汇总文件 `my_population.prepare.json`。本 SNP 芯片
流程会主动剔除 indel、多等位位点和非 A/C/G/T 等位基因。

## 2. 生成设计窗口

```bash
python tools/make_windows.py \
  --reference data/assembly.fa \
  --output data/design_windows.bed \
  --window-size 40000 --step 40000 \
  --exclude-regex '^(MT|chrM)$'
```

若 `assembly.fa.fai` 不存在，工具会自动生成。默认保留每条序列末端不足完整宽度的
窗口，还可使用 `--min-contig-length`、`--min-window-size`、`--include-regex`、
`--exclude-regex` 或重复的 `--exclude-contig` 进行筛选。

设计窗口决定在哪些范围内选择候选 SNP；关注区域 BED 只负责给位点加权，不能替代
设计窗口。

## 3. 准备关注区域与权重

每个文件使用 BED3 或 BED4：

```text
chr1    99999    120000    altitude_region_1
chr1    450000   470000    altitude_region_2
```

起点为 0-based 且包含，终点不包含。空行和以 `#` 开头的行会被忽略。在 YAML 中给
每个文件设置一个权重：

```yaml
scoring:
  profile: custom
  include_maf: true
  include_vep: false
  region_aggregation: sum
  weighted_regions:
    - name: altitude
      path: data/altitude.bed
      weight: 25
    - name: immune
      path: data/immune.bed
      weight: 12
```

权重可以为正数、0 或负数。同一 SNP 命中多个文件时，`sum` 累加全部权重，`max`
取数值最大的权重，`first` 使用 YAML 中第一个命中的文件。

## 4. 选择染色体、预检并运行

复制 `config/custom_regions.example.yaml`，替换路径和物种参数。第一次真实数据测试建议：

```yaml
target_chromosomes:
  - ECA1
```

随后预检并运行到评分断点：

```bash
python tools/preflight.py \
  --config config/my_population.yaml --base-dir . \
  --report preflight.json

snakemake -s workflow/Snakefile \
  --configfile config/my_population.yaml \
  --cores 8 --software-deployment-method conda \
  run_scoring
```

预检会检查 PLINK 二进制结构与记录数、规范 ID、抽样 REF 与 FASTA 是否一致、FAI
边界、BED 坐标语义、染色体名称、区域权重及所有必需文件。错误会阻止运行，警告则
提示用户复核。

检查 `results/02_scoring/selection_ld_tag_set_candidates.tsv` 后，可改为
`target_chromosomes: all`，分配足够核心并显式调用 `run_full`。流程只建立一套共享
Bowtie2/Jellyfish 索引，再按染色体并行执行 OligoMiner。三个断点命令见根目录 README。

## 小型 Demo

运行 `bash demo/run_demo.sh`，即可测试一个包含 6 个 SNP、4 个样本和两个加权区域的
完整小例子。详细说明见 [`demo/README.md`](../demo/README.md)。
