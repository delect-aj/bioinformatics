# 设置工作目录
wd=/home/dongbiao/16S/oral_disease

# 进入工作目录
cd ${wd}
# 激活QIIME2工作环境
module load miniconda/4.9.2
source activate qiime2-2022.2

###分类器
# qiime tools import \
  # --type 'FeatureData[Sequence]' \
  # --input-path reference/99_otus.fasta \
  # --output-path reference/99_otus.qza
# 
# qiime tools import \
#   --type 'FeatureData[Taxonomy]' \
#   --input-format HeaderlessTSVTaxonomyFormat \
#   --input-path 99_otu_taxonomy.txt \
#   --output-path ref-taxonomy.qza
#   
# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads 99_otus.qza \
#   --i-reference-taxonomy ref-taxonomy.qza \
#   --o-classifier classifier.qza

### 项目名称
export project_name=PRJNA666891
mkdir ${project_name}
# 根据metadata生成manifest文件, 双端数据
awk -F '\r' 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} \
  NR>1{print $1"\t$PWD/rawdata/${project_name}/clean_data/"$1"_1.fastq.gz\t$PWD/rawdata/${project_name}/clean_data/"$1"_2.fastq.gz"}' \
  rawdata/${project_name}/sample.txt > manifest
  
# 根据metadata生成manifest文件, 单端/join 数据
awk -F '\r' 'NR==1{print "sample-id\tabsolute-filepath"} \
  NR>1{print $1"\t$PWD/rawdata/${project_name}/"$1".fastq.gz"}' \
  rawdata/${project_name}/sample.txt > manifest


# 数据导入qiime2，格式为双端33格式
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest \
  --output-path ${project_name}/paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
  
# 数据导入qiime2，格式为单端33格式
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest \
  --output-path ${project_name}/single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

### 2. 生成特征表和代表序列

#序列切除引物和质量控制
#切除引物
#双端数据去除引物
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences ${project_name}/paired-end-demux.qza \
        --p-cores 8 \
        --p-front-f GTGCCAGCAGCCGCGGTAA \
		--p-front-r GTGCCAGCAGCCGCGGTAA \
        --o-trimmed-sequences ${project_name}/trimmed-seqs.qza \
        --verbose
#单端数据去除引物
qiime cutadapt trim-paired 
        --i-demultiplexed-sequences ${project_name}/single-end-demux.qza 
        --p-cores 8 --p-front 正向引物碱基 
        --o-trimmed-sequences ${project_name}/trimmed-seqs.qza  
        --verbose

#过滤
qiime quality-filter q-score \
  --i-demux ${project_name}/single-end-demux.qza \
  --o-filtered-sequences ${project_name}/demux-filtered.qza \
  --o-filter-stats ${project_name}/demux-filter-stats.qza


#查看过滤后的数据
qiime demux summarize \
  --i-data ${project_name}/demux-filtered.qza \
  --o-visualization ${project_name}/demux-filtered.qzv

##### Deblur，--p-trim-length n which truncates the sequences at position n. 
##### In general, the Deblur developers recommend setting this value to a length 
##### where the median quality score begins to drop too low.

qiime deblur denoise-16S \
  --i-demultiplexed-seqs ${project_name}/demux-filtered.qza \
  --p-trim-length 220 \
  --p-left-trim-len 20 \
  --o-representative-sequences ${project_name}/rep-seqs-deblur.qza \
  --o-table ${project_name}/table-deblur.qza \
  --p-sample-stats \
  --p-jobs-to-start 16 \
  --o-stats ${project_name}/deblur-stats.qza
  
qiime deblur visualize-stats \
  --i-deblur-stats ${project_name}/deblur-stats.qza \
  --o-visualization ${project_name}/deblur-stats.qzv


#### 导入主要分析流程
cp ${project_name}/table-deblur.qza ${project_name}/table.qza
cp ${project_name}/rep-seqs-deblur.qza ${project_name}/rep-seqs.qza

### meta analysis
qiime feature-table summarize \
  --i-table ${project_name}/table.qza \
  --o-visualization ${project_name}/table.qzv

### closed-reference
qiime vsearch cluster-features-closed-reference \
    --i-sequences ${project_name}/rep-seqs.qza \
    --i-table ${project_name}/table.qza \
    --i-reference-sequences reference/99_otus.qza \
    --p-perc-identity 0.97 --p-threads 12 \
    --o-clustered-table ${project_name}/clustered-table \
    --o-clustered-sequences ${project_name}/clustered-sequences \
    --o-unmatched-sequences ${project_name}/unmatched-sequences

qiime feature-table summarize \
  --i-table ${project_name}/clustered-table.qza \
  --o-visualization ${project_name}/clustered-table.qzv

export project_name=merge
#### 导入主要分析流程
cp ${project_name}/clustered-table.qza ${project_name}/table.qza
cp ${project_name}/clustered-sequences.qza ${project_name}/rep-seqs.qza

### 合并数据集
### 合并table
qiime feature-table merge --i-tables OEP000837/table.qza PRJNA386665/table.qza \
                                     PRJNA700849/table.qza PRJEB37501/table.qza \
                                     PRJNA412445/table.qza PRJNA751046/table.qza \
                                     PRJEB39064/table.qza PRJNA421234/table.qza \
                                     PRJNA756784/table.qza PRJNA666891/table.qza \
                                     PRJNA744870/table.qza \
                          --o-merged-table merge/table.qza
### 合并代表序列
qiime feature-table merge-seqs --i-data OEP000837/rep-seqs.qza PRJNA386665/rep-seqs.qza \
                                        PRJNA700849/rep-seqs.qza PRJEB37501/rep-seqs.qza \
                                        PRJNA412445/rep-seqs.qza PRJNA751046/rep-seqs.qza \
                                        PRJEB39064/rep-seqs.qza PRJNA421234/rep-seqs.qza \
                                        PRJNA756784/rep-seqs.qza PRJNA666891/rep-seqs.qza \
                                        PRJNA744870/rep-seqs.qza \
                                --o-merged-data merge/rep-seqs.qza

qiime feature-table summarize \
  --i-table ${project_name}/table.qza \
  --o-visualization ${project_name}/table.qzv

### alpha-rarefaction
qiime diversity alpha-rarefaction --i-table ${project_name}/table.qza \
                                  --p-max-depth 10000 \
                                  --p-min-depth 100 \
                                  --o-visualization ${project_name}/alpha-rarefaction.qzv

### 3. Alpha和beta多样性分析
# 构建进化树用于多样性分析
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${project_name}/rep-seqs.qza \
  --o-alignment ${project_name}/aligned-rep-seqs.qza \
  --o-masked-alignment ${project_name}/masked-aligned-rep-seqs.qza \
  --o-tree ${project_name}/unrooted-tree.qza \
  --o-rooted-tree ${project_name}/rooted-tree.qza
# 采样深度通常选择最小值，来自table.qzv
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ${project_name}/rooted-tree.qza \
  --i-table ${project_name}/table.qza \
  --p-sampling-depth 2000 \
  --m-metadata-file metadata.txt \
  --output-dir ${project_name}/core-metrics-results

# Alpha多样性
# 可选的alpha指数有 faith_pd、shannon、observed_features、evenness
alpha=(faith_pd shannon observed_features evenness)
for index in ${alpha[*]}; do
    qiime diversity alpha-group-significance \
      --i-alpha-diversity ${project_name}/core-metrics-results/${index}_vector.qza \
      --m-metadata-file metadata.txt \
      --o-visualization ${project_name}/core-metrics-results/${index}-group-significance.qzv
    # 导出alpha多样性
    qiime tools export \
      --input-path ${project_name}/core-metrics-results/${index}_vector.qza \
      --output-path ${project_name}/alpha/${index}
done

# Beta多样性组间显著性分析和可视化
# 可选的beta指数有 unweighted_unifrac、bray_curtis、weighted_unifrac和jaccard
# 指定分组是减少计算量，置换检验较耗时
distance=jaccard
column=clinical_stage
qiime diversity beta-group-significance \
  --i-distance-matrix ${project_name}/core-metrics-results/${distance}_distance_matrix.qza \
  --m-metadata-file metadata1.txt \
  --m-metadata-column ${column} \
  --o-visualization ${project_name}/core-metrics-results/${distance}-${column}-significance.qzv \
  --p-pairwise


# 导出beta多样性表格
beta=(unweighted_unifrac bray_curtis weighted_unifrac jaccard)
for i in ${beta[*]}; do
    qiime tools export \
      --input-path ${project_name}/core-metrics-results/${i}_distance_matrix.qza \
      --output-path ${project_name}/beta/${i}
done

# 导出feature table
qiime tools export \
  --input-path ${project_name}/table.qza \
  --output-path ${project_name}/feature

### 4. 物种注释


# classifier_gg_13_8_99.qza
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads ${project_name}/rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification ${project_name}/taxonomy.qza
#导出物种注释分类水平
level=(2 3 4 5 6)
for i in ${level[*]}; do
	qiime taxa collapse \
	  --i-table ${project_name}/table.qza \
	  --i-taxonomy ${project_name}/taxonomy.qza \
	  --p-level $((i)) \
	  --o-collapsed-table ${project_name}/taxonomy/table-${i}.qza
	# 格式化特征表，添加伪计数
	qiime composition add-pseudocount \
	  --i-table ${project_name}/taxonomy/table-${i}.qza \
	  --o-composition-table ${project_name}/taxonomy/comp-table-${i}.qza
	# 导出物种注释
	qiime tools export \
	  --input-path ${project_name}/taxonomy/comp-table-${i}.qza \
	  --output-path ${project_name}/taxonomy/level_${i}
done
