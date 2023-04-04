module load miniconda/4.9.2
source activate qiime2-2022.2

###分类器
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 99_otus.fasta \
  --output-path 99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 99_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 99_otus.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza

# fastp
mkdir fastp
mkdir data/seq
awk '{system("fastp -i data/rawdata/"$1"_1.fq.gz \
                    -I data/rawdata/"$1"_2.fq.gz \
                    -o data/seq/"$1"_1.fq.gz \
                    -O data/seq/"$1"_2.fq.gz \
                    -j fastp/"$1".json \
                    -h fastp/"$1".html")}' \
                    <(tail -n+2 metadata.txt) 
# 根据metadata生成manifest文件, 双端数据
awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} \
  NR>1{print $1"\t$PWD/data/seq/"$1"_1.fq.gz\t$PWD/data/seq/"$1"_2.fq.gz"}' \
  metadata.txt > manifest
  
# 根据metadata生成manifest文件, 单端/join 数据
# awk -F '\r' 'NR==1{print "sample-id\tabsolute-filepath"} \
#   NR>1{print $1"\t$PWD/rawdata/${project_name}/"$1".fastq.gz"}' \
#   rawdata/${project_name}/sample.txt > manifest


# 数据导入qiime2，格式为双端33格式
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
  
# 数据导入qiime2，格式为单端33格式
# qiime tools import \
#   --type 'SampleData[SequencesWithQuality]' \
#   --input-path manifest \
#   --output-path single-end-demux.qza \
#   --input-format SingleEndFastqManifestPhred33V2

### 2. 生成特征表和代表序列

#序列切除引物和质量控制
#切除引物
#双端数据去除引物
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences paired-end-demux.qza \
        --p-cores 8 \
        --p-front-f AACMGGATTAGATACCCKG \
		    --p-front-r ACGTCATCCCCACCTTCC \
        --o-trimmed-sequences trimmed-seqs.qza \
        --verbose
#单端数据去除引物
# qiime cutadapt trim-paired 
#         --i-demultiplexed-sequences single-end-demux.qza 
#         --p-cores 8 --p-front 正向引物碱基 
#         --o-trimmed-sequences trimmed-seqs.qza  
#         --verbose

#过滤
qiime quality-filter q-score \
  --i-demux trimmed-seqs.qza \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza


#查看过滤后的数据
qiime demux summarize \
  --i-data demux-filtered.qza \
  --o-visualization demux-filtered.qzv

##### Deblur，--p-trim-length n which truncates the sequences at position n. 
##### In general, the Deblur developers recommend setting this value to a length 
##### where the median quality score begins to drop too low.

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 150 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --p-jobs-to-start 16 \
  --o-stats deblur-stats.qza
  
qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv


#### 导入主要分析流程
cp table-deblur.qza table.qza
cp rep-seqs-deblur.qza rep-seqs.qza


### alpha-rarefaction
qiime diversity alpha-rarefaction --i-table table.qza \
                                  --p-max-depth 10000 \
                                  --p-min-depth 100 \
                                  --o-visualization alpha-rarefaction.qzv

### 3. Alpha和beta多样性分析
# 构建进化树用于多样性分析
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
# 采样深度通常选择最小值，来自table.qzv
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 5000 \
  --m-metadata-file metadata.txt \
  --output-dir core-metrics-results

### 4. 物种注释
# classifier_gg_13_8_99.qza
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification taxonomy.qza

### 5. 功能注释
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path rep_seq

qiime tools export \
  --input-path table.qza \
  --output-path feature

source activate picrust2
picrust2_pipeline.py -s rep_seq/dna-sequences.fasta \
    -i feature/feature-table.biom \
    -o picrust2_out_pipeline \
    -p 8

# EC
add_descriptions.py -i picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    -m EC \
                    -o picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
gunzip picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
# KO 
add_descriptions.py -i picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    -m KO \
                    -o picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
gunzip picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
# MetaCyc pathway
add_descriptions.py -i picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz \
                    -m METACYC \
                    -o picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz
gunzip picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz

