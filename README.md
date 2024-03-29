# GenReads
Generate Fake Reads



## 引物测序示意图
假设使用515F/806R引物。双端测序数据对应图示中的R1与R2,三代测序如图示中Pacbio_R。
```
                                     --t_r_R2---->
   -----R1----->                     <-----R2-----
 --*---------------------------------------------*------
  515                                           806
   ------------------Pacbio_R-------------------->

   |--read_len--|
   |---------------insert_size--------------------|
```


## cfg.ini说明
```
[common]
output_folder = output/
each_contig_as_a_genome = True    ## 是否按照gap=insert_size-1链接同一个文件中的contigs
                                  ## False:each_file_as_a_genome   True:each_contig_as_a_genome

[input]                 ## 处理单个输入fasta
input_SampleID = gtdb
input_Amount = 6        ## 为每个contig生成这些数目的fragments，如果contigs已被链接，则算作1个
input_RawFasta = raw/gtdb.all.16s.ncbi-genus.fasta


[inputs]                 ## 处理多个输入fasta
input_cfg = config_Random.xls


[Random]                 ## 随机切割基因组所需参数
read_len = 250
error_rate = 0.01
insert_size = 450


[Primer]                 ## 按primer切割基因组所需参数
read_len = 250
error_rate = 0.01
forward=341F
reverse=806R
insert_size=466
insert_size_max=476
insert_size_min=300
size_punish=0.02
f_treash = 3
r_treash = 3

[Anchor]                  ## TBA
```


config_Random.xls（示例）
```
Sample  Amount  RawGenome
Spe1      1X      raw/GCF_002871005.1_ASM287100v1_genomic.fna
Spe2      50M     raw/GCF_007829915.1_ASM782991v1_genomic.fna
Spe3      2X      raw/GCF_009687845.1_ASM968784v1_genomic.fna


* 以tab分隔
* Amount列中，'X'代表测序深度，B K M G T 以 1024 阶乘；无单位代表reads count；字母大小写不限。
* Amount指为每个contig生成这些数目的fragments，如果contigs已被链接，则算作1个
```

## 使用说明

```
from GenReads import *

cfg_file = 'cfg.ini'
MainCall(cfg_file)
```

自动根据 cfg.ini 中 ```forward``` 与 ```reverse``` 的存在选择模式：

| 模式 | 说明 | -- |
| -- | -- | -- |
| RandomModeGeneration | 使用cfg.ini中```[Random]```设定 | 随机切割，模仿超声波打碎模式 |
| PrimerModeGeneration | 使用cfg.ini中```[Primer]```设定 | 切取F/R引物中间序列，选取Levenshtein小于f/r_treash且两引物直接距离合理者 |
| AnchorModeGeneration | -- | 单侧引物，TBA |


对于引物模式其实blast然后提取序列会更快速，但个人电脑不想装blast。






