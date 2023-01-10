# GenReads
Generate Fake Reads



## 引物测序示意图
假设使用515F/806R引物。双端测序数据对应图示中的R1与t_r_R2,三代测序如图示中Pacbio_R。
```
                                    ---t_r_R2---->
   -----R1---->                     <-----R2------
 --*---------------------------------------------*------
  515                                           806
   ------------------Pacbio_R-------------------->
```


## 对多个输入文件指定指定输入与输出（示例）
```
Sample  Amount  RawGenome
Spe1      1X      raw/GCF_002871005.1_ASM287100v1_genomic.fna
Spe2      50M     raw/GCF_007829915.1_ASM782991v1_genomic.fna
Spe3      2X      raw/GCF_009687845.1_ASM968784v1_genomic.fna


* 以tab分隔
* Amount列中，'X'代表测序深度，B K M G T 以 1024 阶乘；无单位代表B(base pair)；字母大小写不限。
```


## 模式一、随机切割基因组






