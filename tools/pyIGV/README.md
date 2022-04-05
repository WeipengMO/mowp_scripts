<!--
 * @Date         : 2021-08-30 14:43:43
 * @LastEditTime : 2021-08-30 15:03:42
 * @LastEditors  : windz
 * @FilePath     : /tools/pyIGV/README.md
-->

# pyIGV

## What's New

### 2021.08.30

- add [genome browser tracks](./coverage_plot.ipynb) for coverage plot.

### 2021.05.25

- add read_color;bam_title parameter to IGV.add_bam

### 2021.05.17

- add lambda_for_parse_geneId|auto_create_gz to IGV.add_gene_model. (Zhijian)
  - can provide a lambda expression to add_gene_model. It will be used to parse the gene_id column. 
  - add_gene_model will first detect the suffix of anno path. If is .bed, will detect whether corresponding .gz is exists. if not, will automatic build this file and build index. 

### 2021.05.15

- Add type hints to partial function. (Zhijian)
- Add force_tag_check variable to IGV. Related initialization of IGV. If True, reads which don't all have required tags will be ignored. If False, these reads will still be plotted.
- Allow users to customize TAGs. Related initialization of IGV.

