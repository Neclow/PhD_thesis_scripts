# Scripts used for figures in my thesis and defence

## Install dependencies

```bash
pixi install
pixi run post_install
```

## Fusexin sequence similarity

```bash
pixi run clustalo -i data/fusexins2.fa --distmat-out=outputs/fusexins_similarity.txt --percent-id
```

## Animal phylogeny

```bash
pixi run iqtree -s data/animals.phy -m TIM2+I+G -g data/animals.constr0 --prefix outputs/animals.phy -B 1000
```
