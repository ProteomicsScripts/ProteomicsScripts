<h2>Python script to modify descriptions in fasta file. </h2>

<h3> Examples </h3>
For example, the following descriptions:

* `CANCER_sp_Q969R8,ITFG2_HUMAN Integrin-alpha FG-GAP repeat-containing protein 2|c.560C>T|p.S187L|DSLLVTL|Missense|COSMIC|Large intestine(1)`

* `CANCER_sp_O75369,FLNB_HUMAN Filamin-B|c.3631G>A|p.E1211K|YGGKLVP|Missense|COSMIC|Endometrium(1) OR sp|O75369-8,FLNB_HUMAN Isoform 8 of Filamin-B|c.3631G>A|p.E1211K|YGGKLVP|Missense|COSMIC|Endometrium(1)`

will be modified to:
* `CANCER_sp|Q969R8_c.560C>T|ITFG2_HUMAN Integrin-alpha FG-GAP repeat-containing protein 2,,p.S187L,DSLLVTL,Missense,COSMIC,Large intestine(1)`
* `CANCER_sp|O75369_c.3631G>A|FLNB_HUMAN Filamin-B,,p.E1211K,YGGKLVP,Missense,COSMIC,Endometrium(1)`
* `sp|O75369-8_c.3631G>A|FLNB_HUMAN Isoform 8 of Filamin-B,p.E1211K,YGGKLVP,Missense,COSMIC,Endometrium(1)`

Note that the second original description contained ` OR `, and it is split into two descriptions in the modified file (with the same corresponding sequence).

<h3> Example Usage </h3> 

`python modify_database.py example.fasta`

The above command will generate a new fasta file with "_mod" in the file name, i.e. `example_mod.fasta` for the above example
