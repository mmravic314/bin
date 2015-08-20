#!bin/bash
#$ -S /bin/bash                    
#$ -o .                       
#$ -e .                        
#$ -cwd                           
#$ -r y                           
#$ -j y                           
#$ -l mem_free=1G                  
#$ -l arch=linux-x64               
#$ -l netapp=1G,scratch=1G         
#$ -l h_rt=24:00:00
#$ -t 1-151

#python /netapp/home/mmravic/bin/runSadicVolume.py

tasks=(0 1jb0 1kf6 1kqf 1lgh 1nek 1nkz 1orq 1ors 1ots 1p49 1q16 1rh5 1u7g 1yq3 2bhw 2bl2 2bs2 2cfq 2dyr 2hyd 2j58 2j7a 2jln 2nq2 2onk 2qks 2qts 2r9r 2uuh 2vpz 2wsw 2x2v 2xfn 2xq2 2xqu 2xtv 2z73 2zjs 2zxe 3ayf 3b9y 3beh 3cx5 3d9s 3gia 3h90 3h9v 3jqo 3k3f 3kly 3lbw 3ldc 3m73 3mk7 3n5k 3ne5 3o7q 3ob6 3odu 3ouf 3puw 3qe7 3rko 3rlb 3rqw 3rvy 3s8g 3tdp 3tij 3tlw 3tt1 3tui 3tx3 3ug9 3v5u 3vma 3wbn 3wdo 3wmg 3wmm 3wo7 3wu2 3zjz 3zk1 4a01 4bbj 4bpm 4bwz 4c7r 4c9h 4cad 4d2e 4dji 4dve 4dx5 4ev6 4ezc 4g7v 4gc0 4gx0 4h33 4huq 4i0u 4jkv 4jr9 4k5y 4kjs 4knf 4kpp 4ky0 4lds 4lep 4lp8 4lz6 4m48 4mrs 4n7w 4o6m 4o6y 4o9p 4or2 4p02 4p79 4pgr 4phz 4pl0 4pyp 4qnd 4qtn 4quv 4rdq 4rng 4tnw 4tq4 4twk 4u9n 4umw 4us3 4wd8 4wgvS)
input="${tasks[$SGE_TASK_ID]}"

/netapp/home/mmravic/bin/confind --p $input

echo $input
