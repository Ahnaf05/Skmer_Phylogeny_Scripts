mkdir ref_dir_
python3 N_remover.py

rm -r ref_dir
mv ref_dir_ ref_dir

skmer reference ref_dir

mv ref-dist-mat.txt Dist_Matrix.txt

tail -n +2 "Dist_Matrix.txt" > "Dist_Matrix.tmp" && mv "Dist_Matrix.tmp" "Dist_Matrix.txt"

sed -i '1 i\14' Dist_Matrix.txt

fastme -i Dist_Matrix.txt -o Skmer.tre

