model_dir = /directory to downloaded models/
foldseek_dir = /directory to foldseek/

log_dir = "./logs"
test -d ${log_dir} || mkdir -p ${log_dir}

for i in 101;do
python   03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --model_name ${model_dir}/esm2_t33_650M_UR50D --batch_size 16 --num_epochs 20 --top_layers 4 --method top_layer --learning_rate 0.0001 --use_listmle --log_path ${log_dir}/03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.ESM2_epoch20_layer4_listMLE_lr0.0001.seed${i}.avGFP.log

python  03.Prime_690M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --model_name ${model_dir}/Prime_690M --batch_size 16 --num_epochs 20 --top_layers 4 --method top_layer --learning_rate 0.0001 --use_listmle --log_path ${log_dir}/03.Prime_690M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.Prime_epoch20_layer4_listMLE_lr0.0001.seed${i}.avGFP.log

python 03.SaProt_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --foldseek_path ${foldseek_dir}/foldseek/bin/foldseek --foldseek_util_path ${model_dir}/SaProt/utils/foldseek_util.py --pdb_path 02.AF3_avGFP.cif --chain_id A --model_name ${model_dir}/SaProt_650M_PDB --batch_size 16 --num_epochs 20 --top_layers 4 --method top_layer --learning_rate 0.0001 --use_listmle --log_path ${log_dir}/03.SaProt_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.SaProt_epoch20_layer4_listMLE_lr0.0001.seed${i}.avGFP.log

python   03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --model_name ${model_dir}/esm2_t33_650M_UR50D --batch_size 5 --num_epochs 20 --top_layers 4 --method top_layer --learning_rate 0.0001  --log_path ${log_dir}/03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.ESM2_epoch20_layer4_MSEloss_lr0.0001.seed${i}.avGFP.log

python  03.Prime_690M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --model_name ${model_dir}/Prime_690M --batch_size 5 --num_epochs 20 --top_layers 4 --method top_layer --learning_rate 0.0001  --log_path ${log_dir}/03.Prime_690M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.Prime_epoch20_layer4_MSEloss_lr0.0001.seed${i}.avGFP.log

python 03.SaProt_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --foldseek_path ${foldseek_dir}/foldseek/bin/foldseek --foldseek_util_path ${model_dir}/SaProt/utils/foldseek_util.py --pdb_path 02.AF3_avGFP.cif --chain_id A --model_name ${model_dir}/SaProt_650M_PDB --batch_size 5 --num_epochs 20 --top_layers 4 --method top_layer --learning_rate 0.0001 --log_path ${log_dir}/03.SaProt_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.SaProt_epoch20_layer4_MSEloss_lr0.0001.seed${i}.avGFP.log



python   03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --model_name ${model_dir}/esm2_t33_650M_UR50D --batch_size 16 --num_epochs 20 --rank_count 4 --method lora --learning_rate 0.0001 --use_listmle --log_path ${log_dir}/03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.ESM2_epoch20_lora_listMLE_lr0.0001.seed${i}.avGFP.log

python  03.Prime_690M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --model_name ${model_dir}/Prime_690M --batch_size 16 --num_epochs 20 --rank_count 4 --method lora --learning_rate 0.0001 --use_listmle --log_path ${log_dir}/03.Prime_690M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.Prime_epoch20_lora_listMLE_lr0.0001.seed${i}.avGFP.log

python 03.SaProt_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --foldseek_path ${foldseek_dir}/foldseek/bin/foldseek --foldseek_util_path ${model_dir}/SaProt/utils/foldseek_util.py --pdb_path 02.AF3_avGFP.cif --chain_id A --model_name ${model_dir}/SaProt_650M_PDB --batch_size 16 --num_epochs 20 --rank_count 4 --method lora --learning_rate 0.0001 --use_listmle --log_path ${log_dir}/03.SaProt_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.SaProt_epoch20_lora_listMLE_lr0.0001.seed${i}.avGFP.log

python   03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --model_name ${model_dir}/esm2_t33_650M_UR50D --batch_size 5 --num_epochs 20 --rank_count 4 --method lora --learning_rate 0.0001  --log_path ${log_dir}/03.ESM2_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.ESM2_epoch20_lora_MSEloss_lr0.0001.seed${i}.avGFP.log

python  03.Prime_690M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --model_name ${model_dir}/Prime_690M --batch_size 5 --num_epochs 20 --rank_count 4 --method lora --learning_rate 0.0001  --log_path ${log_dir}/03.Prime_690M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.Prime_epoch20_lora_MSEloss_lr0.0001.seed${i}.avGFP.log

python 03.SaProt_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.py --seed ${i} --csv_path 01.avGFP_fluorescence_30202sequence.csv --foldseek_path ${foldseek_dir}/foldseek/bin/foldseek --foldseek_util_path ${model_dir}/SaProt/utils/foldseek_util.py --pdb_path 02.AF3_avGFP.cif --chain_id A --model_name ${model_dir}/SaProt_650M_PDB --batch_size 5 --num_epochs 20 --rank_count 4 --method lora --learning_rate 0.0001 --log_path ${log_dir}/03.SaProt_650M_FineTunning_regression_manual_assembled.predictions_TrainValidationTest.SaProt_epoch20_lora_MSEloss_lr0.0001.seed${i}.avGFP.log
done

