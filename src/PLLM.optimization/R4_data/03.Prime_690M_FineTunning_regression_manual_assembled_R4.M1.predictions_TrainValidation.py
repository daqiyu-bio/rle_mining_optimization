#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------------------#
# Copyright (c) 2026 Beijing Astragenomics Technology Co., Ltd. All rights reserved.        #
# Licensing: see the LICENSE file in the repository root.                                   #
# version : 1.1                                                                             #
#                                                                                           #
# Author:                                                                                   #
#   Daqi Yu <daqi.yu@astragenomics.com>                                                     #
#                                                                                           #
# Description:                                                                              #
#   Using the HuggingFace:AI4Protein/Prime_690M model to solve a regression problem.        #
#   using either full parameter fine tunning, LoRA or top-layers.                           #
# ------------------------------------------------------------------------------------------#

import torch
from transformers import AutoTokenizer, AutoModel, get_linear_schedule_with_warmup
from torch.utils.data import Dataset, DataLoader
import os, re
import pandas as pd
import numpy as np
import torch.nn as nn
import logging
from torch.nn.utils import clip_grad_norm_
from peft import LoraConfig, get_peft_model
from scipy import stats
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Fine-tuning a regression model using HuggingFace's Prime_690M model.")
    parser.add_argument('--seed', type=int, default=42, 
                       help='Random seed for reproducibility, default is 42.')
    parser.add_argument('--csv_path', type=str, required=True, 
                       help='Path to the CSV file containing the dataset.')
    parser.add_argument('--model_name', type=str, required=True, 
                       help='Path to the pre-trained model.')
    parser.add_argument('--batch_size', type=int, default=5, 
                       help='Batch size for training and evaluation.')
    parser.add_argument('--num_epochs', type=int, default=30, 
                       help='Number of epochs for training.')
    # Method-specific arguments
    parser.add_argument('--rank_count', type=int, default=4, 
                       help='Number of LoRA ranks to use (only for LoRA method), default is 4')
    parser.add_argument('--top_layers', type=int, default=4, 
                       help='Number of top layers to fine-tune (only for top_layer method), default is 4')
    parser.add_argument('--method', type=str, choices=['ft', 'lora', 'top_layer'], required=True,
                       help='Fine-tuning method: "ft", "lora" or "top_layer"')
    parser.add_argument('--learning_rate', type=float, default=0.0001,
                       help='Learning rate for optimizer, default is 0.0001')
    parser.add_argument('--warming_up', type=int, default=1000,
                       help='Warming up steps for optimizer, default is 1000')
    parser.add_argument('--use_listmle', action='store_true',
                        help='Use ListMLE as the loss function')
    parser.add_argument('--log_path', type=str, default=None, 
                       help='Path to the log file, default is ./default.log')

    args = parser.parse_args()

    # Validation
    if args.batch_size <= 0:
        raise ValueError("Batch size must be positive")
    if args.num_epochs <= 0:
        raise ValueError("Number of epochs must be positive")
    if args.learning_rate <= 0:
        raise ValueError("Learning rate must be positive")
    if not os.path.isfile(args.csv_path):
        raise FileNotFoundError(f"CSV file not found: {args.csv_path}")
    if not os.path.isdir(args.model_name):
        raise FileNotFoundError(f"Model directory not found: {args.model_name}")

    return args


def configure_logger(log_path=None):
    """Configure the logger for the script."""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    logger.handlers = []
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # File handler
    if log_path is None:
        log_path = os.path.join(os.getcwd(), 'default.log')
    try:
        file_handler = logging.FileHandler(log_path)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    except Exception as e:
        logger.error(f"Failed to create file handler: {e}")
        logger.warning("Continuing with console logging only")

    return logger

class CustomDataset(Dataset):
    def __init__(self, data, tokenizer, split="train"):
        self.data = data
        self.tokenizer = tokenizer
        if split == "train":
            self.input_ids, self.attention_mask, self.labels = self.load_dataset("Train")
        elif split == "valid":
            self.input_ids, self.attention_mask, self.labels = self.load_dataset("Valid")
        else:
            self.input_ids, self.attention_mask, self.labels = self.load_dataset("Test")

    def load_dataset(self, Class_ID):
        try:
            df = self.data.loc[self.data['Class'].isin([ Class_ID ])]
        except:
            raise ValueError( "Class ID error." )
        batch_tokens = self.tokenizer([" ".join(list(j)) for j in df['input']], return_tensors='pt', padding='longest')
        return (batch_tokens['input_ids'], 
                batch_tokens['attention_mask'], 
                torch.tensor([j for j in df['labels']]))

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return (self.input_ids[idx], 
                self.attention_mask[idx], 
                self.labels[idx])

class ESM2Regressor(nn.Module):
    def __init__(self, config):
        super(ESM2Regressor, self).__init__()
        self.method = config['method']
        
        if self.method == 'lora':
            # Initialize LoRA model
            base_model = AutoModel.from_pretrained(config['model_name'], trust_remote_code=True)
            lora_config = LoraConfig(
                r=config['rank_count'],
                lora_alpha=32,
                target_modules=["query", "value"],
                lora_dropout=0.1,
                bias="none"
            )
            self.esm2 = get_peft_model(base_model, lora_config)
        elif self.method == 'top_layer':
            # Initialize top-layer fine-tuning model
            self.esm2 = AutoModel.from_pretrained(config['model_name'], trust_remote_code=True)
            # Freeze all parameters
            for param in self.esm2.parameters():
                param.requires_grad = False
            # Unfreeze top layers
            for param in self.esm2.pro_prime.encoder.layer[-config['top_layers']:].parameters():
                param.requires_grad = True
            # Unfreeze sequence pooling
            for param in self.esm2.sequence_pooling.parameters():
                param.requires_grad = True
        else:
            self.esm2 = AutoModel.from_pretrained(config['model_name'], trust_remote_code=True)

        self.dropout = nn.Dropout(0.1)
        self.classifier = nn.Linear(self.esm2.config.hidden_size, 1)

    def forward(self, input_ids, attention_mask):
        outputs = self.esm2(input_ids=input_ids, attention_mask=attention_mask)
        pooled_output = outputs.sequence_hidden_states
        pooled_output = self.dropout(pooled_output)
        logits = self.classifier(pooled_output)
        return logits.squeeze(-1)

# ListMLE implementation.
class ListMLE(nn.Module):
    def __init__(self):
        super(ListMLE, self).__init__()

    def forward(self, predictions, targets):
        predictions = predictions.view(-1)
        targets = targets.view(-1)
        if len(predictions) == 0 or len(targets) == 0:
            return torch.tensor(0.0, device=predictions.device)
        if len(predictions) == 1:
            # For a single element, use MSE loss
            return (predictions - targets).pow(2).mean()

        #Only the rank in targets is considered.
        #Sort targets in descending order and get sorting indices
        _, sorted_indices = torch.sort(targets, descending=True)
        ordered_preds = predictions[sorted_indices]
        max_pred = ordered_preds.max()
        shifted_preds = ordered_preds - max_pred
        #ListMLE requires reverse cumulative summation
        reversed_shifted = torch.flip(shifted_preds, dims=(-1,))
        reversed_cumsum = torch.logcumsumexp(reversed_shifted, dim=-1)
        logcumsumexp = torch.flip(reversed_cumsum, dims=(-1,))
        log_likelihood = shifted_preds - logcumsumexp

        return -log_likelihood.mean()

def calculate_top_n_recall(y_true, y_pred, n=50):
    """
    Calculate the recall rate for top N predictions.
    
    :param y_true: True labels (ground truth values)
    :param y_pred: Predicted scores
    :param n: Number of top items to consider
    :return: Recall rate
    """
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    
    # Get the indices that would sort y_true and y_pred in descending order
    true_top_n = set(np.argsort(y_true)[::-1][:n])
    pred_top_n = set(np.argsort(y_pred)[::-1][:n])
    
    # Calculate the number of common items
    common_items = len(true_top_n.intersection(pred_top_n))
    
    # Calculate recall
    recall = common_items / n
    
    return recall

def train_and_evaluate(model, train_loader, val_loader, optimizer, scheduler, criterion, num_epochs, logger, device):
    for epoch in range(num_epochs):
        # Training phase
        model.train()
        train_loss = 0
        for i, (input_ids, attention_mask, target) in enumerate(train_loader):
            input_ids = input_ids.to(device)
            attention_mask = attention_mask.to(device)
            target = target.to(device)
            
            output = model(input_ids=input_ids, attention_mask=attention_mask)
            loss = criterion(output, target)
            
            optimizer.zero_grad()
            loss.backward()
            clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            scheduler.step()
            
            train_loss += loss.item()
            if i % 200 == 0:
                logger.info(f'Epoch [{epoch+1}/{num_epochs}] Step [{i+1}/{len(train_loader)}] Loss: {loss.item():.3f}')
        
        logger.info(f"Training loss: {train_loss/len(train_loader)}")

        # Validation phase
        model.eval()
        loader, phase = train_loader, 'Train'
        if loader is None:
            continue
        predictions = []
        targets = []
        
        with torch.no_grad():
            for input_ids, attention_mask, target in loader:
                input_ids = input_ids.to(device)
                attention_mask = attention_mask.to(device)
                target = target.to(device)
                
                output = model(input_ids=input_ids, attention_mask=attention_mask)
                
                predictions.extend(output.cpu().numpy())
                targets.extend(target.cpu().numpy())
            
            correlation = np.corrcoef(predictions, targets)[0, 1]
            spearmancor = stats.spearmanr( predictions, targets).statistic
            logger.info(f'{phase} Pearson correlation: {correlation:.3f}')
            logger.info(f'{phase} Spearman correlation: {spearmancor:.3f}')

        model.eval()
        loader, phase = val_loader, 'Validation'
        if loader is None:
            continue
        predictions = []
        targets = []
        
        with torch.no_grad():
            for input_ids, attention_mask, target in loader:
                input_ids = input_ids.to(device)
                attention_mask = attention_mask.to(device)
                target = target.to(device)
                
                output = model(input_ids=input_ids, attention_mask=attention_mask)
                
                predictions.extend(output.cpu().numpy())
                targets.extend(target.cpu().numpy())
            
            correlation = np.corrcoef(predictions, targets)[0, 1]
            spearmancor = stats.spearmanr( predictions, targets).statistic
            logger.info(f'{phase} Pearson correlation: {correlation:.3f}')
            logger.info(f'{phase} Spearman correlation: {spearmancor:.3f}')


def main():
    args = parse_args()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    logger = configure_logger(args.log_path)

    # Log configuration
    logger.info(f'Random seed: {args.seed}')
    logger.info(f'Model: {args.model_name}')
    logger.info(f'Data: {args.csv_path}')
    logger.info(f'Method: {args.method}')
    if args.method == 'lora':
        logger.info(f'Rank count: {args.rank_count}')
    elif args.method == 'top_layer':
        logger.info(f'Top layers: {args.top_layers}')
    logger.info(f'Batch size: {args.batch_size}')
    logger.info(f'Epochs: {args.num_epochs}')
    logger.info(f'lr: {args.learning_rate}')

    # Set random seed
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)

    # Load and preprocess single mutation data
    Ref_seq = "MITRSKNKNTMQAKVPPTDVETSTASLSSSPRSIFSSITPSSLSFKSRSLRSSTSSTISVDVVPEATVANNPTSATGTTNTRKKWSNEMNEFIWRNYLIITNLETDTNAYLNTLHTKFTEQFPNMTVSKQRIGDQRRAIVNNKLLPTITLDAIRSEVQTLLLTQNGTQQIQTQSQQEYTQINTLAHIQTKTRMKWSNELNETLMKLYFEVTKMETDMTMWRKNLRDKFIECYPHLAHLSEQRISDQKRVIINNKMIPETRLQIIKTVVEQEMHERDATVFQRSNSPILQQTTHYTHTQPQINVENTSLELGSQELDFTDENFDQQPTTINTHTENYERDSDIDELFLQALQYYTETHPTNRCYIPKQRTSKKFAKTVQYINKVILPDHVNTDTEFNQLQTIIYSAAWTAAKINGSKISLSFAQDRNNSQSNRKPRWQKRLEHKIEVIRSKIGRLTYYMQGHRSRKLDREVNKIIRDYKIHTIHEEPNRQLTHYLDTLKQKLTVIKSRLNRYTKCTLRKTQNSQFKNNEKQFYRTLKQSTTNHTTSNNTAHTAFPVEEDLRQFWAGIWENPVCHNTQADWIKTEEINYDTIQQMQFDQIDVETFHDVLNKAHNWKAPGTDHIHNFWLKKFTTIHPHLYNHINNFIKNPTNIPSYLTEGITYMLPKDLNDTANPAKYRPITCLQTMYKIITACITNIIYKHIDTHKILTEEQKGCRKFSQGCKEQLLIDSVILKDVQKYKKDLYTMYIDYKKAYDSVPHSWLLHILKIYKIHPYITRFLENIISKWQTKLHVNTSVGVVETEFIKIQRGIFQGDSLSPLWFCLALNPLSNQLKSLKAGFQLKYRNTVTTTEIRHQHCTINHLLYMDDIKLYAFSPEDLHQLADVTEIFSNDIGMEFGIDKCKINSVKAGKTYTHSYQLQTGDTIDSLDESGTYKYLGYIQTKQILHKEIKSQLTKQFKHRLKIICNTQLNARNTIKAINTYAIPILTYSFGIINWSKTELQNLQRNINTTLTKYRRHHPRSCTQRLTLPRNDGGRGIIDIPNLHNKQITTLRSYFYSKSPSSTLHNAIAQNDTKLTPLNMADQTTQNNETQTDRQTKITAWAHKSLHGRHYHDLNQSNVDKVASNEWLRRGELFPETEGFMLAIQDQVIETRNYQKHIMHIPIPTDICRRCNSASETIQHVTGACKAIVQSDYKHRHDQVANIIHQKLALKYKLISGSPTPYYKYHPEIVLENSAFKLYFDRAILTDKTTHFNRPDITLIDKINKSAQIIDIAIPNTHNLQNTVAEKLSKYTDLKIEITRMWRLNSVTIVPIVLSTTGVIPKQLHQSIKTLELPPYTYRTLQKAAILNTCHIVRKFLNETDKVSFTPTQN"
    def mutaion_to_sequence(sequence, mutation):
        _seq = sequence
        for _mutation in re.split( r';|_|,', mutation ):
            _AA  = re.sub( "[0-9].*", "", _mutation )
            _pos = int( re.sub( "[A-Z]", "", _mutation ) )
            _mut = re.sub( "^.*[0-9]", "", _mutation )
            _seq = _seq[:_pos-1] + _mut + _seq[_pos:]
        return _seq
    
    data = pd.read_csv(args.csv_path, index_col = 0)

    tokenizer = AutoTokenizer.from_pretrained(args.model_name, trust_remote_code=True)

    #Note to remove the N terminal part of protein due to PLLM input length limitation
    data['input']  = data.apply( lambda x:x['input'][335:1355], axis = 1)

    # Create datasets and dataloaders
    train_dataset = CustomDataset(data, tokenizer, "train")
    val_dataset   = CustomDataset(data, tokenizer, "valid")

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
    val_loader   = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False)

    logger.info(f'Train size: { len(train_dataset) }')
    logger.info(f'Valid size: { len(val_dataset) }')
    
    # Initialize model
    config = vars(args)
    model = ESM2Regressor(config).to(device)

    # Training setup
    if args.use_listmle:
        criterion = ListMLE()
        logger.info('ListMLE: Ture')
    else:
        logger.info('MSELoss: Ture')
        criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate, eps=1e-6)
    total_steps = len(train_loader) * args.num_epochs
    if args.warming_up < total_steps:
        num_warmup_steps = args.warming_up
    else:
        num_warmup_steps = int(0.1 * total_steps)  # 10% of total training steps as warmup steps
        logger.warning(f"Provided warm-up steps ({args.warming_up}) exceed total steps ({total_steps}). Using 10% of total steps.")
    logger.info(f'Warming up steps: {num_warmup_steps}')
    scheduler = get_linear_schedule_with_warmup(optimizer, 
                                              num_warmup_steps=num_warmup_steps, 
                                              num_training_steps=total_steps)

    # Train and evaluate
    train_and_evaluate(model, train_loader, val_loader, 
                                    optimizer, scheduler, criterion, args.num_epochs, 
                                    logger, device)

if __name__ == "__main__":
    main()
